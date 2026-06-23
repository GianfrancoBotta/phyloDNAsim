"""
Command-line entry point for phylodnasim.
"""

import argparse
import gzip
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from Bio import SeqIO

from .simulation_utils import (
    getTree,
    getPaths_and_TimeMatrix,
    generateOrder,
    saveMutations,
    getDirichletClone,
    wgsSim_bulk_parallel,
    wgsSim_sc_parallel,
    targetedSim_bulk_parallel,
    targetedSim_sc_parallel,
    aggregate_fastqs,
)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="phylodnasim",
        description="Phylogenetic DNA tumour biopsy simulator",
    )

    # ── I/O ────────────────────────────────────────────────────────────────
    io = parser.add_argument_group("I/O")
    io.add_argument("--outdir",    required=True,  help="Output directory")
    io.add_argument("--genome",    required=True,  help="Reference genome (FASTA)")
    io.add_argument("--bed",       default=None,   help="BED file for targeted sequencing")
    io.add_argument("--mut_infos", default=None,   help="Pre-computed mutation info YAML")
    io.add_argument("--clone_prop",default=None,   help="Pre-computed clone proportions file")
    io.add_argument("--prefix",    default=None,   help="Output filename prefix")
    io.add_argument("--threads",   type=int, default=1)

    # ── Simulation mode ────────────────────────────────────────────────────
    mode = parser.add_argument_group("Simulation mode")
    mode.add_argument("--generate_mutations", action="store_true",
                      help="Only generate & save mutations, then exit")
    mode.add_argument("--paired",   action="store_true", help="Paired-end reads")
    mode.add_argument("--targeted", action="store_true", help="Targeted (panel) sequencing")
    mode.add_argument("--bulk",     action="store_true", help="Bulk sequencing (default: single-cell)")

    # ── Tree / clones ──────────────────────────────────────────────────────
    tree = parser.add_argument_group("Tree / clones")
    tree.add_argument("--num_leaves",       type=int,   required=True)
    tree.add_argument("--num_single_cells", type=int,   default=None,
                      help="Number of single cells to simulate")
    tree.add_argument("--prop_hc",          type=float, default=0.1,
                      help="Proportion of healthy cells")
    tree.add_argument("--alpha",            type=int,   default=10,
                      help="Dirichlet concentration for clone proportions")
    tree.add_argument("--population_size",  type=float, default=8.0e8)

    # ── Sequencing ─────────────────────────────────────────────────────────
    seq = parser.add_argument_group("Sequencing")
    seq.add_argument("--read_len",     type=int,   default=150)
    seq.add_argument("--frag_len",     type=int,   default=400)
    seq.add_argument("--error_rate",   type=float, default=0.001)
    seq.add_argument("--coverage",     type=float, default=None,
                     help="Total coverage (bulk or aggregate single-cell)")
    seq.add_argument("--cell_coverage",type=float, default=None,
                     help="Per-cell coverage (single-cell mode)")
    seq.add_argument("--r",  type=float, default=None,
                     help="NB dispersion r (required in single-cell mode)")
    seq.add_argument("--p",  type=float, default=None,
                     help="NB success prob p (required in single-cell mode)")

    return parser


def validate(args, parser):
    if args.coverage is None and args.cell_coverage is None and not args.generate_mutations:
        parser.error("Either --coverage or --cell_coverage is required.")

    if args.targeted and args.bed is None:
        parser.error("--bed is required when --targeted is set.")

    if not args.bulk and not args.generate_mutations:
        missing = [k for k in ["num_single_cells", "r", "p"]
                   if getattr(args, k) is None]
        if missing:
            parser.error(f"Single-cell mode requires: {', '.join('--' + m for m in missing)}")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    validate(args, parser)

    # ── Convenience aliases ────────────────────────────────────────────────
    working_dir = args.outdir
    full_genome = args.genome
    bed = args.bed
    infos_path = args.mut_infos
    clone_prop = args.clone_prop
    threads = args.threads

    generate_mutations = args.generate_mutations
    paired = args.paired
    targeted = args.targeted
    bulk = args.bulk

    num_clones = args.num_leaves
    num_single_cells = args.num_single_cells

    read_len = args.read_len
    frag_len = args.frag_len
    prop_hc = args.prop_hc
    alpha = args.alpha
    r = args.r
    p = args.p

    error_rate = args.error_rate
    coverage = args.coverage
    cell_coverage = args.cell_coverage
    pop = args.population_size
    prefix = args.prefix

    # ── Load genome ────────────────────────────────────────────────────────
    opener = gzip.open if full_genome.endswith(".gz") else open

    with opener(full_genome, "rt") as handle:
        chrom_names, chroms = map(list, zip(*((rec.id, bytearray(str(rec.seq).upper(), "utf-8")) for rec in SeqIO.parse(handle, "fasta"))))
    numchrommap = dict(enumerate(chrom_names))
    rev_numchrommap = {v: k for k, v in numchrommap.items()}
    tab = bytes.maketrans(b"ACTG", b"TGAC")

    # ── Load BED if targeted ───────────────────────────────────────────────
    regions = None
    if targeted:
        with open(bed) as f:
            lines = [line.strip().split() for line in f]
        lines.sort(key=lambda x: (rev_numchrommap[x[0]], int(x[1]), int(x[2])))
        regions = pd.DataFrame(lines, columns=["chr", "start", "end"])
        regions["start"] = pd.to_numeric(regions["start"], downcast="integer", errors="coerce")
        regions["end"]   = pd.to_numeric(regions["end"],   downcast="integer", errors="coerce")

    # ── Generate or load mutations ─────────────────────────────────────────
    if infos_path is None:
        # mut_events = ["SNV", "CNV", "DEL", "DELSMALL", "INVERSION", "TRANSLOCATION",
        #               "BFB", "CHROMOTHRIP", "INSERTIONSMALL", "KATAEGIS", "ANEUPLOIDY"]

        ultralow = [7.0e-10, 5.0e-10, 3.0e-10]
        low = [5.0e-9,  2.0e-9,  9.0e-10]
        medium = [1.0e-8,  9.0e-9,  7.0e-9]
        high = [7.0e-8,  5.0e-8,  3.0e-8]

        list_of_rates = {
            "SNV": high,
            "CNV": high,
            "DEL": medium,
            "DELSMALL": medium,
            "INVERSION": low,
            "TRANSLOCATION": low,
            "BFB": ultralow,
            "CHROMOTHRIP": ultralow,
            "INSERTIONSMALL": medium,
            "KATAEGIS": ultralow,
            "ANEUPLOIDY": ultralow,
        }
        if targeted:
            list_of_rates = {k: [x / 10 for x in v] for k, v in list_of_rates.items()}

        tot_nodes = 2 * num_clones - 1
        os.makedirs(working_dir, exist_ok=True)

        tree = getTree(num_clones, pop, working_dir)
        list_of_paths, time_matrix, _ = getPaths_and_TimeMatrix(tree, num_clones)
        mutationedge_list, _ = generateOrder(tree, time_matrix, list_of_rates)

        if targeted:
            infos = saveMutations(chroms, tot_nodes, list_of_paths,
                                  mutationedge_list, numchrommap, targeted, regions)
        else:
            infos = saveMutations(chroms, tot_nodes, list_of_paths,
                                  mutationedge_list, numchrommap)

        fname = (prefix + "_information_list.yaml") if prefix else "information_list.yaml"
        with open(os.path.join(working_dir, fname), "w") as f:
            json.dump(infos, f, indent=4)
    else:
        with open(infos_path) as f:
            infos = {int(k): v for k, v in json.load(f).items()}

    if generate_mutations:
        sys.exit(0)

    # ── Clone proportions ──────────────────────────────────────────────────
    tot_nodes = 2 * num_clones - 1
    int_nodes = tot_nodes - 1

    os.makedirs(working_dir, exist_ok=True)
    print("Starting simulation.")
    ts = time.time()

    if clone_prop is None:
        clone_prop_arr = getDirichletClone(int_nodes, alpha)
        fname = (prefix + "_clones_proportions.txt") if prefix else "clones_proportions.txt"
        np.savetxt(os.path.join(working_dir, fname), clone_prop_arr)
    else:
        clone_prop_arr = np.loadtxt(clone_prop)

    # ── Common kwargs for simulation functions ─────────────────────────────
    common = dict(
        ls=chroms, rl=read_len, fl=frag_len,
        floc=working_dir, erate=error_rate, tab=tab,
        infos=infos, paired=paired,
    )

    # ── Write parameter file ───────────────────────────────────────────────
    pfname = (prefix + "_parameter_list.yaml") if prefix else "parameter_list.yaml"
    with open(os.path.join(working_dir, pfname), "w") as f:
        f.write(f"num leaves: {num_clones}\n")
        f.write(f"num internal nodes: {int_nodes}\n")
        f.write(f"dir_conc: {alpha}\n")
        f.write(f"cell pop: {pop}\n")
        f.write(f"coverage: {coverage}\n")
        f.write(f"read len: {read_len}\n")
        f.write(f"frag len: {frag_len}\n")
        f.write(f"paired: {paired}\n")
        f.write(f"targeted: {targeted}\n")
        f.write(f"error rate: {error_rate}\n")
        f.write(f"proportions of healthy cells: {prop_hc}\n")
        if not bulk:
            f.write(f"num single cells: {num_single_cells}\n")
            f.write(f"NB parameters: r={r} p={p}\n")

    # ── Run simulation ─────────────────────────────────────────────────────
    if bulk:
        current_cov = coverage
        if targeted:
            targetedSim_bulk_parallel(
                prop_hc=prop_hc, coverage=current_cov,
                num_clones=int_nodes, raw_clone_prop=clone_prop_arr,
                threads=threads, regions=regions,
                rev_numchrommap=rev_numchrommap, **common,
            )
        else:
            wgsSim_bulk_parallel(
                prop_hc=prop_hc, coverage=current_cov,
                num_clones=int_nodes, raw_clone_prop=clone_prop_arr,
                threads=threads, **common,
            )

        lname = (prefix + "_bulkleft.fq.gz")  if prefix else "bulkleft.fq.gz"
        rname = (prefix + "_bulkright.fq.gz") if prefix else "bulkright.fq.gz"

    else:
        current_cov = cell_coverage * num_single_cells if coverage is None else coverage
        if targeted:
            targetedSim_sc_parallel(
                num_single_cells=num_single_cells,
                prop_hc=prop_hc, coverage=current_cov,
                num_clones=int_nodes, clone_prop=clone_prop_arr,
                r=r, p=p, threads=threads,
                regions=regions, rev_numchrommap=rev_numchrommap,
                **common,
            )
        else:
            wgsSim_sc_parallel(
                num_single_cells=num_single_cells,
                prop_hc=prop_hc, coverage=current_cov,
                num_clones=int_nodes, clone_prop=clone_prop_arr,
                r=r, p=p, threads=threads, **common,
            )

        lname = (prefix + "_scleft.fq.gz")  if prefix else "scleft.fq.gz"
        rname = (prefix + "_scright.fq.gz") if prefix else "scright.fq.gz"

    aggregate_fastqs(
        working_dir,
        os.path.join(working_dir, lname),
        os.path.join(working_dir, rname),
        paired,
    )

    te = time.time()
    mode_label = "bulk" if bulk else "single-cell"
    print(f"Time elapsed for {mode_label} biopsy simulation: {te - ts:.1f}s")


if __name__ == "__main__":
    main()
