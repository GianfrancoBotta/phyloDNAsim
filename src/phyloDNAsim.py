import argparse
import glob
import gzip
import math
import msprime
import numpy as np
import os
import pandas as pd
import pickle
import psutil
import random
import re
import string
import sys
import time
import tskit

from Bio import SeqIO
from datasketch import *
from simulation_utils import *

### Parse arguments
parser = argparse.ArgumentParser(description="Simulation parameters")

# Inputs / outputs
parser.add_argument("--outdir", required=True)
parser.add_argument("--genome", required=True)
parser.add_argument("--bed", default=None)
parser.add_argument("--threads", type=int, default=1)
# parser.add_argument("--log", default=None)

# Tumor biopsy simulation parameters
parser.add_argument("--paired", action="store_true")
parser.add_argument("--targeted", action="store_true")
parser.add_argument("--bulk", action="store_true")

parser.add_argument("--num_samples", type=int, required=True)
parser.add_argument("--num_leaves", type=int, required=True)
parser.add_argument("--num_single_cells_list", type=int, nargs='+', default=None, help="List of single cell counts")

parser.add_argument("--read_len_tb", type=int, required=True)
parser.add_argument("--frag_len_tb", type=int, required=True)

parser.add_argument("--prop_hc_tb", type=float, required=True)
parser.add_argument("--alpha", type=int, required=True)
parser.add_argument("--r", type=float, default=None, help="Parameter r (required if not bulk)")
parser.add_argument("--p", type=float, default=None, help="Parameter p (required if not bulk)")

parser.add_argument("--error_rate", type=float, default=0.001)
parser.add_argument("--coverage_tb", type=float, default=None, help="Total coverage")
parser.add_argument("--cell_coverage_tb", type=float, default=None, help="Per cell coverage")
parser.add_argument("--genome_length", type=float, default=3095693983.0)
parser.add_argument("--population_size", type=float, default=8.0e+8)


# Liquid biopsies simulation parameters
parser.add_argument("--tumor_liquid_biopsy", action="store_true")
parser.add_argument("--read_len_lb", type=int, default=None)
parser.add_argument("--frag_len_lb", type=int, default=None)
parser.add_argument("--coverage_lb", type=float, default=None)
parser.add_argument("--prop_hc_tlb", type=float, default=None)

parser.add_argument("--plasma_liquid_biopsy", action="store_true")
parser.add_argument("--prop_hc_plb", type=float, default=None)

args = parser.parse_args()

# Check coverage
if args.coverage_tb is None and args.cell_coverage_tb is None:
    parser.error(f"either --coverage_tb or --cell_coverage_tb (in single-cell mode) is required")

# Check targeted sequencing parameters
if args.targeted:
    if getattr(args, "bed") is None:
        parser.error(f"--bed is required when --bulk is not set")
    
# Check single-cell parameters
if not args.bulk:
    for k in ["num_single_cells_list", "r", "p", "prop_hc_tlb"]:
        if getattr(args, k) is None:
            parser.error(f"--{k} is required when --bulk is not set")

# Check tumor liquid biopsy parameters
if args.tumor_liquid_biopsy:
    for k in ["read_len_lb", "frag_len_lb", "coverage_lb", "prop_hc_tlb"]:
        if getattr(args, k) is None:
            parser.error(f"--{k} is required when --tumor_liquid_biopsy is set")

# Check plasma liquid biopsy parameters
if args.plasma_liquid_biopsy:
    for k in ["read_len_lb", "frag_len_lb", "coverage_lb", "prop_hc_plb"]:
        if getattr(args, k) is None:
            parser.error(f"--{k} is required when --plasma_liquid_biopsy is set")

# # Redirect stdout to log
# if args.log is not None:
#     log = open(args.log, "a")
#     sys.stdout = log

# Assign variables for tumor biopsy data
working_dir = args.outdir
full_genome = args.genome
bed = args.bed
threads = args.threads

paired = args.paired
targeted = args.targeted
bulk = args.bulk

num_samples = args.num_samples
num_clones = args.num_leaves
num_single_cells_list = args.num_single_cells_list

read_len_tb = args.read_len_tb
frag_len_tb = args.frag_len_tb

prop_hc_tb = args.prop_hc_tb
alpha = args.alpha
r = args.r
p = args.p

error_rate = args.error_rate
coverage_tb = args.coverage_tb
cell_coverage_tb = args.cell_coverage_tb
genome_length = args.genome_length
pop = args.population_size

# Load genome and bed files
chrom_names, chroms = map(list, zip(*((record.id, bytearray(str(record.seq).upper(), 'utf-8')) for record in SeqIO.parse(full_genome, "fasta")))) # Parse sequences and names together
numchrommap = dict(zip(range(len(chrom_names)), chrom_names))
rev_numchrommap = {v: k for k, v in numchrommap.items()}
tab = bytes.maketrans(b"ACTG", b"TGAC")

if targeted:
    EXON_FILE = bed
    total_num_intervals = 0
    exonDict = {}
    panel_chroms = [bytearray() for _ in chroms]
    strings_to_idx = []
    for i in numchrommap.values():
        exonDict[i] = []
    with open(EXON_FILE, 'r') as f:
        lines = [line.strip().split() for line in f]
    lines.sort(key=lambda x: (rev_numchrommap[x[0]], int(x[1]), int(x[2])))
    regions = pd.DataFrame(lines, columns = ['chr', 'start', 'end'])
    regions['start'] = pd.to_numeric(regions['start'], downcast='integer', errors='coerce')
    regions['end'] = pd.to_numeric(regions['end'], downcast='integer', errors='coerce')
    
# Mutation process rate lists (dependent on the length of the regions we are analysing)
mut_events =  ["SNV", "CNV", "DEL", "DELSMALL", "INVERSION", "TRANSLOCATION", "BFB", "CHROMOTHRIP", "INSERTIONSMALL", "KATAEGIS", "ANEUPLOIDY"]
# Set rates for mutational events
ultralow_rates_list = [7.0e-10, 5.0e-10, 3.0e-10]
low_rates_list = [5.0e-9, 2.0e-9, 9.0e-10]
medium_rates_list = [1.0e-8, 9.0e-9, 7.0e-9]
high_rates_list = [7.0e-8, 5.0e-8, 3.0e-8]
list_of_rates = {
    "SNV": high_rates_list,
    "CNV": high_rates_list,
    "DEL": medium_rates_list,
    "DELSMALL": medium_rates_list,
    "INVERSION": low_rates_list,
    "TRANSLOCATION": low_rates_list,
    "BFB": ultralow_rates_list,
    "CHROMOTHRIP": ultralow_rates_list,
    # "CHROMOPLEX": ultralow_rates_list,
    "INSERTIONSMALL": medium_rates_list,
    "KATAEGIS": ultralow_rates_list,
    "ANEUPLOIDY": ultralow_rates_list
}
if targeted:
    list_of_rates = {k: [float(r) / 10 for r in v] for k, v in list_of_rates.items()}

# Create mutations and store them
tot_nodes = 2*num_clones - 1
int_nodes = tot_nodes - 1
os.makedirs(working_dir, exist_ok=True)
tree = getTree(num_clones, pop, working_dir)
list_of_paths, time_matrix, depth  = getPaths_and_TimeMatrix(tree, num_clones)
mutationedge_list, avg_rate_list = generateOrder(tree, time_matrix, list_of_rates)
infos = saveMutations(
    chroms,
    tot_nodes,
    list_of_paths,
    mutationedge_list,
    numchrommap,
    targeted,
    regions)
with open(os.path.join(working_dir, 'information_list.json'), 'w') as f:
    json.dump(infos, f, indent=4)
    
# Compute clone proportions
clone_prop = getDirichletClone(int_nodes, alpha)


print('Start simulation for biopsy data.')
ts_sb = time.time()
for num_single_cells in num_single_cells_list:
    # Compute coverage if coverage per cell is passed
    if coverage_tb is None:
        current_cov = cell_coverage_tb * num_single_cells
    else:
        current_cov = coverage_tb
    for sample in range(num_samples):
        sample_working_dir = os.path.join(working_dir, f'sample_{sample+1}', f'cell_{num_single_cells}')
        os.makedirs(sample_working_dir, exist_ok=True)
        with open(os.path.join(sample_working_dir, 'parameter_list_tb.yaml'), 'w') as f:
            f.write('num leaves: ' + str(num_clones)+'\n')
            f.write('num internal nodes: ' + str(int_nodes)+'\n')
            f.write('dir_conc: ' + str(alpha)+'\n')
            f.write('cell pop: ' + str(pop)+'\n')
            f.write('coverage: ' + str(coverage_tb)+'\n')
            f.write('num single cells: ' + str(num_single_cells)+'\n')
            f.write('read len: ' + str(read_len_tb)+'\n')
            f.write('frag len: ' + str(frag_len_tb)+'\n')
            f.write('paired: ' + str(paired)+'\n')
            f.write('targeted: ' + str(targeted)+'\n')
            f.write('rates of variants: ' + str(avg_rate_list)+'\n')
            f.write('full poisson time: ' + str(depth) + '\n')
            f.write('error rate: ' + str(error_rate) + '\n')
            f.write('proportions of healthy cells: ' + str(prop_hc_tb) + '\n')
            if(not bulk):
                f.write('NB parameters: r=' + str(r) + ' p=' + str(p) + '\n')

        if targeted:
            if bulk: # Bulk simulation
                targetedSim_bulk_parallel(prop_hc = prop_hc_tb,
                                                       coverage = current_cov,
                                                       num_clones = int_nodes,
                                                       clone_prop = clone_prop,
                                                       threads = threads,
                                                       ls = chroms,
                                                       rl = read_len_tb,
                                                       fl = frag_len_tb,
                                                       floc = sample_working_dir,
                                                       regions = regions,
                                                       rev_numchrommap = rev_numchrommap,
                                                       erate = error_rate,
                                                       tab = tab,
                                                       infos = infos,
                                                       paired = paired)
                
            else:
                targetedSim_sc_parallel(num_single_cells = num_single_cells,
                                                     prop_hc = prop_hc_tb,
                                                     coverage = current_cov,
                                                     num_clones = int_nodes,
                                                     clone_prop = clone_prop,
                                                     r = r,
                                                     p = p,
                                                     threads = threads,
                                                     ls = chroms,
                                                     rl = read_len_tb,
                                                     fl = frag_len_tb,
                                                     floc = sample_working_dir,
                                                     regions = regions,
                                                     rev_numchrommap = rev_numchrommap,
                                                     erate = error_rate,
                                                     tab = tab,
                                                     infos = infos,
                                                     paired = paired)
                
        else:            
            if bulk: # Bulk simulation
                wgsSim(ls = chroms,
                    num_clones = int_nodes,
                    coverage = current_cov,
                    rl = read_len_tb,
                    fl = frag_len_tb,
                    floc = sample_working_dir,
                    alpha = alpha,
                    erate = error_rate,
                    tab = tab,
                    infos = infos,
                    num_single_cells=1,
                    flag=0,
                    paired = paired)
            else:
                wgsSim(ls = chroms,
                    num_clones = int_nodes,
                    coverage = current_cov,
                    rl = read_len_tb,
                    fl = frag_len_tb,
                    floc = sample_working_dir,
                    alpha = alpha,
                    erate = error_rate,
                    tab = tab,
                    infos = infos,
                    r = r,
                    p = p,
                    num_single_cells=num_single_cells,
                    flag=0,
                    paired = paired)
                
        # Aggregate all fastq.gz files from different cells/threads
        if(bulk):
            aggregate_fastqs(sample_working_dir, os.path.join(sample_working_dir, "bulkleft.fq.gz"), os.path.join(sample_working_dir, "bulkright.fq.gz"), paired)
        else:
            aggregate_fastqs(sample_working_dir, os.path.join(sample_working_dir, "scleft.fq.gz"), os.path.join(sample_working_dir, "scright.fq.gz"), paired)

te_sb = time.time()
if(bulk):
    print('Time elapsed for bulk biopsy simulation', te_sb-ts_sb)
else:
    print('Time elapsed for single-cell biopsy simulation', te_sb-ts_sb)
print('Finished simulation for biopsy data.')



# Assign variables for tumor liquid biopsy data
tumor_liquid_biopsy = args.tumor_liquid_biopsy
read_len_lb = args.read_len_lb
frag_len_lb = args.frag_len_lb
coverage_lb = args.coverage_lb
prop_hc_tlb = args.prop_hc_tlb

if(tumor_liquid_biopsy):
    print('Start simulation for tumor liquid biopsy data.')
    ts_tlb = time.time()
    for num_single_cells in num_single_cells_list:
        # Compute coverage if coverage per cell is passed
        if coverage_tb is None:
            coverage_tb = cell_coverage_tb * num_single_cells
        for sample in range(num_samples):
                sample_working_dir = os.path.join(working_dir, f'sample_{sample+1}', f'cell_{num_single_cells}')
                os.makedirs(sample_working_dir, exist_ok=True)
                with open(os.path.join(sample_working_dir, 'parameter_list_tlb.yaml'), 'w') as f:
                    f.write('num leaves: ' + str(num_clones)+'\n')
                    f.write('num internal nodes: ' + str(int_nodes)+'\n')
                    f.write('dir_conc: ' + str(alpha)+'\n')
                    f.write('cell pop: ' + str(pop)+'\n')
                    f.write('coverage: ' + str(coverage_lb)+'\n')
                    f.write('num single cells: ' + str(num_single_cells)+'\n')
                    f.write('read len: ' + str(read_len_lb)+'\n')
                    f.write('frag len: ' + str(frag_len_lb)+'\n')
                    f.write('paired: ' + str(paired)+'\n')
                    f.write('targeted: ' + str(targeted)+'\n')
                    f.write('rates of variants: ' + str(avg_rate_list)+'\n')
                    f.write('full poisson time: ' + str(depth) + '\n')
                    f.write('error rate: ' + str(error_rate) + '\n')
                    f.write('proportions of healthy cells: ' + str(prop_hc_tlb) + '\n')

                if targeted:
                    targetedSim_bulk_parallel(prop_hc = prop_hc_tlb,
                                                            coverage = coverage_lb,
                                                            num_clones = int_nodes,
                                                            alpha = alpha,
                                                            clone_prop = clone_prop,
                                                            threads = threads,
                                                            ls = chroms,
                                                            rl = read_len_lb,
                                                            fl = frag_len_lb,
                                                            floc = sample_working_dir,
                                                            regions = regions,
                                                            rev_numchrommap = rev_numchrommap,
                                                            erate = error_rate,
                                                            tab = tab,
                                                            infos = infos,
                                                            paired = paired)
                
                aggregate_fastqs(sample_working_dir, os.path.join(sample_working_dir, "T.ctleft.fq.gz"), os.path.join(sample_working_dir, "T.ctright.fq.gz"), paired)
                
    te_tlb = time.time()
    print('Time elapsed for tumor liquid biopsy simulation', te_tlb-ts_tlb)



# Assign variables for plasma liquid biopsy data
plasma_liquid_biopsy = args.plasma_liquid_biopsy
prop_hc_plb = args.prop_hc_plb

if(plasma_liquid_biopsy):
    print('Start simulation for plasma liquid biopsy data.')
    ts_plb = time.time()
    for num_single_cells in num_single_cells_list:
        # Compute coverage if coverage per cell is passed
        if coverage_tb is None:
            coverage_tb = cell_coverage_tb * num_single_cells
        for sample in range(num_samples):
            sample_working_dir = os.path.join(working_dir, f'sample_{sample+1}', f'cell_{num_single_cells}')
            os.makedirs(sample_working_dir, exist_ok=True)
            with open(os.path.join(sample_working_dir, 'parameter_list_plb.yaml'), 'w') as f:
                f.write('num leaves: ' + str(num_clones)+'\n')
                f.write('num internal nodes: ' + str(int_nodes)+'\n')
                f.write('dir_conc: ' + str(alpha)+'\n')
                f.write('cell pop: ' + str(pop)+'\n')
                f.write('coverage: ' + str(coverage_lb)+'\n')
                f.write('num single cells: ' + str(num_single_cells)+'\n')
                f.write('read len: ' + str(read_len_lb)+'\n')
                f.write('frag len: ' + str(frag_len_lb)+'\n')
                f.write('paired: ' + str(paired)+'\n')
                f.write('targeted: ' + str(targeted)+'\n')
                f.write('rates of variants: ' + str(avg_rate_list)+'\n')
                f.write('full poisson time: ' + str(depth) + '\n')
                f.write('error rate: ' + str(error_rate) + '\n')
                f.write('proportions of healthy cells: ' + str(prop_hc_plb) + '\n')

            if targeted:
                targetedSim_bulk_parallel(prop_hc = prop_hc_plb,
                                                        coverage = coverage_lb,
                                                        num_clones = int_nodes,
                                                        alpha = alpha,
                                                        clone_prop = clone_prop,
                                                        threads = threads,
                                                        ls = chroms,
                                                        rl = read_len_lb,
                                                        fl = frag_len_lb,
                                                        floc = sample_working_dir,
                                                        regions = regions,
                                                        rev_numchrommap = rev_numchrommap,
                                                        erate = error_rate,
                                                        tab = tab,
                                                        infos = infos,
                                                        paired = paired)
            
            aggregate_fastqs(sample_working_dir, os.path.join(sample_working_dir, "P.ctleft.fq.gz"), os.path.join(sample_working_dir, "P.ctright.fq.gz"), paired)

    te_plb = time.time()
    print('Time elapsed for plasma liquid biopsy simulation', te_plb-ts_plb)
