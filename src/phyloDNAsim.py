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
parser.add_argument("--mut_infos", default=None)
parser.add_argument("--clone_prop", default=None)
parser.add_argument("--threads", type=int, default=1)
# parser.add_argument("--log", default=None)

# Tumor biopsy simulation parameters
parser.add_argument("--generate_mutations", action="store_true")
parser.add_argument("--paired", action="store_true")
parser.add_argument("--targeted", action="store_true")
parser.add_argument("--bulk", action="store_true")

parser.add_argument("--num_leaves", type=int, required=True)
parser.add_argument("--num_single_cells", type=int, default=None, help="Cells to simulate")

parser.add_argument("--read_len", type=int, default=150)
parser.add_argument("--frag_len", type=int, default=400)

parser.add_argument("--prop_hc", type=float, default=0.1)
parser.add_argument("--alpha", type=int, default=10)
parser.add_argument("--r", type=float, default=None, help="Parameter r (required if not bulk)")
parser.add_argument("--p", type=float, default=None, help="Parameter p (required if not bulk)")

parser.add_argument("--error_rate", type=float, default=0.001)
parser.add_argument("--coverage", type=float, default=None, help="Total coverage")
parser.add_argument("--cell_coverage", type=float, default=None, help="Per cell coverage")
parser.add_argument("--population_size", type=float, default=8.0e+8)

parser.add_argument("--prefix", default=None)

args = parser.parse_args()
            
# Check coverage
if args.coverage is None and args.cell_coverage is None and not args.generate_mutations:
    parser.error(f"either --coverage or --cell_coverage (in single-cell mode) is required")

# Check targeted sequencing parameters
if args.targeted:
    if getattr(args, "bed") is None:
        parser.error(f"--bed is required when --bulk is not set")
    
# Check single-cell parameters
if not args.bulk and not args.generate_mutations:
    for k in ["num_single_cells_list", "r", "p", "prop_hc"]:
        if getattr(args, k) is None:
            parser.error(f"--{k} is required when --bulk is not set")

# Assign variables for tumor biopsy data
working_dir = args.outdir
full_genome = args.genome
bed = args.bed
infos = args.mut_infos
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

prop_hc_tb = args.prop_hc
alpha = args.alpha
r = args.r
p = args.p

error_rate = args.error_rate
coverage = args.coverage
cell_coverage = args.cell_coverage
pop = args.population_size

prefix = args.prefix

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
    
if infos is None:
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
    if targeted:
        infos = saveMutations(
            chroms,
            tot_nodes,
            list_of_paths,
            mutationedge_list,
            numchrommap,
            targeted,
            regions)
    else:
        infos = saveMutations(
            chroms,
            tot_nodes,
            list_of_paths,
            mutationedge_list,
            numchrommap)
        
    if prefix is not None:
            info_fname = prefix + '_information_list.yaml'
    else:
        info_fname = 'information_list.yaml'
            
    with open(os.path.join(working_dir, info_fname), 'w') as f:
        json.dump(infos, f, indent=4)
else:
    # Load infos
    with open(infos, 'r') as f:
        infos = {int(k): v for k, v in json.load(f).items()}

# Exit if the goal was to generate the mutations
if generate_mutations:
    sys.exit(0)

tot_nodes = 2*num_clones - 1
int_nodes = tot_nodes - 1

print('Start simulation.')
ts = time.time()
os.makedirs(working_dir, exist_ok=True)

if clone_prop is None:
    # Compute clone proportions
    clone_prop = getDirichletClone(int_nodes, alpha)
    
    if prefix is not None:
        props_fname = prefix + '_clones_proportions.txt'
    else:
        props_fname = 'clones_proportions.txt'
            
    np.savetxt(os.path.join(working_dir, props_fname), clone_prop)
    
else:
    # Load clone proportions
    clone_prop = np.loadtxt(clone_prop)

##################### BULK SIMULATION
if bulk: # Bulk simulation
    current_cov = coverage
    
    if prefix is not None:
        params_fname = prefix + '_parameter_list.yaml'
    else:
        params_fname = 'parameter_list.yaml'
        
    with open(os.path.join(working_dir, params_fname), 'w') as f:
        f.write('num leaves: ' + str(num_clones)+'\n')
        f.write('num internal nodes: ' + str(int_nodes)+'\n')
        f.write('dir_conc: ' + str(alpha)+'\n')
        f.write('cell pop: ' + str(pop)+'\n')
        f.write('coverage: ' + str(coverage)+'\n')
        f.write('read len: ' + str(read_len)+'\n')
        f.write('frag len: ' + str(frag_len)+'\n')
        f.write('paired: ' + str(paired)+'\n')
        f.write('targeted: ' + str(targeted)+'\n')
        f.write('error rate: ' + str(error_rate) + '\n')
        f.write('proportions of healthy cells: ' + str(prop_hc_tb) + '\n')
    
    if targeted:
        targetedSim_bulk_parallel(prop_hc = prop_hc_tb,
                                    coverage = current_cov,
                                    num_clones = int_nodes,
                                    raw_clone_prop = clone_prop,
                                    threads = threads,
                                    ls = chroms,
                                    rl = read_len,
                                    fl = frag_len,
                                    floc = working_dir,
                                    regions = regions,
                                    rev_numchrommap = rev_numchrommap,
                                    erate = error_rate,
                                    tab = tab,
                                    infos = infos,
                                    paired = paired)
    else:
        wgsSim_bulk_parallel(prop_hc = prop_hc_tb,
                            coverage = current_cov,
                            num_clones = int_nodes,
                            raw_clone_prop = clone_prop,
                            threads = threads,
                            ls = chroms,
                            rl = read_len,
                            fl = frag_len,
                            floc = working_dir,
                            erate = error_rate,
                            tab = tab,
                            infos = infos,
                            paired = paired)
        
    # Aggregate all fastq.gz files from different threads
    if prefix is not None:
        leftfq_fname = prefix + '_bulkleft.fq.gz'
        rightfq_fname = prefix + '_bulkright.fq.gz'
    else:
        leftfq_fname = 'bulkleft.fq.gz'
        rightfq_fname = 'bulkright.fq.gz'
        
    aggregate_fastqs(working_dir, os.path.join(working_dir, leftfq_fname), os.path.join(working_dir, rightfq_fname), paired)

##################### SINGLE-CELL SIMULATION
else:
    # Compute coverage if coverage per cell is passed
    if coverage is None:
        current_cov = cell_coverage * num_single_cells
    else:
        current_cov = coverage
    if prefix is not None:
        params_fname = prefix + '_parameter_list.yaml'
    else:
        params_fname = 'parameter_list.yaml'
        
    with open(os.path.join(working_dir, 'parameter_list.yaml'), 'w') as f:
        f.write('num leaves: ' + str(num_clones)+'\n')
        f.write('num internal nodes: ' + str(int_nodes)+'\n')
        f.write('dir_conc: ' + str(alpha)+'\n')
        f.write('cell pop: ' + str(pop)+'\n')
        f.write('coverage: ' + str(coverage)+'\n')
        f.write('num single cells: ' + str(num_single_cells)+'\n')
        f.write('read len: ' + str(read_len)+'\n')
        f.write('frag len: ' + str(frag_len)+'\n')
        f.write('paired: ' + str(paired)+'\n')
        f.write('targeted: ' + str(targeted)+'\n')
        f.write('error rate: ' + str(error_rate) + '\n')
        f.write('proportions of healthy cells: ' + str(prop_hc_tb) + '\n')
        f.write('NB parameters: r=' + str(r) + ' p=' + str(p) + '\n')

    if targeted:
        targetedSim_sc_parallel(num_single_cells = num_single_cells,
                                prop_hc = prop_hc_tb,
                                coverage = current_cov,
                                num_clones = int_nodes,
                                clone_prop = clone_prop,
                                r = r,
                                p = p,
                                threads = threads,
                                ls = chroms,
                                rl = read_len,
                                fl = frag_len,
                                floc = working_dir,
                                regions = regions,
                                rev_numchrommap = rev_numchrommap,
                                erate = error_rate,
                                tab = tab,
                                infos = infos,
                                paired = paired)
            
    else:
        wgsSim_sc_parallel(num_single_cells = num_single_cells,
                            prop_hc = prop_hc_tb,
                            coverage = current_cov,
                            num_clones = int_nodes,
                            clone_prop = clone_prop,
                            r = r,
                            p = p,
                            threads = threads,
                            ls = chroms,
                            rl = read_len,
                            fl = frag_len,
                            floc = working_dir,
                            erate = error_rate,
                            tab = tab,
                            infos = infos,
                            paired = paired)
        
    # Aggregate all fastq.gz files from different cells
    if prefix is not None:
        leftfq_fname = prefix + '_scleft.fq.gz'
        rightfq_fname = prefix + '_scright.fq.gz'
    else:
        leftfq_fname = 'scleft.fq.gz'
        rightfq_fname = 'scright.fq.gz'
        
    aggregate_fastqs(working_dir, os.path.join(working_dir, leftfq_fname), os.path.join(working_dir, rightfq_fname), paired)

te = time.time()
if(bulk):
    print('Time elapsed for bulk biopsy simulation.', te-ts)
else:
    print('Time elapsed for single-cell biopsy simulation.', te-ts)