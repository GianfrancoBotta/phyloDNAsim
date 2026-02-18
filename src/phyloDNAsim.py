import psutil
import os
import msprime
import tskit
import numpy as np
import gzip
import glob
import string
import random
import time
import pickle
import sys
import pandas as pd
from Bio import SeqIO
from datasketch import *
import math
import re
import yaml

from simulation_utils import *

# Redirect stdout to log
log = open(snakemake.log[0], "a")
sys.stdout = log


print('Start simulation for biopsy data.')
ts_sb = time.time()
with open(snakemake.params['yaml_config']) as f:
    params = yaml.safe_load(f)

storage_dir = snakemake.output['outdir']
full_genome = snakemake.input['genome']
signature_distributions = [float(1/params['num_signatures'])]*params['num_signatures']
sig_df = pd.read_csv(snakemake.input['signatures'], sep = '\t', header = None)
sig_df = sig_df.iloc[1:]
sig_df = sig_df.iloc[:,1:]
signatures_matrix = sig_df.to_numpy()
tab = bytes.maketrans(b"ACTG", b"TGAC")
read_len = params['read_len']
frag_len = params['frag_len']
coverage = params['coverage']
num_single_cells = params['num_single_cells']
paired = params['paired']
targeted = params['targeted']
error_rate = random.choice(params['error_rate_list'])
r = params['r']
p = params['p']
bulk = params['bulk']
liquid_biopsy = params['liquid_biopsy']
prop_hc = params['prop_hc_sb']

chrom_names, chroms = map(list, zip(*((record.id, bytearray(str(record.seq).upper(), 'utf-8')) for record in SeqIO.parse(full_genome, "fasta")))) # Parse sequences and names together
numchrommap = dict(zip(range(len(chrom_names)), chrom_names))
rev_numchrommap = {v: k for k, v in numchrommap.items()}

if targeted:
    EXON_FILE = snakemake.input['bed']
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
    
os.makedirs(storage_dir, exist_ok=True)

# Mutation process rate lists (dependent on the length of the regions we are analysing)
genome_length = params['genome_length']
mut_events = params['mutational_events']
list_of_rates = {
    "SNV": params['high_rates_list'],
    "CNV": params['high_rates_list'],
    "DEL": params['low_rates_list'],
    "DELSMALL": params['medium_rates_list'],
    "INVERSION": params['low_rates_list'],
    "TRANSLOCATION": params['low_rates_list'],
    "BFB": params['ultralow_rates_list'],
    "CHROMOTHRIP": params['ultralow_rates_list'],
    # "CHROMOPLEX": params['ultralow_rates_list'],
    "INSERTIONSMALL": params['medium_rates_list'],
    "KATAEGIS": params['ultralow_rates_list'],
    "ANEUPLOIDY": params['ultralow_rates_list']
}
if targeted:
    list_of_rates = {k: [float(r) / 10 for r in v] for k, v in list_of_rates.items()}

num_tumors = params['num_tumors']
num_samples = params['num_samples']
num_clones_list = params['clone_list']
for num_clones in num_clones_list:
    alpha = params['alpha_list'][0]
    tot_nodes = 2*num_clones - 1
    int_nodes = tot_nodes - 1
    pop = params['pop_size']
    clone_number_dir = f'clone_{num_clones}'
    working_dir = os.path.join(storage_dir , clone_number_dir)
    os.makedirs(working_dir, exist_ok=True)
    tree = getTree(num_clones, pop, working_dir)
    list_of_paths, time_matrix, depth  = getPaths_and_TimeMatrix(tree, num_clones)
    mutationedge_list, avg_rate_list = generateOrder(tree, time_matrix, list_of_rates)
    infos = saveMutations(
        chroms,
        tot_nodes,
        list_of_paths,
        params['use_signatures'],
        mutationedge_list,
        params['num_signatures'],
        params['signature_alpha'],
        signature_distributions,
        signatures_matrix,
        numchrommap,
        params['list_of_bases'],
        params['list_of_pairs'],
        tab,
        targeted,
        regions)
    with open(os.path.join(working_dir, 'information_list.json'), 'w') as f:
        json.dump(infos, f, indent=4)

    for sample in range(num_samples):
        sample_working_dir = os.path.join(working_dir, f'sample_{sample+1}')
        os.makedirs(sample_working_dir, exist_ok=True)
        with open(os.path.join(sample_working_dir, 'parameter_list.yaml'), 'w') as f:
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
            f.write('rates of variants: ' + str(avg_rate_list)+'\n')
            f.write('full poisson time: ' + str(depth) + '\n')
            f.write('error rate: ' + str(error_rate) + '\n')
            f.write('proportions of healthy cells: ' + str(prop_hc) + '\n')
            if(not bulk):
                f.write('NB parameters: r=' + str(r) + ' p=' + str(p) + '\n')

        if targeted:
            if bulk: # Bulk simulation
                clone_prop = targetedSim_bulk_parallel(prop_hc = prop_hc,
                                                       coverage = coverage,
                                                       num_clones = int_nodes,
                                                       alpha = alpha,
                                                       threads = snakemake.threads,
                                                       ls = chroms,
                                                       rl = read_len,
                                                       fl = frag_len,
                                                       floc = sample_working_dir,
                                                       regions = regions,
                                                       rev_numchrommap = rev_numchrommap,
                                                       erate = error_rate,
                                                       tab = tab,
                                                       infos = infos,
                                                       paired = paired)
            else:
                clone_prop = targetedSim_sc_parallel(num_single_cells = num_single_cells,
                                                     prop_hc = prop_hc,
                                                     coverage = coverage,
                                                     num_clones = int_nodes,
                                                     alpha = alpha,
                                                     r = r,
                                                     p = p,
                                                     threads = snakemake.threads,
                                                     ls = chroms,
                                                     rl = read_len,
                                                     fl = frag_len,
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
                    coverage = coverage,
                    rl = read_len,
                    fl = frag_len,
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
                    coverage = coverage,
                    rl = read_len,
                    fl = frag_len,
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

print('Start simulation for liquid biopsy data.')
ts_lb = time.time()

read_len_lb = params['read_len_lb']
frag_len_lb = params['frag_len_lb']
prop_hc = params['prop_hc_lb']

if(liquid_biopsy):
    for num_clones in num_clones_list:
        alpha = params['alpha_list'][0]
        tot_nodes = 2*num_clones - 1
        int_nodes = tot_nodes - 1
        clone_number_dir = f'clone_{num_clones}'
        working_dir = os.path.join(storage_dir , clone_number_dir)
        for sample in range(num_samples):
                sample_working_dir = os.path.join(working_dir, f'sample_{sample+1}')
                os.makedirs(sample_working_dir, exist_ok=True)
                with open(os.path.join(sample_working_dir, 'parameter_list_lb.yaml'), 'w') as f:
                    f.write('num leaves: ' + str(num_clones)+'\n')
                    f.write('num internal nodes: ' + str(int_nodes)+'\n')
                    f.write('dir_conc: ' + str(alpha)+'\n')
                    f.write('cell pop: ' + str(pop)+'\n')
                    f.write('coverage: ' + str(coverage)+'\n')
                    f.write('num single cells: ' + str(num_single_cells)+'\n')
                    f.write('read len: ' + str(read_len_lb)+'\n')
                    f.write('frag len: ' + str(frag_len_lb)+'\n')
                    f.write('paired: ' + str(paired)+'\n')
                    f.write('targeted: ' + str(targeted)+'\n')
                    f.write('rates of variants: ' + str(avg_rate_list)+'\n')
                    f.write('full poisson time: ' + str(depth) + '\n')
                    f.write('error rate: ' + str(error_rate) + '\n')
                    f.write('proportions of healthy cells: ' + str(prop_hc) + '\n')

                if targeted:
                    clone_prop = targetedSim_bulk_parallel(prop_hc = prop_hc,
                                                           coverage = coverage,
                                                           num_clones = int_nodes,
                                                           alpha = alpha,
                                                           threads = snakemake.threads,
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
                
                aggregate_fastqs(sample_working_dir, os.path.join(sample_working_dir, "ctleft.fq.gz"), os.path.join(sample_working_dir, "ctright.fq.gz"), paired)

te_lb = time.time()
print('Time elapsed for liquid biopsy simulation', te_lb-ts_lb)
