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

print('start')
ts = time.time()
with open(snakemake.params['yaml_config']) as f:
    params = yaml.safe_load(f)

storage_dir = snakemake.output['outdir']
reference_working_dir = os.path.join(storage_dir , 'reference')
full_genome = snakemake.input['genome']
signature_distributions = [float(1/params['num_signatures'])]*params['num_signatures']
sig_df = pd.read_csv(snakemake.input['signatures'], sep = '\t', header = None)
sig_df = sig_df.iloc[1:]
sig_df = sig_df.iloc[:,1:]
signatures_matrix = sig_df.to_numpy()
tab = bytes.maketrans(b"ACTG", b"TGAC")
read_len = params['read_len']
frag_len = params['frag_len']
ref_tot_nodes = 2*params['ref_clones']-1
ref_root_node = ref_tot_nodes-1
ref_int_nodes = ref_root_node-1
ref_alpha = random.choice(params['alpha_list'])
coverage = params['coverage_list'][snakemake.params['region']]
num_single_cells = params['num_single_cells']
paired = params['paired']
WES = params['WES']
error_rate = random.choice(params['error_rate_list'])
r = params['r']
p = params['p']
bulk = params['bulk']

chrom_names, chroms = map(list, zip(*((record.id, bytearray(str(record.seq).upper(), 'utf-8')) for record in SeqIO.parse(full_genome, "fasta")))) # Parse sequences and names together
numchrommap = dict(zip(range(len(chrom_names)), chrom_names))
rev_numchrommap = {v: k for k, v in numchrommap.items()}

if(WES):
    EXON_FILE = snakemake.input['binned_bed']
    total_num_intervals = 0
    exonDict = {}
    panel_chroms = [bytearray() for _ in chroms]
    strings_to_idx = []
    for i in numchrommap.values():
        exonDict[i] = []
    with open(EXON_FILE, 'r') as f:
        lines = [line.strip().split() for line in f]
    lines.sort(key=lambda x: (rev_numchrommap[x[0]], int(x[1]), int(x[2])))
    end_i = 0
    old_chrom = 'chr1'
    shift = 0
    for chrom, start, end in lines: # BED has to have 3 columns
        if chrom in numchrommap.values():
            if(chrom != old_chrom):
                panel_chroms[rev_numchrommap[old_chrom]] = bytearray().join(strings_to_idx)
                strings_to_idx.clear()
                end_i = 0
                old_chrom = chrom
                shift = 0
            old_end_i = end_i
            start_i = int(start) - read_len
            end_i = int(end) + read_len
            if(start_i < old_end_i):
                start_i = old_end_i
            if(end_i < old_end_i):
                end_i = old_end_i
                continue
            shift += start_i - old_end_i
            interval = [start_i - shift, end_i - shift]
            total_num_intervals += 1
            exonDict[chrom].append(interval)
            strings_to_idx.append(chroms[chrom_names.index(chrom)][start_i:end_i]) # Extract sequences corresponding to intervals in the BED file (add read length on the sides of the interval for sampling reads)
    # Finalize last chromosome
    if strings_to_idx:
        panel_chroms[rev_numchrommap[chrom]] = bytearray().join(strings_to_idx)
    chroms = panel_chroms
    del panel_chroms
    
os.makedirs(storage_dir, exist_ok=True)

# Mutation process rate lists (dependent on the length of the regions we are analysing)
genome_length = params['genome_length']
current_genome_length = sum(len(c) for c in chroms)
# corr_factor = (current_genome_length / genome_length) * (10 if snakemake.params['region'] == "on_target" else 1) # Model rate proportional to the panel length, since mutational events are considered independent across sites
corr_factor = 1.0
mut_events = params['mutational_events']
list_of_rates = {
    "SNV": [x * corr_factor for x in params['high_rates_list']],
    "CNV": [x * corr_factor for x in params['high_rates_list']],
    "DEL": [x * corr_factor for x in params['medium_rates_list']],
    "DELSMALL": [x * corr_factor for x in params['medium_rates_list']],
    "INVERSION": [x * corr_factor for x in params['low_rates_list']],
    "TRANSLOCATION": [x * corr_factor for x in params['low_rates_list']],
    "BFB": [x * corr_factor for x in params['ultralow_rates_list']],
    "CHROMOTHRIP": [x * corr_factor for x in params['ultralow_rates_list']],
    # "CHROMOPLEX": [x * corr_factor for x in params['ultralow_rates_list']],
    "INSERTIONSMALL": [x * corr_factor for x in params['medium_rates_list']],
    "KATAEGIS": [x * corr_factor for x in params['ultralow_rates_list']],
    "ANEUPLOIDY": [x * corr_factor for x in params['low_rates_list']]
}

if(False):
    # do reference first
    print('Making reference reads')
    os.makedirs(reference_working_dir, exist_ok=True)
    # clear_dir(reference_working_dir)
    # wgz(reference_working_dir + str(ref_root_node) + '.gz', chroms)
    with open(os.path.join(reference_working_dir, 'parameter_list.txt'), 'w') as f:
        f.write('coverage: ' + str(params['ref_coverage'])+'\n')
        f.write('read len: ' + str(ref_read_len)+'\n')
        f.write('frag len: ' + str(ref_frag_len)+'\n')
        f.write('paired: ' + str(params['ref_paired'])+'\n')
        f.write('WES: ' + str(params['ref_WES'])+'\n')
        f.write('error rate: ' + str(params['ref_erate']) + '\n')

    if(params['ref_paired']):
        if(params['ref_WES']):
            exonrunPairedSim(chroms, ref_int_nodes, params['ref_coverage'], ref_read_len, ref_frag_len, params['reference_working_dir'],
                            reference_working_dir, params['batch_size'], ref_root_node, exonDict, numchrommap, params['subblock_size'], ref_alpha, params['ref_erate'], tab, flag=1)
        else:
            runPairedSim(chroms, ref_int_nodes, params['ref_coverage'], ref_read_len, ref_frag_len, reference_working_dir,
                        reference_working_dir, params['batch_size'], ref_root_node, ref_alpha, params['ref_erate'], tab, flag=1)
    else:
        if(params['ref_WES']):
            exonrunSim(ref_int_nodes, params['ref_coverage'], ref_read_len, reference_working_dir, reference_working_dir,
                    params['batch_size'], ref_root_node, exonDict, numchrommap, params['subblock_size'], ref_alpha, params['ref_erate'], tab, flag=1)
        else:
            runSim(ref_int_nodes, params['ref_coverage'], ref_read_len, reference_working_dir,
                reference_working_dir, params['batch_size'], ref_root_node, ref_alpha, params['ref_erate'], flag=1)
    getmemory()


running_clone_list = []
num_tumors = params['num_tumors']
num_samples = params['num_samples']
num_clones_list = params['clone_list']
for num_clones in num_clones_list:
    alpha = params['alpha_list'][0]
    running_clone_list.append(num_clones)
    tot_nodes = 2*num_clones - 1
    root_node = tot_nodes - 1
    int_nodes = root_node - 1
    if(params['use_leaf_only']): 
        use_nodes = num_clones
    else: 
        use_nodes = int_nodes
    pop = params['pop_size']
    clone_number_dir = f'clone_{num_clones}'
    working_dir = os.path.join(storage_dir , str(clone_number_dir))
    os.makedirs(working_dir, exist_ok=True)
    tree = getTree(num_clones, pop, working_dir)
    list_of_paths, time_matrix, depth  = getPaths_and_TimeMatrix(tree, num_clones)
    mutationedge_list, avg_rate_list = generateOrder(tree, time_matrix, list_of_rates)
    infos, muts = saveMutations(
        chroms,
        tot_nodes,
        working_dir,
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
        tab)
    with open(os.path.join(working_dir, 'mutation_list.json'), 'w') as f:
        json.dump(muts, f, indent=4)
    with open(os.path.join(working_dir, 'information_list.json'), 'w') as f:
        json.dump(infos, f, indent=4)
    approx_len = len(muts[0])

    for sample in range(num_samples):
        sample_working_dir = os.path.join(working_dir, f'sample_{sample+1}')
        os.makedirs(sample_working_dir, exist_ok=True)
        with open(os.path.join(sample_working_dir, 'parameter_list.txt'), 'w') as f:
            f.write('num leaves: ' + str(num_clones)+'\n')
            f.write('dir_conc: ' + str(alpha)+'\n')
            f.write('cell pop: ' + str(pop)+'\n')
            f.write('coverage: ' + str(coverage)+'\n')
            f.write('num single cells: ' + str(num_single_cells)+'\n')
            f.write('read len: ' + str(read_len)+'\n')
            f.write('frag len: ' + str(frag_len)+'\n')
            f.write('paired: ' + str(paired)+'\n')
            f.write('WES: ' + str(WES)+'\n')
            f.write('rates of variants: ' + str(avg_rate_list)+'\n')
            f.write('full poisson time: ' + str(depth) + '\n')
            f.write('error rate: ' + str(error_rate) + '\n')
            if(num_single_cells > 1):
                f.write('NB parameters: r=' + str(r) + ' p=' + str(p) + '\n')

        if(WES):
            if bulk: # Bulk simulation
                targetedSim_bulk(ls = chroms,
                                 num_clones = use_nodes,
                                 coverage = coverage,
                                 rl = read_len,
                                 fl = frag_len,
                                 floc = sample_working_dir,
                                 batch = params['batch_size'],
                                 exonDict = exonDict,
                                 numchrommap = numchrommap,
                                 alpha = alpha,
                                 erate = error_rate,
                                 tab = tab,
                                 infos = infos,
                                 muts = muts,
                                 paired = paired)
            else:
                targetedSim_sc_parallel(num_single_cells,
                                        coverage,
                                        r,
                                        p,
                                        snakemake.threads,
                                        chroms,
                                        use_nodes,
                                        read_len,
                                        frag_len,
                                        sample_working_dir,
                                        params['batch_size'],
                                        exonDict,
                                        numchrommap,
                                        alpha,
                                        error_rate,
                                        tab,
                                        infos,
                                        muts,
                                        paired)
        else:            
            if bulk: # Bulk simulation
                wgsSim(ls = chroms,
                       num_clones = use_nodes,
                       coverage = coverage,
                       rl = read_len,
                       fl = frag_len,
                       floc = sample_working_dir,
                       batch = params['batch_size'],
                       alpha = alpha,
                       erate = error_rate,
                       tab = tab,
                       infos = infos,
                       muts = muts,
                       num_single_cells=1,
                       flag=0,
                       paired = paired)
            else:
                wgsSim(ls = chroms,
                       num_clones = use_nodes,
                       coverage = coverage,
                       rl = read_len,
                       fl = frag_len,
                       floc = sample_working_dir,
                       batch = params['batch_size'],
                       alpha = alpha,
                       erate = error_rate,
                       tab = tab,
                       infos = infos,
                       muts = muts,
                       r = r,
                       p = p,
                       num_single_cells=num_single_cells,
                       flag=0,
                       paired = paired)
                
        # Aggregate all fastq.gz files form different cell types
        if(not bulk):
            aggregate_fastqs(sample_working_dir, os.path.join(sample_working_dir, "bulkleft.fq.gz"), os.path.join(sample_working_dir, "bulkright.fq.gz"), paired)

print('finished tumors')
liquid_biopsy = params['liquid_biopsy']
if(liquid_biopsy):
    lb_dir = os.path.dir(storage_dir, 'liquid_biopsy')
    os.makedirs(lb_dir, exist_ok=True)
    # clear_dir(lb_dir)
    lb_alpha = random.choice(params['alpha_list'])
    lb_coverage = random.choice(params['coverage_list'])
    ctdna_frac = random.choice(params['ctdna_frac_list'])
    lbrunSim(num_tumors, running_clone_list, lb_coverage,
             storage_dir, lb_dir, ref_root_node, lb_alpha, ctdna_frac, 6)

te = time.time()
print('time elapsed', te-ts)
