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
# randomid = sys.argv[1]
with open(snakemake.params['yaml_config']) as f:
    params = yaml.safe_load(f)

storage_dir = snakemake.output['outdir']
# base_working_dir = storage_dir+'/' 
reference_working_dir = os.path.join(storage_dir , 'reference')
full_genome = snakemake.input['genome']
EXON_FILE = snakemake.input['binned_bed']
signature_distributions = [float(1/params['num_signatures'])]*params['num_signatures']
sig_df = pd.read_csv(snakemake.input['signatures'], sep = '\t', header = None)
sig_df = sig_df.iloc[1:]
sig_df = sig_df.iloc[:,1:]
signatures_matrix = sig_df.to_numpy()
tab = bytes.maketrans(b"ACTG", b"TGAC")
read_length_index = random.randint(0, len(params['read_len_list'])-1)
ref_read_len = params['read_len_list'][read_length_index]
ref_frag_len = params['frag_len_list'][read_length_index]
ref_tot_nodes = 2*params['ref_clones']-1
ref_root_node = ref_tot_nodes-1
ref_int_nodes = ref_root_node-1
ref_alpha = random.choice(params['alpha_list'])
# Mutation process rate lists
mut_events = params['mutational_events']
list_of_rates = {
    "SNV": params['ultrahigh_rates_list'],
    "CNV": params['high_rates_list'],
    "DEL": params['medium_rates_list'],
    "DELSMALL": params['medium_rates_list'],
    "INVERSION": params['low_rates_list'],
    "TRANSLOCATION": params['low_rates_list'],
    "BFB": params['ultralow_rates_list'],
    "CHROMOTHRIP": params['ultralow_rates_list'],
    # "CHROMOPLEX": params['ultralow_rates_list'],
    "INSERTIONSMALL": params['medium_rates_list'],
    "KATAEGIS": params['ultralow_rates_list'],
    "ANEUPLOIDY": params['high_rates_list']
}

# reversemap = {'chr1':0, 'chr10':2, 'chr11':4, 'chr12':6, 'chr13':8, 'chr14':10, 'chr15':12, 'chr16':14, 'chr17':16, 'chr18':18, 'chr19':20, 'chr2':22, 'chr20':24, 'chr21':26, 'chr22':28, 'chr3':30, 'chr4':32, 'chr5':34, 'chr6':36, 'chr7':38, 'chr8':40, 'chr9':42, 'chrX':44, 'chrY':45} 
numchrommap = {0: 'chr1', 1: 'chr1', 2: 'chr10', 3: 'chr10', 4: 'chr11', 5: 'chr11', 6: 'chr12', 7: 'chr12', 8: 'chr13', 9: 'chr13', 10: 'chr14', 11: 'chr14', 12: 'chr15', 13: 'chr15', 14: 'chr16', 15: 'chr16', 16: 'chr17', 17: 'chr17', 18: 'chr18', 19: 'chr18', 20: 'chr19', 21: 'chr19', 22: 'chr2', 23: 'chr2', 24: 'chr20', 25: 'chr20', 26: 'chr21', 27: 'chr21', 28: 'chr22', 29: 'chr22', 30: 'chr3', 31: 'chr3', 32: 'chr4', 33: 'chr4', 34: 'chr5', 35: 'chr5', 36: 'chr6', 37: 'chr6', 38: 'chr7', 39: 'chr7', 40: 'chr8', 41: 'chr8', 42: 'chr9', 43: 'chr9', 44: 'chrX', 45: 'chrY'}
# chrom_dict = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
#               'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY']
# reduced_chrom_dict = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
#                       'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']
# sex_chroms = ['chrX', 'chrY']
# chroms = [str(record.seq) for record in SeqIO.parse(full_genome, "fasta")]
# chrom_names = [record.id for record in SeqIO.parse(full_genome, "fasta")]
chrom_names, chroms = map(list, zip(*((record.id, bytearray(str(record.seq), 'utf-8')) for record in SeqIO.parse(full_genome, "fasta")))) # Parse sequences and names together

##### There is no need to iterate over the genome and rewrite it, you can prepare it before
# fasta_sequences = SeqIO.parse(open(full_genome), 'fasta')
# for fasta in fasta_sequences:
#     name = fasta.id
#     if (name in reduced_chrom_dict):
#         # print(fasta.name)
#         chroms.append(str(fasta.seq).upper())
#         chroms.append(str(fasta.seq).upper())
#     elif(name in sex_chroms):
#         # print(fasta.name)
#         chroms.append(str(fasta.seq).upper())
#     else:
#         continue
# print('chroms loaded')
# getmemory()

total_num_intervals = 0
exonDict = {}
strings_to_idx = []
for i in numchrommap.values():
    exonDict[i] = []
with open(EXON_FILE, 'r') as f:
    for line in f:
        chrom, start, end = line.split() # BED has to have 3 columns
        if chrom in numchrommap.values():
            start_i = int(start)
            end_i = int(end)
            interval = [start_i, end_i]
            total_num_intervals += 1
            exonDict[chrom].append(interval)
            strings_to_idx.append(chroms[chrom_names.index(chrom)][start_i:end_i]) # Extract sequences corresponding to intervals in the BED file
# del chroms
getmemory()
# print(len(strings_to_idx))

os.makedirs(storage_dir, exist_ok=True)
# clear_dir(base_working_dir)
# print(base_working_dir)

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
# num_tumors = random.choice(params['num_tumors_list'])
# num_samples = random.choice(params['num_samples_list'])
num_tumors = params['num_tumors_list'][0]
num_samples = params['num_samples_list'][0]
num_clones_list = params['clone_list']
for num_clones in num_clones_list:
    getmemory()
    # alpha = random.choice(params['alpha_list'])
    # num_clones = random.choice(params['clone_list'])
    alpha = params['alpha_list'][0]
    # num_clones = params['clone_list'][0]
    running_clone_list.append(num_clones)
    tot_nodes = 2*num_clones - 1
    root_node = tot_nodes - 1
    int_nodes = root_node - 1
    if(params['use_leaf_only']): 
        use_nodes = num_clones
    else: 
        use_nodes = int_nodes
    # pop = random.choice(params['pop_list'])
    pop = params['pop_list'][0]
    clone_number_dir = f'clone_{num_clones}'
    working_dir = os.path.join(storage_dir , str(clone_number_dir))
    os.makedirs(working_dir, exist_ok=True)
    # clear_dir(working_dir)
    # baseline_chroms = rgz(reference_working_dir + str(ref_root_node) + '.gz')
    # wgz(working_dir + str(root_node) + '.gz', baseline_chroms)
    # del baseline_chroms
    # print('written reg')
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
    # with open(os.path.join(working_dir, 'mutation_list.txt'), 'w') as f:
    #     f.write(str(muts))
    # with open(os.path.join(working_dir, 'information_list.txt'), 'w') as f:
    #     f.write(str(infos))
    approx_len = len(muts[0])
    # print('approx mutations', approx_len)
    # print(muts)
    # print('mutated genomes stored')
    getmemory()

    for sample in range(num_samples):
        print('starting sample')
        exonDictR = exonDict
        real_working_dir = os.path.join(working_dir, f'sample_{sample}')
        os.makedirs(real_working_dir, exist_ok=True)
        # clear_dir(real_working_dir)
        # coverage = random.choice(params['coverage_list'])
        coverage = params['coverage_list'][snakemake.params['region']]
        num_single_cells = params['num_single_cell_list'][0]
        read_length_index = random.randint(0, len(params['read_len_list'])-1)
        read_len = params['read_len_list'][read_length_index]
        frag_len = params['frag_len_list'][read_length_index]
        # paired = random.choice(params['paired_list'])
        paired = params['paired_list'][0]
        # WES = random.choice(params['WES_list'])
        WES = params['WES_list'][0]
        error_rate = random.choice(params['error_rate_list'])
        with open(os.path.join(real_working_dir, 'parameter_list.txt'), 'w') as f:
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
        getmemory()
        if(paired):
            if(WES):
                if num_single_cells == 0: # Bulk simulation
                    exonrunPairedSim(chroms, use_nodes, coverage, read_len, frag_len, working_dir, real_working_dir,
                                    params['batch_size'], root_node, exonDictR, numchrommap, params['subblock_size'], alpha, error_rate, tab, infos, muts, num_single_cells=1, flag=0)
                for i in range(num_single_cells):
                    # single_cell_dir = working_dir+f'samplenum_{sample}_singlecell_{i}/'
                    # os.makedirs(single_cell_dir, exist_ok=True)
                    # clear_dir(single_cell_dir)
                    exonrunPairedSim(chroms, use_nodes, coverage, read_len, frag_len, working_dir, real_working_dir,
                                     params['batch_size'], root_node, exonDictR, numchrommap, params['subblock_size'], alpha, error_rate, tab, infos, muts, num_single_cells=num_single_cells, flag=2)
                    # print(f"Memory for cell {i}")
                    # getmemory()
                    # sys.stdout.flush()
                    
            else:
                if num_single_cells == 0: # Bulk simulation
                    runPairedSim(use_nodes, coverage, read_len, frag_len, working_dir,
                                real_working_dir, params['batch_size'], root_node, alpha, error_rate, tab, flag=0)
                for i in range(num_single_cells):
                    # single_cell_dir = working_dir+f'samplenum_{sample}_singlecell_{i}/'
                    # os.makedirs(single_cell_dir, exist_ok=True)
                    # clear_dir(single_cell_dir) 
                    runPairedSim(use_nodes, coverage, read_len, frag_len, working_dir,
                                 real_working_dir, params['batch_size'], root_node, alpha, error_rate, tab, flag=2)
        else:
            if(WES):
                if num_single_cells == 0: # Bulk simulation
                    exonrunSim(use_nodes, coverage, read_len, working_dir, real_working_dir,
                            params['batch_size'], root_node, exonDictR, numchrommap, params['subblock_size'], alpha, error_rate, tab, flag=0)
                for i in range(num_single_cells):
                    # single_cell_dir = working_dir+f'samplenum_{sample}_singlecell_{i}/'
                    # os.makedirs(single_cell_dir, exist_ok=True)
                    # clear_dir(single_cell_dir) 
                    exonrunSim(use_nodes, coverage, read_len, working_dir, real_working_dir,
                               params['batch_size'], root_node, exonDictR, numchrommap, params['subblock_size'], alpha, error_rate, tab, flag=2)

            else:
                if num_single_cells == 0: # Bulk simulation
                    runSim(use_nodes, coverage, read_len, working_dir,
                        real_working_dir, params['batch_size'], root_node, alpha, error_rate, tab, flag=0)
                for i in range(num_single_cells):
                    # os.makedirs(single_cell_dir, exist_ok=True)
                    # clear_dir(single_cell_dir) 
                    runSim(use_nodes, coverage, read_len, working_dir,
                           real_working_dir, params['batch_size'], root_node, alpha, error_rate, tab, flag=2)
print('finished tumors')
liquid_biopsy = random.choice(params['liquid_biopsy_list'])
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
