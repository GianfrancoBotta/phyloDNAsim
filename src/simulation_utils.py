import copy
import gc
import glob
import gzip
import io
import json
import math
import msprime
import multiprocessing as mp
import numpy as np
import os
import psutil
import pandas as pd
import random
import re
import shutil
import statistics
import string
import tskit

from create_mutations import *
from apply_mutations import *
from adapt_targeted_regions import *

def generateMatrix(the_matrix, list_of_weights, time_matrix):
    rate_matrix = getRateMatrix(time_matrix, list_of_weights)
    pop = time_matrix.shape[0]
    avg_rate = 0.0
    num_rates = 0
    for i in range(pop):
        for j in range(pop):
            if(time_matrix[i, j] != 0):
                time = 0
                list_of_times = []
                while(time < time_matrix[i, j]):
                    num_rates += 1
                    avg_rate += rate_matrix[i, j]
                    mutations_per_branch = rate_matrix[i, j]*time_matrix[i, j]
                    wait_time = np.random.exponential(1/rate_matrix[i, j])
                    time += wait_time
                    if(time < time_matrix[i, j]):
                        list_of_times.append(time)
                the_matrix[i, j] = list_of_times
    avg_rate = avg_rate/num_rates
    return the_matrix, avg_rate


def getTree(num_clones, pop, working_dir):
    # tree_sequence = msprime.sim_ancestry(samples=num_clones, population_size=pop, recombination_rate=0, sequence_length=seq_len, random_seed=9, ploidy=1)
    tree_sequence = msprime.simulate(sample_size=num_clones, Ne=pop, recombination_rate=0)
    tree_sequence.dump(working_dir + '/tree_sequence.tree')
    tree = tree_sequence.first()
    return tree


def getPaths_and_TimeMatrix(tree, num_clones):
    '''
    Returns a list of paths (also lists) from leaf to root.
    '''
    list_of_paths = []
    time_matrix = np.zeros((tree.root + 1, tree.root+1))
    dep = 0
    for u in range(num_clones):
        path = []
        while u != tskit.NULL:
            path.insert(0, u)
            time_matrix[u, tree.parent(u)] = tree.get_branch_length(u)
            dep += tree.get_branch_length(u)
            u = tree.parent(u)
        list_of_paths.append(path)
    time_matrix = np.transpose(time_matrix)
    return list(list_of_paths), time_matrix, dep

def getRateMatrix(time_matrix, list_of_weights):
    the_shape = np.shape(time_matrix)[0]
    rate_matrix = np.zeros((the_shape, the_shape))
    for i in range(the_shape):
        for j in range(the_shape):
            if(time_matrix[i, j] > 0):
                rate_matrix[i, j] = np.random.choice(list_of_weights)
    return rate_matrix

def generateOrder(tree, time_matrix, list_of_rates):
    pop = time_matrix.shape[0]
    init_matrix = np.zeros((pop, pop), dtype=object) # pass a copy to avoid shared-object aliasing
    SNV_matrix, snv_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['SNV'], time_matrix)
    CNV_matrix, cnv_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['CNV'], time_matrix)
    DELETION_matrix, deletion_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['DEL'], time_matrix)
    DELETIONSMALL_matrix, delsmall_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['DELSMALL'], time_matrix)
    INVERSION_matrix, inv_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['INVERSION'], time_matrix)
    TRANSLOCATION_matrix, trans_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['TRANSLOCATION'], time_matrix)
    BFB_matrix, bfb_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['BFB'], time_matrix)
    CHROMOTHRIPSIS_matrix, chromothripsis_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['CHROMOTHRIP'], time_matrix)
    # CHROMOPLEXY_matrix, chromoplexy_rate = generateMatrix(
    #     init_matrix.copy(), list_of_rates['CHROMOPLEX'], time_matrix)
    INSERTIONSMALL_matrix, inssmall_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['INSERTIONSMALL'], time_matrix)
    KATAEGIS_matrix, kat_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['KATAEGIS'], time_matrix)
    ANEUPLOIDY_matrix, an_rate = generateMatrix(
        init_matrix.copy(), list_of_rates['ANEUPLOIDY'], time_matrix)
    mutationedge_list = []
    fin_rate_list = [snv_rate, cnv_rate, deletion_rate, delsmall_rate, inv_rate, trans_rate,
                     bfb_rate, chromothripsis_rate, inssmall_rate, kat_rate, an_rate] # add chromoplexy_rate
    for i in range(pop-1):
        # Gives mutations from i to its parent ONLY
        parent_node = tree.parent(i)
        mutation_lists = [
            SNV_matrix[parent_node, i],
            CNV_matrix[parent_node, i],
            DELETION_matrix[parent_node, i],
            DELETIONSMALL_matrix[parent_node, i],
            INVERSION_matrix[parent_node, i],
            TRANSLOCATION_matrix[parent_node, i],
            BFB_matrix[parent_node, i],
            CHROMOTHRIPSIS_matrix[parent_node, i],
            # CHROMOPLEXY_matrix[parent_node, i],
            INSERTIONSMALL_matrix[parent_node, i],
            KATAEGIS_matrix[parent_node, i],
            ANEUPLOIDY_matrix[parent_node, i],
        ]
        ordered_muts = {}
        for idx, lst in enumerate(mutation_lists):
            ordered_muts.update({list(list_of_rates.keys())[idx]: lst})
        mutationedge_list.append(ordered_muts)
    return mutationedge_list, fin_rate_list

def saveMutations(current_genome, tot_nodes, list_of_paths, use_signatures, mutationedge_list, num_signatures, signature_alpha, signature_distributions, signatures_matrix, numchrommap, list_of_bases, list_of_pairs, tab, targeted, regions):
    save_functions = {
        'SNV': (create_SNPSig) if use_signatures else (create_speedSNP),
        'CNV': create_CNV,
        'DEL': create_deletion,
        'DELSMALL': create_deletionsmall,
        'INVERSION': create_inversion,
        'TRANSLOCATION': create_translocation,
        'BFB': create_BFB,
        'CHROMOTHRIP': create_chromothripsis,
        # 'CHROMOPLEX': create_chromoplexy,
        'INSERTIONSMALL': create_insertionsmall,
        'KATAEGIS': create_kataegis,
        'ANEUPLOIDY': create_aneuploidy
    }
    apply_functions = {
        'SNV': apply_SNP,
        'CNV': apply_CNV,
        'DEL': apply_deletion,
        'DELSMALL': apply_deletion,
        'INVERSION': apply_inversion,
        'TRANSLOCATION': apply_translocation,
        'BFB': apply_BFB,
        'CHROMOTHRIP':  apply_chromothripsis,
        # 'CHROMOPLEX': apply_chromoplexy,
        'INSERTIONSMALL': apply_insertion,
        'KATAEGIS': apply_kataegis,
        'ANEUPLOIDY': apply_aneuploidy
    }
    infos = {}
    infos[tot_nodes - 1] = []
    for path in list_of_paths:
        mutated_genome = copy.deepcopy(current_genome) # Create a copy not to modify directly the genome
        for i in range(len(path)-1):
            if(path[i+1] != (tot_nodes-1)):
                c_infos = infos[path[i]].copy()
                current_muts = mutationedge_list[path[i+1]]
                
                for m_type in current_muts.keys():
                    for _ in current_muts[m_type]:
                        if(m_type == "SNV" and use_signatures):
                            info = save_functions[m_type](mutated_genome, num_signatures, signature_alpha, signature_distributions, signatures_matrix, numchrommap, list_of_bases, list_of_pairs, tab, targeted, regions)
                        else:
                            info = save_functions[m_type](mutated_genome, numchrommap, targeted, regions)
                        if info is None:
                            print("Invalid information to apply mutations.")
                        c_infos.append(info)
                        # Adjust mutations info for the mutated genomes (shift indices if there are two or more events on the same chromosome)
                        mutated_genome = apply_functions[m_type](mutated_genome, info)
                infos[path[i+1]] = c_infos
    return infos

def applyMutations(seqs, infos, clone):
    '''
    Creates a new genome for the given clone.
    '''
    info = infos[clone]
    functions = {
        'SNV': apply_SNP,
        'CNV': apply_CNV,
        'DEL': apply_deletion,
        'DELSMALL': apply_deletion,
        'INVERSION': apply_inversion,
        'TRANSLOCATION': apply_translocation,
        'BFB': apply_BFB,
        'CHROMOTHRIP':  apply_chromothripsis,
        # 'CHROMOPLEX': apply_chromoplexy,
        'INSERTIONSMALL': apply_insertion,
        'KATAEGIS': apply_kataegis,
        'ANEUPLOIDY': apply_aneuploidy
    }
    
    clone_seqs = copy.deepcopy(seqs)
    for i in info:
        m = i['event']
        clone_seqs = functions[m](clone_seqs, i)
        
    return clone_seqs

def adapt_targeted_regions(regions, infos, clone):
    '''
    Adapts the targeted regions based on the mutations.
    '''
    info = infos[clone]
    functions = {
        'CNV': adapt_CNV,
        'DEL': adapt_deletion,
        'DELSMALL': adapt_deletion,
        'TRANSLOCATION': adapt_translocation,
        'BFB': adapt_BFB,
        'CHROMOTHRIP':  adapt_chromothripsis,
        # 'CHROMOPLEX': adapt_chromoplexy,
        'INSERTIONSMALL': adapt_insertion,
        'ANEUPLOIDY': adapt_aneuploidy
    }
    
    clone_regions = copy.deepcopy(regions)
    for i in info:
        m = i['event']
        if m in ['SNV', 'INVERSION', 'KATAEGIS']:
            continue
        clone_regions = functions[m](clone_regions, i)
        
    return clone_regions
    
    
def getfrag(r):
    return_val = 0
    while(return_val < r-20):
        return_val = np.random.negative_binomial(r, 0.5)
    return return_val

def drawPoisson(r):
    return np.random.poisson(r)

def getDNAchunk(length, seqs):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    seq = seqs[chrom[0]]
    n = len(seq)
    stidx = random.randint(0, n-1)
    edidx = stidx + length
    subseq = seq[stidx:edidx]
    return subseq

def mutateFrag(frag, err_rate, bases=['A', 'C', 'T', 'G']):
    num_muts = 0
    n = len(frag)
    q = 1-err_rate
    if (err_rate > 0.1):
        num_muts = int(random.gauss(n*err_rate, n*err_rate*q))
    else:
        num_muts = drawPoisson(n*err_rate)
    s = frag
    for i in range(num_muts):
        pos = random.randint(0, n-1)
        s[pos] = ord(random.choice(list(set(bases)-set(list(chr(frag[pos]))))))
    return s

def split(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def revc(sequence, tab):
    return sequence.translate(tab)[::-1]

def wgsSim(ls, num_clones, coverage, rl, fl, floc, batch, alpha, erate, tab, infos, r=None, p=None, flag=0, num_single_cells = 1, paired=False):
    # Open files to write
    if(flag == 0):
        f1 = gzip.open(os.path.join(floc, 'bulkleft.fq.gz'), 'wt')
        if(paired):
            f2 = gzip.open(os.path.join(floc, 'bulkright.fq.gz'), 'wt')
    elif(flag == 1):
        f1 = gzip.open(os.path.join(floc, 'refleft.fq.gz'), 'wt')
        if(paired):
            f2 = gzip.open(os.path.join(floc, 'refright.fq.gz'), 'wt')
    else:
        f1 = gzip.open(os.path.join(floc, 'singlecellleft.fq.gz'), 'wt')
        if(paired):
            f2 = gzip.open(os.path.join(floc, 'singlecellright.fq.gz'), 'wt')
    # Initialize coverage
    cov = 0.0
    if(paired):
        ratio = 2*rl/fl
    else:
        fl = rl
    if(flag == 0):
        chroms = ls # Duplicate healthy chromosomes to modify them
    if(flag == 2):
        if(any(x is None for x in [p, r])):
            print("Negative binomial parameters for single-cell depths are needed in single-cell mode.")
        target_cov = np.random.negative_binomial(r, p, size=num_single_cells) # Sample cell depths from a negative binomial
        target_cov = target_cov / target_cov.sum() * coverage # Scale to match total coverage
    else:
        target_cov = list(coverage)
    # Initialize lists to store reads from each single-cell
    r1 = []
    if(paired):
        r2 = []
    for i in range(num_single_cells):
        if(flag == 2): # Pick a clone for every single-cell if it is a single-cell simultation
            distn = getDirichletClone(num_clones, alpha)
            clone = pickdclone(distn, num_clones)
            chroms = applyMutations(ls, infos, clone)
        while(cov < target_cov[i]):
            if(flag == 1): # Pick a clone for every added coverage if it is a bulk simulation
                distn = getDirichletClone(num_clones, alpha)
                clone = pickdclone(distn, num_clones)
                chroms = applyMutations(ls, infos, clone)
            chromnum = 0
            for j in range(batch):
                for chrom in chroms:
                    if(len(chrom)==0): # Deleted chromosomes
                        continue
                    # frag_len = 0
                    # while(frag_len <= rl):
                    #     frag_len = getfrag(fl)
                    for seq in split(chrom, frag_len):
                        sub = seq
                        if random.random() > 0.5:
                            sub = revc(seq, tab)
                        random_str = ''.join(random.choices(
                            string.ascii_letters, k=15))
                        pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
                        qual1 = 'K'*len(pair1)
                        if(paired):
                            pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
                            qual2 = 'K'*len(pair2)
                        if(flag == 2):
                            r1.append('\n'.join([f'@cell{clone}_{random_str}', pair1, '+', qual1]) + '\n')
                            if(paired):    
                                r2.append('\n'.join([f'@cell{clone}_{random_str}', pair2, '+', qual2]) + '\n')
                        else:
                            r1.append('\n'.join([f'@{random_str}', pair1, '+', qual1]) + '\n')
                            if(paired):    
                                r2.append('\n'.join([f'@{random_str}', pair2, '+', qual2]) + '\n')
                    chromnum += 1
            # Update coverage after each iteration
            cov += (2*batch*ratio) if paired else 2*batch
        # Write to files
        f1.write("".join(r1))
        if(paired):
            f2.write("".join(r2))
        r1.clear()
        if(paired):
            r2.clear()
    # Close files
    f1.close()
    if(paired):
        f2.close()
    return(0)

def targetedSim_bulk(thread_id, clone_prop, target_cov, num_clones, ls, rl, fl, floc, regions, rev_numchrommap, erate, tab, infos, paired=False):
    # Per-thread output directory or prefix
    thread_dir = os.path.join(floc, f"thread_{thread_id}")
    os.makedirs(thread_dir, exist_ok=True)
    # Open files to write
    f1 = io.BufferedWriter(gzip.open(os.path.join(thread_dir, 'threadleft.fq.gz'), 'wb'), buffer_size = 4 * 1024**2)
    if(paired):
        f2 = io.BufferedWriter(gzip.open(os.path.join(thread_dir, 'threadright.fq.gz'), 'wb'), buffer_size = 4 * 1024**2)
    else:
        fl = rl
        
    # Initialize coverage
    cov = 0.0
    
    for clone in range(num_clones+1):
        clone_target_cov = target_cov * clone_prop[clone]
        chroms = applyMutations(ls, infos, clone)
        mod_regions = adapt_targeted_regions(regions, infos, clone)
        
        # Compute panel length to scale coverage and read quality and sampling weights
        panel_size = (mod_regions['end'] - mod_regions['start']).sum()
        weights = mod_regions['end'] - mod_regions['start']
        
        while cov < clone_target_cov:
            region = mod_regions.sample(weights=weights).iloc[0].tolist()
            chrom = chroms[rev_numchrommap[region[0]]]
            if len(chrom) == 0:
                continue
            if((region[2]-region[1]) < fl):
                startindex = region[1]
            else:
                startindex = random.randint(region[1], region[2]-fl)
                            
            # Update coverage in each iteration
            if((region[2]-startindex) < rl):
                cov += (2 * (region[2]-startindex) / panel_size) if paired else ((region[2]-startindex) / panel_size)
            else:
                cov += (2 * rl / panel_size) if paired else (rl / panel_size)
        
            # Extrct the read
            if((region[2]-region[1]) < fl and (region[2]-region[1]) > rl):
                sub = chrom[startindex:startindex + region[2]-region[1]]
            elif((region[2]-region[1]) < rl):
                sub = chrom[startindex:startindex + rl]
            else:
                sub = chrom[startindex:startindex + fl]
                
            if random.random() > 0.5:
                sub = revc(sub, tab)
            random_str = ''.join(random.choices(string.ascii_letters, k=15))
            pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
            qual1 = 'K'*len(pair1)
            if(paired):
                pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
                qual2 = 'K'*len(pair2)
            f1.write(('\n'.join([f'@clone{clone}_{random_str}', pair1, '+', qual1]) + '\n').encode('utf-8'))
            if(paired):
                f2.write(('\n'.join([f'@clone{clone}_{random_str}', pair2, '+', qual2]) + '\n').encode('utf-8'))
    # Close files
    f1.close()
    if(paired):
        f2.close()
    return(0)

def targetedSim_sc(cell_id, clone, target_cov, ls, rl, fl, floc, regions, rev_numchrommap, erate, tab, infos, paired=False):
    # Per-cell output directory or prefix
    cell_dir = os.path.join(floc, f"cell_{cell_id}")
    os.makedirs(cell_dir, exist_ok=True)
    # Open files to write
    f1 = io.BufferedWriter(gzip.open(os.path.join(cell_dir, 'scleft.fq.gz'), 'wb'), buffer_size = 4 * 1024**2)
    if paired:
        f2 = io.BufferedWriter(gzip.open(os.path.join(cell_dir, 'scright.fq.gz'), 'wb'), buffer_size = 4 * 1024**2)
        
    # Initialize coverage
    cov = 0.0
    chroms = applyMutations(ls, infos, clone)
    mod_regions = adapt_targeted_regions(regions, infos, clone)
    
    # Compute panel length to scale coverage and read quality and sampling weights
    panel_size = (mod_regions['end'] - mod_regions['start']).sum()
    weights = mod_regions['end'] - mod_regions['start']
    
    while cov < target_cov:
        region = mod_regions.sample(weights=weights).iloc[0].tolist()
        chrom = chroms[rev_numchrommap[region[0]]]
        if len(chrom) == 0:
            continue
        if((region[2]-region[1]) < fl):
            startindex = region[1]
        else:
            startindex = random.randint(region[1], region[2]-fl)
            
        # Update coverage in each iteration
        if((region[2]-startindex) < rl):
            cov += (2 * (region[2]-startindex) / panel_size) if paired else ((region[2]-startindex) / panel_size)
        else:
            cov += (2 * rl / panel_size) if paired else (rl / panel_size)
        
        # Extrct the read
        if((region[2]-region[1]) < fl and (region[2]-region[1]) > rl):
            sub = chrom[startindex:startindex + region[2]-region[1]]
        elif((region[2]-region[1]) < rl):
            sub = chrom[startindex:startindex + rl]
        else:
            sub = chrom[startindex:startindex + fl]
        
        if random.random() > 0.5:
            sub = revc(sub, tab)
        random_str = ''.join(random.choices(string.ascii_letters, k=15))
        pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
        qual1 = 'K'*len(pair1)
        if(paired):
            pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
            qual2 = 'K'*len(pair2)
        f1.write(('\n'.join([f'@cell{cell_id}_clone{clone}_{random_str}', pair1, '+', qual1]) + '\n').encode('utf-8'))
        if(paired):
            f2.write(('\n'.join([f'@cell{cell_id}_clone{clone}_{random_str}', pair2, '+', qual2]) + '\n').encode('utf-8'))
    f1.close()
    if paired:
        f2.close()
    return(0)

def targetedSim_bulk_parallel(prop_hc, coverage, num_clones, alpha=None, clone_prop=None, threads=None, **kwargs):
    if(threads == None):
        print("Using 4 cores as number of cores to use is not specified.")
        threads = 4
    thread_cov = coverage / threads # Scale to match total coverage
    
    # Compute tumor clones proportions
    if clone_prop is None:
        raw_clone_prop = getDirichletClone(num_clones, alpha)
    else:
        raw_clone_prop = clone_prop
    
    clone_prop = [clone_p * (1-prop_hc) for clone_p in raw_clone_prop]
    clone_prop.append(1-sum(clone_prop))
    
    # Bring kwargs to args as starmap accepts only positional arguments
    extra_args = tuple(kwargs.values())
    
    with mp.Pool(processes=threads) as pool:
        pool.starmap(
            targetedSim_bulk,
            [(i+1, clone_prop, thread_cov, num_clones, *extra_args) for i in range(threads)]
        )
    
    return raw_clone_prop

def targetedSim_sc_parallel(num_single_cells, prop_hc, coverage, num_clones, alpha=None, clone_prop=None, r=None, p=None, threads=None, **kwargs):
    # Simulate negative binomial coverage
    if(any(x is None for x in [p, r])):
            print("Negative binomial parameters for single-cell depths are needed in single-cell mode.")
    cell_cov = np.random.negative_binomial(r, p, size=num_single_cells) # Sample cell depths from a negative binomial
    cell_cov = cell_cov / cell_cov.sum() * coverage # Scale to match total coverage
    if(threads == None):
        print("Using 4 cores as number of cores to use is not specified.")
        threads = 4
    
    # Compute healthy and tumor cells to simulate
    hc = int(num_single_cells * prop_hc // 1)
    tc = num_single_cells - hc
    
    # Compute tumor clones proportions
    if clone_prop is None:
        clone_prop = getDirichletClone(num_clones, alpha)
        
    # Bring kwargs to args as starmap accepts only positional arguments
    extra_args = tuple(kwargs.values())
    
    with mp.Pool(processes=threads) as pool:
        pool.starmap(
            targetedSim_sc,
            [(i+1, num_clones, cell_cov[i], *extra_args) for i in range(hc)]
        )
    with mp.Pool(processes=threads) as pool:
        pool.starmap(
            targetedSim_sc,
            [(i+1+hc, pickdclone(clone_prop, num_clones), cell_cov[i+hc], *extra_args) for i in range(tc)]
        )
        
    return clone_prop

def aggregate_fastqs(fastq_dir, output_left_fastq, output_right_fastq=None, paired=False):
    left_fastqs = sorted(glob.glob(f"{fastq_dir}/*/*left.fq.gz"))
    if(paired):
        right_fastqs = sorted(glob.glob(f"{fastq_dir}/*/*right.fq.gz"))

    with open(output_left_fastq, "ab") as out_f:
        for fq in left_fastqs:
            with open(fq, "rb") as in_f:
                shutil.copyfileobj(in_f, out_f)
    if(paired):
        with open(output_right_fastq, "ab") as out_f:
            for fq in right_fastqs:
                with open(fq, "rb") as in_f:
                    shutil.copyfileobj(in_f, out_f)

    cleanup_dir(fastq_dir)
    return(0)

def cleanup_dir(fastq_dir):
    for path in glob.glob(os.path.join(fastq_dir, "*")):
        if os.path.isdir(path):
            shutil.rmtree(path)
