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
import random
import re
import shutil
import statistics
import string
import tskit

from create_mutations import *
from apply_mutations import *

def getmemory():
    process = psutil.Process(os.getpid())
    print(process.memory_info().rss)

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


def getTree(num_clones, pop, working_dir, seq_len=1):
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
                     bfb_rate, chromothripsis_rate, inssmall_rate, kat_rate, an_rate] ##### ADD chromoplexy_rate
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

def saveMutations(current_genome, tot_nodes, working_dir, list_of_paths, use_signatures, mutationedge_list, num_signatures, signature_alpha, signature_distributions, signatures_matrix, numchrommap, list_of_bases, list_of_pairs, tab):
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
        'CHROMOTHRIP': create_chromothripsis,
        # 'CHROMOPLEX': create_chromoplexy,
        'INSERTIONSMALL': apply_insertion,
        'KATAEGIS': apply_kataegis,
        'ANEUPLOIDY': apply_aneuploidy
    }
    infos = {}
    muts = {}
    infos[tot_nodes - 1] = []
    muts[tot_nodes - 1] = []
    for path in list_of_paths:
        mutated_genome = copy.deepcopy(current_genome) # Create a copy not to modify directly the genome
        for i in range(len(path)-1):
            if(path[i+1] != (tot_nodes-1)):
                c_infos = infos[path[i]].copy()
                c_muts = muts[path[i]].copy()
                current_muts = mutationedge_list[path[i+1]]
                
                for m_type in current_muts.keys():
                    for m in current_muts[m_type]:
                        if(m_type == "SNV" and use_signatures):
                            info = save_functions[m_type](mutated_genome, num_signatures, signature_alpha, signature_distributions, signatures_matrix, numchrommap, list_of_bases, list_of_pairs, tab)
                        else:
                            info = save_functions[m_type](mutated_genome, numchrommap)
                        if info is None:
                            print("invalid")
                        c_infos.append(info)
                        c_muts.append(m_type)
                        # Adjust mutations info for the mutated genomes (shift indices if there are two or more events on the same chromosome)
                        mutated_genome = apply_functions[m_type](mutated_genome, info)
                infos[path[i+1]] = c_infos
                muts[path[i+1]] = c_muts
    return infos, muts

def applyMutations(seqs, infos, muts, clone):
    '''
    Creates a new genome for the given clone.
    '''
    mut = muts[clone]
    info = infos[clone]
    functions = {
        'SNV': apply_SNP,
        'CNV': apply_CNV,
        'DEL': apply_deletion,
        'DELSMALL': apply_deletion,
        'INVERSION': apply_inversion,
        'TRANSLOCATION': apply_translocation,
        'BFB': apply_BFB,
        'CHROMOTHRIP': create_chromothripsis,
        # 'CHROMOPLEX': create_chromoplexy,
        'INSERTIONSMALL': apply_insertion,
        'KATAEGIS': apply_kataegis,
        'ANEUPLOIDY': apply_aneuploidy
    }
    
    clone_seqs = copy.deepcopy(seqs)
    for i, m in enumerate(mut):
        mut_info = info[i]
        clone_seqs = functions[m](clone_seqs, mut_info)
        
    return clone_seqs

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

def wgsSim(ls, num_clones, coverage, rl, fl, floc, batch, alpha, erate, tab, infos, muts, r=None, p=None, flag=0, num_single_cells = 1, paired=False):
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
            chroms = applyMutations(ls, infos, muts, clone)
        while(cov < target_cov[i]):
            if(flag == 1): # Pick a clone for every added coverage if it is a bulk simulation
                distn = getDirichletClone(num_clones, alpha)
                clone = pickdclone(distn, num_clones)
                chroms = applyMutations(ls, infos, muts, clone)
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
                        qual1 = 'K'*len(sub[:rl])
                        pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
                        if(paired):
                            qual2 = 'K'*len(sub[-rl:])
                            pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
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

def targetedSim_bulk(ls, num_clones, coverage, rl, fl, floc, batch, exonDict, numchrommap, alpha, erate, tab, infos, muts, paired=False):
    os.makedirs(floc, exist_ok=True)
    # Open files to write
    f1 = io.BufferedWriter(gzip.open(os.path.join(floc, 'bulkleft.fq.gz'), 'wb'), buffer_size = 4 * 1024**2)
    if(paired):
        f2 = io.BufferedWriter(gzip.open(os.path.join(floc, 'bulkright.fq.gz'), 'wb'), buffer_size = 4 * 1024**2)
    else:
        fl = rl
    # Compute panel length to scale coverage and read quality
    panel_size = float(sum((interval[1] - interval[0]) for cnum in range(25) for interval in exonDict[numchrommap[cnum]]))
    qual = 'K'*rl
    # Initialize coverage
    cov = 0.0
    while(cov < coverage):
        # Pick a clone for every added coverage in a bulk simulation
        distn = getDirichletClone(num_clones, alpha)
        clone = pickdclone(distn, num_clones)
        chroms = applyMutations(ls, infos, muts, clone)
        chromnum = 0
        for j in range(batch):
            for chrom in chroms:
                if(len(chrom)==0): # Deleted chromosomes
                    continue
                # frag_len = 0
                # while(frag_len <= rl):
                #     frag_len = getfrag(fl)
                # Sample intervals consistent with the coverage of each cell (for very low coverage)
                k = math.ceil(len(exonDict[numchrommap[chromnum]]) * min(coverage, 1))
                sampled_intervals = random.sample(exonDict[numchrommap[chromnum]], k)
                for interval in sampled_intervals:
                    # Update coverage after each iteration
                    cov += (2*batch*rl / panel_size) if paired else (batch*rl / panel_size)
                    startindex = random.randint(interval[0], interval[1])
                    sub = chrom[startindex:startindex+fl]
                    if random.random() > 0.5:
                        sub = revc(sub, tab)
                    random_str = ''.join(random.choices(
                        string.ascii_letters, k=15))
                    pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
                    if(paired):
                        pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
                    f1.write(('\n'.join([f'@{random_str}', pair1, '+', qual]) + '\n').encode('utf-8'))
                    if(paired):
                        f2.write(('\n'.join([f'@{random_str}', pair2, '+', qual]) + '\n').encode('utf-8'))
                chromnum += 1
    # Close files
    f1.close()
    if(paired):
        f2.close()
    return(0)


def targetedSim_sc(cell_id, target_cov, ls, num_clones, rl, fl, floc, batch, exonDict, numchrommap, alpha, erate, tab, infos, muts, paired=False):
    # per-cell output directory or prefix
    cell_dir = os.path.join(floc, f"cell_{cell_id}")
    os.makedirs(cell_dir, exist_ok=True)
    f1 = io.BufferedWriter(gzip.open(os.path.join(cell_dir, 'left.fq.gz'), 'wb'), buffer_size=4 * 1024**2)
    if paired:
        f2 = io.BufferedWriter(gzip.open(os.path.join(cell_dir, 'right.fq.gz'), 'wb'), buffer_size=4 * 1024**2)
    
    # Compute panel length to scale coverage and read quality
    panel_size = float(sum((interval[1] - interval[0]) for cnum in range(25) for interval in exonDict[numchrommap[cnum]]))
    qual = 'K' * rl
    cov = 0.0

    distn = getDirichletClone(num_clones, alpha)
    clone = pickdclone(distn, num_clones)
    chroms = applyMutations(ls, infos, muts, clone)

    while cov < target_cov:
        cell_cov = 0.0
        chromnum = 0

        for j in range(batch):
            for chrom in chroms:
                if len(chrom) == 0:
                    continue
                k = math.ceil(len(exonDict[numchrommap[chromnum]]) * min(target_cov, 1))
                sampled_intervals = random.sample(exonDict[numchrommap[chromnum]], k)
                for interval in sampled_intervals:
                    cell_cov += (2 * batch * rl) if paired else (batch * rl)
                    startindex = random.randint(interval[0], interval[1])
                    sub = chrom[startindex:startindex + fl]
                    if random.random() > 0.5:
                        sub = revc(sub, tab)
                    random_str = ''.join(random.choices(string.ascii_letters, k=15))
                    pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
                    if(paired):
                        pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
                    f1.write(('\n'.join([f'@cell{cell_id}_clone{clone}_{random_str}', pair1, '+', qual]) + '\n').encode('utf-8'))
                    if(paired):
                        f2.write(('\n'.join([f'@cell{cell_id}_clone{clone}_{random_str}', pair2, '+', qual]) + '\n').encode('utf-8'))
                chromnum += 1
        cov += cell_cov / panel_size
    f1.close()
    if paired:
        f2.close()
    return(0)

def targetedSim_sc_parallel(num_single_cells, coverage, r=None, p=None, threads=None, *args):
    # Simulate negative binomial coverage
    if(any(x is None for x in [p, r])):
            print("Negative binomial parameters for single-cell depths are needed in single-cell mode.")
    cell_cov = np.random.negative_binomial(r, p, size=num_single_cells) # Sample cell depths from a negative binomial
    cell_cov = cell_cov / cell_cov.sum() * coverage # Scale to match total coverage
    if(threads == None):
        print("Using 4 cores as number of cores to use is not specified.")
        threads = 4
    with mp.Pool(processes=threads) as pool:
        pool.starmap(
            targetedSim_sc,
            [(i+1, cell_cov[i], *args) for i in range(num_single_cells)]
        )

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

def lbrunSim(num_tumors, num_clones_list, coverage, base_dir, floc, root, alpha, ctdna_frac, batchsize):
    f = open(floc + 'liquid_biopsyfull.fasta', 'w')
    ls = []
    ls = rgz(f'{base_dir}reference/{root}.gz')
    cfdna_frag = []
    ctdna_frag = []
    tot_clones = 0
    for i in num_clones_list:
        tot_clones += i
    distn = getDirichletClone(tot_clones, alpha)
    frags_per_clone = int(10e7/tot_clones)
    for i in range(tot_clones):
        tnum = random.randint(0, num_tumors-1)
        num_clones = num_clones_list[tnum]
        clone = pickdclone(distn, num_clones)
        thechrom = rgz(base_dir + f'tumor_{tnum}/{clone}.gz')
        for j in range(frags_per_clone):
            l = getfrag(133)
            ctdna_frag.append(getDNAchunk(l, thechrom))
    thechrom = rgz(base_dir + f'reference/{root}.gz')
    for i in range(int(9e7)):
        l = getfrag(166)
        cfdna_frag.append(getDNAchunk(l, thechrom))
    total_frags_needed = int(3e9/150)*coverage
    frags_from_tumor = int(ctdna_frac*total_frags_needed)
    frags_from_normal = total_frags_needed - frags_from_tumor
    perbatch_tumor = int(frags_from_tumor/batchsize)
    perbatch_normal = int(frags_from_normal/batchsize)
    for i in range(batchsize):
        print(i)
        alltumorfrags = random.choices(ctdna_frag, k=perbatch_tumor)
        allnormalfrags = random.choices(cfdna_frag, k=perbatch_normal)
        for frag in alltumorfrags:
            random_str = ''.join(random.choices(string.ascii_letters, k=15))
            f.write(f'@{random_str}\n')
            f.write(f'{frag}\n')
            f.write('+\n')
            qual = 'K'*len(frag)
            f.write(f'{qual}\n')
        for frag in allnormalfrags:
            random_str = ''.join(random.choices(string.ascii_letters, k=15))
            f.write(f'@{random_str}\n')
            f.write(f'{frag}\n')
            f.write('+\n')
            qual = 'K'*len(frag)
            f.write(f'{qual}\n')
        print('first batch write done')