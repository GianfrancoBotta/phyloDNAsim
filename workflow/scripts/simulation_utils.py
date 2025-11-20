import copy
import gc
import glob
import gzip
import json
import msprime
import numpy as np
import os
import psutil
import random
import re
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


def getTree(num_clones, pop, working_dir):
    tree_sequence = msprime.simulate(
        sample_size=num_clones, Ne=pop, recombination_rate=0)
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


# def getTimeMatrix(tree, num_clones):
#     time_matrix = np.zeros((tree.root + 1, tree.root+1))
#     dep = 0
#     for u in range(num_clones):
#         while u != tskit.NULL:
#             time_matrix[u, tree.parent(u)] = tree.get_branch_length(u)
#             dep += tree.get_branch_length(u)
#             u = tree.parent(u)
#     time_matrix = np.transpose(time_matrix)
#     return time_matrix, dep


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
        # gives mutations from i to its parent ONLY
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
        merged_list = sorted(sum(mutation_lists, []))
        ordered_muts = {}
        for idx, lst in enumerate(mutation_lists):
            ordered_muts.update({list(list_of_rates.keys())[idx]: event for event in lst})
        # ordered_muts = [lookup[list(list_of_rates.keys())[idx]] for event in merged_list]
        mutationedge_list.append(ordered_muts)
    return mutationedge_list, fin_rate_list

def wgz(floc, ob):
    with gzip.open(floc, 'wt') as f:
        for i in ob:
            f.write(f'{i}\n')

def rgz(floc):
    return_list = []
    with gzip.open(floc, 'rt') as file:
        return_list = [line.strip() for line in file]
    return return_list

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
    d = {}
    infos = {}
    muts = {}
    for i in range(tot_nodes):
        d[i] = 'na'
    d[tot_nodes-1] = os.path.join(working_dir, 'node' + str(tot_nodes - 1) + '_mutations.json') # Save mutations to json file to load and apply to the genome
    infos[tot_nodes - 1] = []
    muts[tot_nodes - 1] = []
    for path in list_of_paths:
        for i in range(len(path)-1):
            if(d[path[i+1]] == 'na'):
                c_infos = infos[path[i]].copy()
                c_muts = muts[path[i]].copy()
                current_muts = mutationedge_list[path[i+1]]
                
                for m in current_muts:
                    if(m == "SNV" and use_signatures):
                        info = save_functions[m](current_genome, num_signatures, signature_alpha, signature_distributions, signatures_matrix, numchrommap, list_of_bases, list_of_pairs, tab)
                    else:
                        info = save_functions[m](current_genome, numchrommap)
                    if info is None:
                        print("invalid")
                    c_infos.append(info)
                    c_muts.append(m)
                    # Adjust mutations info for the mutated genomes (shift indices if there are two or more events on the same chromosome)
                    current_genome = apply_functions[m](current_genome, info)
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

# def clear_dir(working_dir):
#     files = glob.glob(working_dir+'*')
#     for f in files:
#         os.remove(f)


# def pickClone(num_clones, dir_conc):
#     alpha_dir = np.tile(dir_conc, num_clones)
#     distribution = np.random.dirichlet(alpha_dir)
#     clone_picked = np.random.choice(range(0, num_clones), p=distribution)
#     return int(clone_picked)


# def writeToFasta(file_loc, curr_string):
#     import random
#     with open(file_loc, 'at') as f:
#         f.write(curr_string)


def getfrag(r):
    return_val = 0
    while(return_val < r-20):
        return_val = np.random.negative_binomial(r, 0.5)
    return return_val


def drawPoisson(r):
    return_val = 0
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
    # s = str(frag)
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


def runSim(num_clones, coverage, rl, read_loc, floc, batch, root, alpha, erate, tab, flag=0):
    if(flag == 0):
        f = open(floc + 'bulk.fasta', 'w')
    elif(flag == 1):
        f = open(floc + 'ref.fasta', 'w')
    else:
        randid = random_str = ''.join(random.choice(
            string.ascii_lowercase) for _ in range(3))
        f = open(floc + f'singlecell.fasta', 'w')
    cov = 0.0
    ratio = rl
    #random_str = 0
    ls = []
    if(flag == 2):
        ri = random.randint(0, num_clones)
        ls = rgz(f'{read_loc}{ri}.gz')
    if(flag == 1):
        ls = rgz(f'{read_loc}{root}.gz')
    giga_list = []
    while(cov < coverage):
        print(giga_list)
        print(cov)
        for i in range(batch):
            distn = getDirichletClone(num_clones, alpha)
            clone = pickdclone(distn, num_clones)
            if(flag == 0):
                ls = rgz(f'{read_loc}{clone}.gz')
            for chrom in ls:
                frag_len = 0
                while(frag_len <= rl):
                    frag_len = getfrag(rl)
                if random.random() < 0.5:
                    altchrom = chrom
                else:
                    altchrom = revc(chrom, tab)
                for sub in split(altchrom, frag_len):
                    random_str = ''.join(random.choices(
                        string.ascii_letters, k=15))
                    sub = mutateFrag(sub, erate)  # check type here
                    qual = 'K'*len(sub)
                    giga_list.extend([f'@{random_str}', str(sub), '+', qual])
        for x in giga_list:
            f.write(x)
            f.write('\n')
        giga_list.clear()
        print('first batch write done')
        # f.write('\n'.join(giga_list))
        cov += 2*batch
        #del giga_list
    del giga_list
    del ls

def runPairedSim(ls, num_clones, coverage, rl, fl, read_loc, floc, batch, root, alpha, erate, tab, flag=0):
    if(flag == 0):
        f1 = gzip.open(os.path.join(floc, 'bulkleft.fq.gz'), 'wt')
        f2 = gzip.open(os.path.join(floc, 'bulkright.fq.gz'), 'wt')
    elif(flag == 1):
        f1 = gzip.open(os.path.join(floc, 'refleft.fq.gz'), 'wt')
        f2 = gzip.open(os.path.join(floc, 'refright.fq.gz'), 'wt')
    else:
        randid = random_str = ''.join(
        random.choices(string.ascii_lowercase, k=4))
        f1 = gzip.open(os.path.join(floc, 'singlecellleft.fq.gz'), 'wt')
        f2 = gzip.open(os.path.join(floc, 'singlecellright.fq.gz'), 'wt')
    if(flag == 2):
        ri = random.randint(0, num_clones)
        ls = rgz(f'{read_loc}{ri}.gz')
    cov = 0.0
    ratio = 2*rl/fl
    while(cov < coverage):
        for i in range(batch):
            distn = getDirichletClone(num_clones, alpha)
            clone = pickdclone(distn, num_clones)
            if(flag == 0):
                ls = rgz(f'{read_loc}{clone}.gz')
            for chrom in ls:
                frag_len = 0
                while (frag_len <= rl):
                    frag_len = getfrag(fl)
                if random.random() < 0.5:
                    altchrom = chrom
                else:
                    altchrom = revc(chrom, tab)
                for sub in split(altchrom, frag_len):
                    random_str = ''.join(random.choices(
                        string.ascii_letters, k=15))
                    #sub = mutateFrag(sub, erate)
                    pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
                    pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
                    qual1 = 'K'*len(pair1)
                    qual2 = 'K'*len(pair2)
                    f1.write('\n'.join([f'@{random_str}', pair1, '+', qual1]) + '\n')
                    f2.write('\n'.join([f'@{random_str}', pair2, '+', qual2]) + '\n')
                    # giga_list.extend([f'@{random_str}', write1, '+', qual])
                    # giga_list2.extend(
                    #     [f'@{random_str}', write2, '+', qual2])
        # for x in giga_list:
        #     f.write(x)
        #     f.write('\n')
        # for x in giga_list2:
        #     f2.write(x)
        #     f2.write('\n')
        # giga_list.clear()
        # giga_list2.clear()
        # print('first batch write done')
        cov += 2*batch*ratio
    # del giga_list
    # del giga_list2
    # del ls
    return(0)

def get_random_str(main_str, substr_len):
    if(len(main_str) < substr_len): 
        return main_str, 0
    else: 
        idx = random.randrange(0, len(main_str) - substr_len + 1)
        return main_str[idx : (idx+substr_len)], idx


def exonrunSim(num_clones, coverage, rl, rloc, floc, batch, root, exonDict, numchrommap, subblock, alpha, erate, tab, flag=0):
    if(flag == 0):
        f = open(floc + 'bulk.fasta', 'w')
    elif(flag == 1):
        f = open(floc + 'ref.fasta', 'w')
    else:
        randid = random_str = ''.join(random.choice(
            string.ascii_lowercase) for _ in range(3))
        f = open(floc + f'singlecell.fasta', 'w')
    cov = 0.0
    ratio = rl
    #random_str = 0
    ls = []
    if(flag == 2):
        ri = random.randint(0, num_clones)
        ls = rgz(f'{rloc}{ri}.gz')
    if(flag == 1):
        ls = rgz(f'{rloc}{root}.gz')
    giga_list = []
    while(cov < coverage):
        print(giga_list)
        print(cov)
        for i in range(batch):
            distn = getDirichletClone(num_clones, alpha)
            clone = pickdclone(distn, num_clones)
            if(flag == 0):
                ls = rgz(f'{rloc}{clone}.gz')
            chromnum = 0
            for chrom in ls:
                frag_len = 0
                while(frag_len <= rl):
                    frag_len = getfrag(rl)
                if random.random() < 0.5:
                    altchrom = chrom
                else:
                    altchrom = revc(chrom, tab)
                for interval in exonDict[numchrommap[chromnum]]:
                    for i in range(subblock):
                        startindex = random.randint(interval[0], interval[1])
                        sub = altchrom[startindex:startindex+rl]
                        qual = 'K'*len(sub)
                        random_str = ''.join(random.choices(
                            string.ascii_letters, k=15))
                        sub = mutateFrag(sub, erate)
                        giga_list.extend([f'@{random_str}', sub, '+', qual])
                chromnum += 1
        for x in giga_list:
            f.write(x)
            f.write('\n')
        giga_list.clear()
        print('first batch write done')
        cov += 2*batch*subblock
    del giga_list


def exonrunPairedSim(ls, num_clones, coverage, rl, fl, rloc, floc, batch, root, exonDict, numchrommap, subblock, alpha, erate, tab, infos, muts, num_single_cells=1, flag=0):
    os.makedirs(floc, exist_ok=True)
    if(flag == 0):
        f1 = gzip.open(os.path.join(floc, 'bulkleft.fq.gz'), 'wt')
        f2 = gzip.open(os.path.join(floc, 'bulkright.fq.gz'), 'wt')
    elif(flag == 1):
        f1 = gzip.open(os.path.join(floc, 'refleft.fq.gz'), 'wt')
        f2 = gzip.open(os.path.join(floc, 'refright.fq.gz'), 'wt')
    else:
        f1 = gzip.open(os.path.join(floc, 'singlecellleft.fq.gz'), 'at') # Append mode since the function gets called for each single-cell if flag==2
        f2 = gzip.open(os.path.join(floc, 'singlecellright.fq.gz'), 'at')
    cov = 0.0
    ratio = 2*rl/fl
    # if(flag == 1):
    #     ls = rgz(f'{rloc}{root}.gz')
    # giga_list = []
    # giga_list2 = []
    if(flag == 2):
        ri = random.randint(0, num_clones)
        chroms = applyMutations(ls, infos, muts, ri)
    elif(flag == 0):
        distn = getDirichletClone(num_clones, alpha)
        clone = pickdclone(distn, num_clones)
        chroms = applyMutations(ls, infos, muts, clone)
    else:
        chroms = ls
    while(cov < coverage/num_single_cells):
        # print(giga_list)
        # print(giga_list2)
        # print(cov)
        chromnum = 0
        for i in range(batch):
            for chrom in chroms:
                if(len(chrom)==0): # Deleted chromosomes
                    continue
                frag_len = 0
                while(frag_len <= rl):
                    frag_len = getfrag(fl)
                # if random.random() < 0.5:
                #     altchrom = chrom
                # else:
                #     altchrom = revc(chrom, tab)
                for interval in exonDict[numchrommap[chromnum]]:
                    for i in range(subblock):
                        startindex = random.randint(interval[0], interval[1])
                        sub = chrom[startindex:startindex+fl]
                        if random.random() > 0.5:
                            sub = revc(sub, tab)
                        random_str = ''.join(random.choices(
                            string.ascii_letters, k=15))
                        qual1 = 'K'*len(sub[:rl])
                        qual2 = 'K'*len(sub[-rl:])
                        #sub = mutateFrag(sub, erate)
                        pair1 = mutateFrag(sub[:rl], erate).decode("utf-8")
                        pair2 = mutateFrag(revc(sub[-rl:], tab), erate).decode("utf-8")
                        if(flag == 2):
                            f1.write('\n'.join([f'@{ri}_{random_str}', pair1, '+', qual1]) + '\n')
                            f2.write('\n'.join([f'@{ri}_{random_str}', pair2, '+', qual2]) + '\n')
                        else:
                            f1.write('\n'.join([f'@{random_str}', pair1, '+', qual1]) + '\n')
                            f2.write('\n'.join([f'@{random_str}', pair2, '+', qual2]) + '\n')
                        # giga_list.extend(
                        #     [f'@{random_str}', pair1, '+', qual1])
                        # giga_list2.extend(
                        #     [f'@{random_str}', pair2, '+', qual2])
                chromnum += 1
        # for x in giga_list:
        #     f.write(x)
        #     f.write('\n')
        # for x in giga_list2:
        #     f2.write(x)
        #     f2.write('\n')
        # giga_list.clear()
        # giga_list2.clear()
        # print('first batch write done')
        cov += 2*batch*subblock*ratio
    # del giga_list
    # del giga_list2
    f1.close()
    f2.close()
    # del chroms
    # gc.collect()
    return(0)

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