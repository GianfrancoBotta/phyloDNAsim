import numpy as np
import random
import re

def create_SNPSig(seqs, num_signatures, signature_alpha, signature_distributions, signatures_matrix, numchrommap, list_of_bases, list_of_pairs, tab):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    # seq = seqs[chrom[0]]
    signature_drawn_dist = getDirichletCloneFromDistn(
        num_signatures-1, signature_alpha, signature_distributions)
    signature_drawn = pickdclone(
        signature_drawn_dist, len(signature_drawn_dist)-1)
    distribution_over_mutations = list(signatures_matrix[:, signature_drawn])
    position = pickdclone(distribution_over_mutations,
                          len(distribution_over_mutations)-1)
    pos = 0
    first_letter = int(position/24)
    last_letter = position % 24 % 4  # if == 0 > A if 1 > C if 2 > G if 3 >T
    middle_two = int((position % 24)/4)
    if random.random() < 0.5:
        fullstring = list_of_bases[first_letter] + \
            list_of_pairs[middle_two][0] + list_of_bases[last_letter]
        mutated_base = list_of_pairs[middle_two][1]
        indices = [m.start() for m in re.finditer(fullstring, seqs[chrom[0]])]
        pos = random.choice(indices)+1 # Change the base in the middle, not the one detected
        # seq = seq[:pos] + mutated_base + seq[pos+1:]
    else:
        altstring = list_of_bases[first_letter] + \
            list_of_pairs[middle_two][0].translate(
                tab) + list_of_bases[last_letter]
        mutated_base = list_of_pairs[middle_two][1].translate(tab) # Reverse strand
        indices = [m.start() for m in re.finditer(altstring, seqs[chrom[0]])]
        pos = random.choice(indices)+1
        # seq = seq[:pos] + mutated_base + seq[pos+1:]
    # seqs[chrom[0]] = seq
    return {'char': mutated_base, 'pos': pos, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}

def create_speedSNP(seqs, numchrommap, bases=['A', 'C', 'T', 'G']):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    pos = random.randint(0, n-1)
    # original_char = seqs[chrom[0]][pos]
    mut_bases = list(set(bases) - set(list(chr(seqs[chrom[0]][pos]))))
    # mut_bases = [c for c in bases if c not in list(original_char)]
    char = random.choice(mut_bases)
    return {'char': char, 'pos': pos, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}

def create_insertionsmall(seqs, numchrommap, bases=['A', 'C', 'T', 'G']):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    pos = random.randint(0, n-1)
    number_of_bases_to_insert = random.randint(1, 15)
    ins_str = ''
    for i in range(number_of_bases_to_insert):
        ins_str = ins_str + random.choice(bases)
    # new_seq = seq[:pos] + ins_str + seq[pos:]
    # seqs[chrom[0]] = new_seq
    return {'insertion': ins_str, 'pos': pos, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}

def create_deletionsmall(seqs, numchrommap):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    # pos = 0
    # if (n == 0):
    #     pos = -1
    # elif (n == 1):
    #     pos = 0
    # else:
    stidx = random.randint(0, n-1)
    edidx = stidx + random.randint(1, 15)
    # new_seq = seq[:pos] + seq[pos+deletion_range:]
    # seqs[chrom[0]] = new_seq
    # position = int(pos)
    # chromosome = int(chrom[0])
    return {'start': stidx, 'end': edidx, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}

def create_CNV(seqs, numchrommap):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    # seq = seqs[chrom[0]]
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    dist = n
    while(dist > 0.1*n):
        dist = getSVSize(chromosome_size=n) # Insert a CN which is shorter than 1/10 the length of the chromosome
    stidx = random.randint(0, n-1)
    edidx = min(stidx + dist, n-1)
    max_reps = 10 # It was 5
    rep_num = random.randint(2, max_reps)
    # subseq = seq[stidx:edidx]
    # repeatedsubseq = subseq*rep_num
    # start_string = seq[:stidx]
    # end_string = seq[edidx:]
    # full_string = start_string + repeatedsubseq + end_string
    # seqs[chrom[0]] = full_string
    return {'start': stidx, 'end': edidx, 'rep_num': rep_num, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}


def create_aneuploidy(seqs, numchrommap):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    max_reps = 3
    if (random.random() < 0.2): # Deletion with probabilty 0.2
        rep_num = 0
    else: # Insertion with probability 0.8
        rep_num = random.randint(2, max_reps)
    return {'rep_num': rep_num, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}


def create_deletion(seqs, numchrommap):
    '''
    For long deletions.
    '''
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    # seq = seqs[chrom[0]]
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    # if(n == 0):
    #     return chrom[0]
    # elif(n == 1):
    #     seqs[chrom[0]] = ''
    #     return [0, chrom[0], target]
    dist = n
    while(dist > 0.1*n): # Delete a fragent which is shorter than 1/10 the length of the chromosome
        dist = getSVSize() 
    stidx = random.randint(0, n-1)
    edidx = min(stidx + dist, n-1)
    # start_string = seq[:stidx]
    # end_string = seq[edidx:]
    # full_string = start_string + end_string
    # seqs[chrom[0]] = full_string
    return {'start': stidx, 'end': edidx, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}


def create_inversion(seqs, numchrommap):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    # seq = seqs[chrom[0]]
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    # if (n == 0 or n == 1):
    #     return [-1, chrom[0]]
    dist = n
    while(dist > 0.1*n): # Revert a fragent which is shorter than 1/10 the length of the chromosome
        dist = getSVSize() 
    stidx = random.randint(0, n-1)
    edidx = min(stidx + dist, n-1)
    # subseq = seq[stidx:edidx]
    # reversed_subseq = subseq[::-1]
    # start_string = seq[:stidx]
    # end_string = seq[edidx:]
    # full_string = start_string + reversed_subseq + end_string
    # seqs[chrom[0]] = full_string
    return {'start': stidx, 'end': edidx, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}


def create_kataegis(seqs, numchrommap):
    '''
    Simulate regions of the genomes that are hypermutated.
    '''
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    # if (n == 0 or n == 1):
    #     return [-1, chrom[0]]
    dist = n
    while(dist > 0.1*n):
        dist = getSVSize() 
    stidx = random.randint(0, n-1)
    edidx = min(stidx + dist, n-1)
    type_of_variation = random.randint(0, 2) # Consider picking 1 less frequently
    # subseq = seq[stidx:edidx]
    # new_seq = ''
    # if(type_of_variation == 0):
    #     def rand_replace(c): return 'T' if random.random(
    #     ) < 0.8 and c == 'C' else c
    #     new_seq = ''.join([rand_replace(c) for c in subseq])
    # elif(type_of_variation == 1):
    #     def rand_replace(c): return 'A' if random.random(
    #     ) < 0.8 and c == 'C' else c
    #     new_seq = ''.join([rand_replace(c) for c in subseq])
    # else:
    #     def rand_replace(c): return 'G' if random.random(
    #     ) < 0.8 and c == 'C' else c
    #     new_seq = ''.join([rand_replace(c) for c in subseq])
    # start_string = seq[:stidx]
    # end_string = seq[edidx:]
    # full_string = start_string + new_seq + end_string
    
    # Pick mutational type
    # seqs[chrom[0]] = full_string
    idx = stidx
    res = []
    for c in seqs[chrom[0]][stidx:edidx]:
        if random.random() < 0.8 and c == ord('C'):
            if(type_of_variation == 0):
                mut = 'T'
            elif(type_of_variation == 1): 
                mut = 'A'
            else:
                mut = 'G'
            res.append({'chrom_num': chrom[0], 'pos': idx, 'char': mut})
        idx = idx + 1
    return {'start': stidx, 'end': edidx, 'mut': res, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}

def create_translocation(seqs, numchrommap):
    c = list(range(len(seqs)))
    chrom = random.sample(c, 2)
    while(len(seqs[chrom[0]]) == 0 or len(seqs[chrom[1]]) == 0): # Sample until you find chromosomes that were not deleted before
        chrom = random.sample(c, 2)
    l1 = len(seqs[chrom[0]])
    l2 = len(seqs[chrom[1]])
    bkpt1 = random.randint(0, l1-1)
    bkpt2 = random.randint(0, l1-1)
    rng_recip = random.random()
    if(rng_recip < 0.95):
        # first_part = seq1[:bkpt1]
        # second_part = seq1[bkpt1:]
        # fp2 = seq2[:bkpt2]
        # sp2 = seq2[bkpt2:]
        # string1 = first_part + sp2
        # string2 = fp2 + second_part
        # seqs[chrom[0]] = string1
        # seqs[chrom[1]] = string2
        res = {'chrom_num1': chrom[0], 'chrom_name1': numchrommap[chrom[0]], 'chrom_num2': chrom[1], 'chrom_name2': numchrommap[chrom[1]], 'bkpt1': bkpt1, 'bkpt2': bkpt2, 'normal': True}
    else:
        # first_part = seq1[:bkpt1]
        # second_part = seq1[bkpt1:]
        # fp2 = seq2[:bkpt2]
        # sp2 = seq2[bkpt2:]
        # string1 = first_part
        # string2 = fp2 + sp2 + second_part
        # seqs[chrom[0]] = string1
        # seqs[chrom[1]] = string2
        res = {'chrom_num1': chrom[0], 'chrom_name1': numchrommap[chrom[0]], 'chrom_num2': chrom[1], 'chrom_name2': numchrommap[chrom[1]], 'bkpt1': bkpt1, 'bkpt2': bkpt2, 'normal': False}
    return res

def create_chromothripsis(seqs, numchrommap):
    keep_frequency = random.random()
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    # if (n == 0 or n == 1):
    #     return [-1], [-1, chrom[0]]
    NotLongEnough = True
    attempts = 0
    stidx = 0
    edidx = 0
    splits = 0
    while(NotLongEnough):
        splits = random.randint(2, 100)
        dist = n
        while(dist > 0.1*n):
            dist = getSVSize()
        stidx = random.randint(0, n-1)
        edidx = min(stidx + dist, n-1)
        dist = edidx - stidx
        attempts += 1
        if(dist > splits):
            NotLongEnough = False
        # if(attempts > 10):
        #     return [-1], [-1, chrom[0]]
    # first_part = orig_seq[:stidx]
    # last_part = orig_seq[edidx:]
    # seq = orig_seq[stidx:edidx]
    # seq_len = len(seq)
    breakpoints = numpy_choices(len(seqs[chrom[0]])[stidx:edidx], splits) # Returns locations of breakpoints
    # subseq = []
    # curridx = 0
    # for i in breakpoints:
    #     subseq.append(seq[curridx:i])
    #     curridx = i
    # subseq.append(seq[breakpoints[-1]:])
    rearrange = np.random.permutation(splits).tolist()
    n_to_select = int(splits * keep_frequency)
    rearrange = random.sample(list(rearrange), n_to_select)
    # build_seq = ''
    # for i in rearrange:
    #     build_seq += subseq[i]
    # seqs[chrom[0]] = first_part + build_seq + last_part
    # breakpoints = list(breakpoints)
    return {'start': stidx, 'end': edidx, 'bkpts': list(breakpoints), 'order': rearrange, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}

def create_BFB(seqs, numchrommap):
    n_bkpts = random.randint(2, 4)
    c = list(range(len(seqs)))
    chrom = random.sample(c, 1)
    while(len(seqs[chrom[0]]) == 0):
        chrom = random.sample(c, 1)
    n = len(seqs[chrom[0]])
    # seq = seqs[chrom[0]]
    # curr_seq = seq
    bkpts = []
    for i in range(n_bkpts):
        bkpoint = random.randint(1, len(seqs[chrom[0]])-2)
        bkpts.append(bkpoint)
        # first_portion = curr_seq[:bkpoint]
        # next_portion = first_portion[::-1]
        # curr_seq = first_portion + next_portion
    # seqs[chrom[0]] = curr_seq
    return {'bkpts': bkpts, 'chrom_num': chrom[0], 'chrom_name': numchrommap[chrom[0]]}

def create_chromoplexy(seqs, numchrommap):
    n_chroms = random.randint(2, 6)
    c = list(range(len(seqs)))
    chrom = random.sample(c, n_chroms)
    while any(len(seqs[ch]) == 0 for ch in chrom):
        chrom = random.sample(c, n_chroms)
    n = len(seqs[chrom[0]])
    list_of_start_tels = []
    list_of_end_tels = []
    middle_segs = []
    total_splits = 0
    for i in range(n_chroms):
        num_splits = random.randint(2, 20)
        breakpoints = numpy_choices(len(seqs[chrom[i]]), num_splits) # Returns locations of breakpoints
        total_splits += num_splits
        curridx = 0
        # for i in breakpoints[1:]:
        #     middle_segs.append(curr_sequence[curridx:i])
        #     curridx = i
        # list_of_end_tels.append(curr_sequence[breakpoints[-1]:])
        # list_of_start_tels.append(curr_sequence[0:breakpoints[0]])
    usable_splits = total_splits - 2*(n_chroms)
    splitter = getsummation(usable_splits, n_chroms) # Divides the usable splits in the chromosomes
    # for i in range(len(splitter)):
    #     build_str = list_of_start_tels.pop(
    #         random.randrange(len(list_of_start_tels)))
    #     for segs in range(splitter[i]):
    #         build_str = build_str + \
    #             middle_segs.pop(random.randrange(len(middle_segs)))
    #     build_str = build_str + \
    #         list_of_end_tels.pop(random.randrange(len(list_of_end_tels)))
    #     seqs[chrom[i]] = build_str
    return {'splits': splitter, 'chrom_num': chrom, 'chrom_name': numchrommap[chrom]}

# Functions to handle chromosomes and clones' proportions    
    
def getSVSize(chromosome_size = 1e8):
    '''
    Generates the size of a chromosome with a structural variant, returning minimum 50bp, and maximum 2e7bp
    '''
    return_val1 = np.random.negative_binomial(0.1, 0.0001) # extremely variable, small mean
    return_val2 = np.random.negative_binomial(1, 0.00001) # larger mean, more typical SV sizes
    return_val3 = np.random.negative_binomial(200, 0.00002) # very large mean, rare but huge events
    z = random.random()
    if(z < 0.5): 
        o_return_val = return_val1
    elif(z > 0.5 and z < 0.995):
        o_return_val = return_val2
    else: 
        o_return_val = return_val3
    return(min(2e7, max(50, int((chromosome_size/1e8)*o_return_val))))

def getDirichletCloneFromDistn(num_clones, alpha_dir, the_distn):
    cumulate_product = 1.0
    current_weights = []
    ep = 0.001
    epsil = 1-ep
    tot = 0.0
    while(tot < epsil):
        beta = random.betavariate(1, alpha_dir)
        pi_app = beta*cumulate_product
        cumulate_product = cumulate_product * (1-beta)
        current_weights.append(pi_app)
        tot += pi_app
    sum_current_weights = sum(current_weights)
    weights = []
    for i in current_weights:
        weights.append(i/sum_current_weights)
    print(weights)
    distn = [0.0]*(num_clones+1)
    for i in weights:
        all_choices = list(range(num_clones+1))
        pos = random.choices(all_choices, weights=the_distn)[0]
        distn[pos] = distn[pos] + i
    return distn

def getKmers(sequence, size):
    return {sequence[x:x+size].upper() for x in range(len(sequence) - size + 1)}

def getDirichletClone(num_clones, alpha_dir):
    cumulate_product = 1.0
    current_weights = []
    ep = 0.001
    epsil = 1-ep
    tot = 0.0
    while(tot < epsil):
        beta = random.betavariate(1, alpha_dir)
        pi_app = beta*cumulate_product
        cumulate_product = cumulate_product * (1-beta)
        current_weights.append(pi_app)
        tot += pi_app
    sum_current_weights = sum(current_weights)
    weights = []
    for i in current_weights:
        weights.append(i/sum_current_weights)
    distn = [0.0]*(num_clones+1)
    for i in weights:
        pos = random.randint(0, num_clones)
        distn[pos] = distn[pos] + i
    return distn

def pickdclone(prob_list, num_clones):
    return np.random.choice(np.arange(num_clones+1), p=prob_list)

def numpy_choices(seq_len, splits):
    breakpoints = np.random.choice(seq_len, splits-1, replace=False)
    breakpoints.sort()
    return breakpoints.tolist()

def getsummation(seq_len, splits):
    breakpoints = np.random.choice(seq_len, splits-1, replace=False)
    breakpoints.sort()
    zed = breakpoints.tolist()
    new_list = []
    idx = 0
    for i in zed:
        new_list.append(i-idx)
        idx = i
    new_list.append(seq_len - idx)
    return new_list