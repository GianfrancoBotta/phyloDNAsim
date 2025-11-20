def apply_SNP(seqs, info):
    # seqs[info['chrom_num']] = seqs[info['chrom_num']][:info['pos']] + ord(info['char']) + seqs[info['chrom_num']][info['pos']+1:]
    seqs[info['chrom_num']][info['pos']] = ord(info['char'])
    return seqs

def apply_insertion(seqs, info):
    # seqs[info['chrom_num']] = seqs[info['chrom_num']][:info['pos']] + info['insertion'].encode() + seqs[info['chrom_num']][info['pos']:]
    seqs[info['chrom_num']][info['pos']:info['pos']] = info['insertion'].encode()
    return seqs

def apply_deletion(seqs, info):
    # seqs[info['chrom_num']] = seqs[info['chrom_num']][:info['start']] + seqs[info['chrom_num']][info['end']:]
    del seqs[info['chrom_num']][info['start']:info['end']]
    return seqs

def apply_CNV(seqs, info):
    seqs[info['chrom_num']][info['end']:info['end']] = seqs[info['chrom_num']][info['start']:info['end']]*(info['rep_num']-1)
    return seqs

def apply_aneuploidy(seqs, info):
    # seqs[info['chrom_num']] = seqs[info['chrom_num']] * info['rep_num']
    seqs[info['chrom_num']][:] = seqs[info['chrom_num']] * info['rep_num']
    return seqs

def apply_inversion(seqs, info):
    # seqs[info['chrom_num']] = seqs[info['chrom_num']][:info['start']] + seqs[info['chrom_num']][info['start']:info['end']][::-1] + seqs[info['chrom_num']][info['end']:]
    seqs[info['chrom_num']][info['start']:info['end']] = seqs[info['chrom_num']][info['start']:info['end']][::-1]
    return seqs

def apply_kataegis(seqs, info):
    seq = seqs[info['chrom_num']][info['start']:info['end']]
    for res in info['res']:
        seqs[info['chrom_num']] = apply_SNP(seq, res)
    return seqs

def apply_translocation(seqs, info):
    seq1 = seqs[info['chrom_num1']]
    seq2 = seqs[info['chrom_num2']]
    if(info['normal']):
        fp1 = seq1[:info['bkpt1']]
        sp1 = seq1[info['bkpt1']:]
        fp2 = seq2[:info['bkpt2']]
        sp2 = seq2[info['bkpt2']:]
        seq1 = fp1 + sp2
        seq2 = fp2 + sp1
    else:
        fp1 = seq1[:info['bkpt1']]
        sp1 = seq1[info['bkpt1']:]
        fp2 = seq2[:info['bkpt2']]
        sp2 = seq2[info['bkpt2']:]
        seq1 = fp1
        seq2 = fp2 + sp2 + sp1
    seqs[info['chrom_num1']] = seq1
    seqs[info['chrom_num2']] = seq2
    return seqs

def apply_chromothripsis(seqs, info):
    fp = seqs[info['chrom_num']][:info['start']]
    sp = seqs[info['chrom_num']][info['end']:]
    mutp = seqs[info['chrom_num']][info['stidx']:info['end']]
    subseq = []
    curridx = 0
    for i in info['bkpts']:
        subseq.append(mutp[curridx:i])
        curridx = i
    subseq.append(mutp[info['bkpts'][-1]:])
    mutp = b''
    for i in info['order']:
        mutp += subseq[i]
    seqs[info['chrom_num']] = fp + mutp + sp
    return seqs

def apply_BFB(seqs, info):
    for bkpt in info['bkpts']:
        fp = seqs[info['chrom_num']][:bkpt]
        np = fp[::-1]
        curr_seq = fp + np
    seqs[info['chrom_num']] = curr_seq
    return seqs

# def apply_chromoplexy(seqs, info):
#     n_chroms = random.randint(2, 6)
#     c = list(range(len(seqs)))
#     chrom = random.sample(c, n_chroms)
#     list_of_start_tels = []
#     list_of_end_tels = []
#     middle_segs = []
#     total_splits = 0
#     for i in range(n_chroms):
#         num_splits = random.randint(2, 20)
#         breakpoints = numpy_choices(len(seqs[chrom[i]]), num_splits) # Returns locations of breakpoints
#         total_splits += num_splits
#         curridx = 0
#         # for i in breakpoints[1:]:
#         #     middle_segs.append(curr_sequence[curridx:i])
#         #     curridx = i
#         # list_of_end_tels.append(curr_sequence[breakpoints[-1]:])
#         # list_of_start_tels.append(curr_sequence[0:breakpoints[0]])
#     usable_splits = total_splits - 2*(n_chroms)
#     splitter = getsummation(usable_splits, n_chroms) # Divides the usable splits in the chromosomes
    
    
#     for i in range(len(info['splits'])):
#         build_str = list_of_start_tels.pop(
#             random.randrange(len(list_of_start_tels)))
#         for segs in range(splitter[i]):
#             build_str = build_str + middle_segs.pop(random.randrange(len(middle_segs)))
#         build_str = build_str + list_of_end_tels.pop(random.randrange(len(list_of_end_tels)))
#         seqs[chrom[i]] = build_str
#     return {'splits': splitter, 'chrom_num': chrom, 'chrom_name': numchrommap[chrom]}