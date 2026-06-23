import copy

def apply_SNP(seqs, info):
    seqs[info['chrom_num']][info['pos']] = ord(info['char'])
    return seqs

def apply_insertion(seqs, info):
    seqs[info['chrom_num']] = seqs[info['chrom_num']][:info['pos']] + info['insertion'].encode() + seqs[info['chrom_num']][info['pos']:]
    return seqs

def apply_deletion(seqs, info):
    seqs[info['chrom_num']] = seqs[info['chrom_num']][:info['start']] + seqs[info['chrom_num']][info['end']:]
    return seqs

def apply_CNV(seqs, info):
    seqs[info['chrom_num']] = seqs[info['chrom_num']][:info['start']] + seqs[info['chrom_num']][info['start']:info['end']] * info['rep_num'] + seqs[info['chrom_num']][info['end']:]
    return seqs

def apply_aneuploidy(seqs, info):
    seqs[info['chrom_num']] = seqs[info['chrom_num']] * info['rep_num']
    return seqs
    
def apply_inversion(seqs, info):
    seqs[info['chrom_num']][info['start']:info['end']] = seqs[info['chrom_num']][info['start']:info['end']][::-1]
    return seqs

def apply_kataegis(seqs, info):
    for res in info['mut']:
        seqs = apply_SNP(seqs, res)
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
        seq2 = fp2 + sp1 + sp2
    seqs[info['chrom_num1']] = seq1
    seqs[info['chrom_num2']] = seq2
    return seqs

def apply_chromothripsis(seqs, info):
    fp = seqs[info['chrom_num']][:info['start']]
    sp = seqs[info['chrom_num']][info['end']:]
    mutp = seqs[info['chrom_num']][info['start']:info['end']]
    subseqs = []
    curridx = 0
    for i, bkpt in enumerate(info['bkpts']):
        if info['reverse'][i]:
            subseqs.append(mutp[curridx:bkpt])
        else:
            subseqs.append(mutp[curridx:bkpt][::-1])
        curridx = bkpt
    if info['reverse'][-1]:
        subseqs.append(mutp[info['bkpts'][-1]:])
    else:
        subseqs.append(mutp[info['bkpts'][-1]:][::-1])
    mutp = b''
    for i in info['order']:
        mutp += subseqs[i]
    seqs[info['chrom_num']] = fp + mutp + sp
    return seqs

def apply_BFB(seqs, info):
    mod_seq = bytearray()
    for i in range(len(info['bkpts'])):
        if i == 0:
            fp = seqs[info['chrom_num']][:info['bkpts'][i]]
        else:
            fp = seqs[info['chrom_num']][info['bkpts'][i-1]:info['bkpts'][i]]
        mod_seq += fp + fp[::-1]
    mod_seq += seqs[info['chrom_num']][info['bkpts'][-1]:]
    seqs[info['chrom_num']] = mod_seq
    return seqs
    
# def apply_chromoplexy(seqs, info):
#     original_seqs = [copy.deepcopy(seqs[i]) for i in info['chrom_num']]
#     mod_seqs = [bytearray() for _ in range(len(seqs))]
#     for bkpt in info['bkpts']:
#         mod_seqs[info['chrom_num'][bkpt[0]]] += original_seqs[info['chrom_num'][bkpt[0]]][bkpt[1]:bkpt[2]]
#     for chrom in info['chrom_num']:
#         seqs[chrom] = mod_seqs[chrom]
#     return seqs