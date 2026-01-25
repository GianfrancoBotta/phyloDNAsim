from pydantic import BaseModel, Field, conint
from typing import List
import copy

class ApplySNPInfo(BaseModel):
    chrom_num: conint(ge=0)
    pos: conint(ge=0)
    char: str = Field(min_length=1, max_length=1)

def apply_SNP(seqs, info: ApplySNPInfo):
    seqs[info.chrom_num][info.pos] = ord(info.char)
    return seqs

class ApplyInsertionInfo(BaseModel):
    chrom_num: conint(ge=0)
    pos: conint(ge=0)
    insertion: str = Field(min_length=1)

def apply_insertion(seqs, info: ApplyInsertionInfo):
    seqs[info.chrom_num] = seqs[info.chrom_num][:info.pos] + info.insertion.encode() + seqs[info.chrom_num][info.pos:]
    return seqs

class ApplyDeletionInfo(BaseModel):
    chrom_num: conint(ge=0)
    start: conint(ge=0)
    end: conint(gt=0)

def apply_deletion(seqs, info: ApplyDeletionInfo):
    seqs[info.chrom_num] = seqs[info.chrom_num][:info.start] + seqs[info.chrom_num][info.end:]
    return seqs

class ApplyCNVInfo(BaseModel):
    chrom_num: conint(ge=0)
    start: conint(ge=0)
    end: conint(gt=0)
    rep_num: conint(ge=2)

def apply_CNV(seqs, info: ApplyCNVInfo):
    seqs[info.chrom_num] = seqs[info.chrom_num][:info.start] + seqs[info.chrom_num][info.start:info.end]*info.rep_num + seqs[info.chrom_num][info.end:]
    return seqs

class ApplyAneuploidyInfo(BaseModel):
    chrom_num: conint(ge=0)
    rep_num: conint(ge=0)
    
def apply_aneuploidy(seqs, info: ApplyAneuploidyInfo):
    seqs[info.chrom_num] = seqs[info.chrom_num] * info.rep_num
    return seqs

class ApplyInversionInfo(BaseModel):
    chrom_num: conint(ge=0)
    start: conint(ge=0)
    end: conint(gt=0)
    
def apply_inversion(seqs, info: ApplyInversionInfo):
    seqs[info.chrom_num][info.start:info.end] = seqs[info.chrom_num][info.start:info.end][::-1]
    return seqs

class ApplyKataegisInfo(BaseModel):
    chrom_num: conint(ge=0)
    start: conint(ge=0)
    end: conint(gt=0)
    res: List[ApplySNPInfo]

def apply_kataegis(seqs, info: ApplyKataegisInfo):
    seq = seqs[info.chrom_num][info.start:info.end]
    for res in info['res']:
        res_val = ApplySNPInfo(res)
        seqs[info.chrom_num] = apply_SNP(seq, res_val)
    return seqs

class ApplyTranslocationInfo(BaseModel):
    chrom_num1: conint(ge=0)
    chrom_num2: conint(ge=0)
    bkpt1: conint(ge=0)
    bkpt2: conint(ge=0)
    normal: bool

def apply_translocation(seqs, info: ApplyTranslocationInfo):
    seq1 = seqs[info.chrom_num1]
    seq2 = seqs[info.chrom_num2]
    if(info.normal):
        fp1 = seq1[:info.bkpt1]
        sp1 = seq1[info.bkpt1:]
        fp2 = seq2[:info.bkpt2]
        sp2 = seq2[info.bkpt2:]
        seq1 = fp1 + sp2
        seq2 = fp2 + sp1
    else:
        fp1 = seq1[:info.bkpt1]
        sp1 = seq1[info.bkpt1:]
        fp2 = seq2[:info.bkpt2]
        sp2 = seq2[info.bkpt2:]
        seq1 = fp1
        seq2 = fp2 + sp2 + sp1
    seqs[info.chrom_num1] = seq1
    seqs[info.chrom_num2] = seq2
    return seqs

class ApplyChromothripsisInfo(BaseModel):
    chrom_num: conint(ge=0)
    start: conint(ge=0)
    end: conint(gt=0)
    stidx: conint(ge=0)
    bkpts: List[conint(gt=0)]
    order: List[int]

def apply_chromothripsis(seqs, info: ApplyChromothripsisInfo):
    fp = seqs[info.chrom_num][:info.start]
    sp = seqs[info.chrom_num][info.end:]
    mutp = seqs[info.chrom_num][info.start:info.end]
    subseq = []
    curridx = 0
    for i in range(info.bkpts):
        bkpt = info.bkpts[i]
        if info.reverse[i]:
            subseq.append(mutp[curridx:bkpt])
        else:
            subseq.append(mutp[curridx:bkpt][::-1])
        curridx = bkpt
    if info.reverse[-1]:
        subseq.append(mutp[info.bkpts[-1]:])
    else:
        subseq.append(mutp[info.bkpts[-1]:][::-1])
    mutp = b''
    for i in info.order:
        mutp += subseq[i]
    seqs[info.chrom_num] = fp + mutp + sp
    return seqs

class ApplyBFBInfo(BaseModel):
    chrom_num: conint(ge=0)
    bkpts: List[conint(gt=0)]

def apply_BFB(seqs, info: ApplyBFBInfo):
    for bkpt in info.bkpts:
        fp = seqs[info.chrom_num][:bkpt]
        np = fp[::-1]
        curr_seq = fp + np
    seqs[info.chrom_num] = curr_seq
    return seqs

class ApplyChromoplexyInfo(BaseModel):
    chrom_num: conint(ge=0)
    bkpts: List[conint(ge=0)]
    
def apply_chromoplexy(seqs, info: ApplyChromoplexyInfo):
    original_seqs = [copy.deepcopy(seqs[i]) for i in info.chrom_num]
    mod_seqs = [bytearray() for _ in range(len(seqs))]
    for bkpt in info.bkpts:
        mod_seqs[info.chrom_num[bkpt[0]]] += original_seqs[info.chrom_num[bkpt[0]]][bkpt[1]:bkpt[2]]
    for chrom in info.chrom_num:
        seqs[chrom] = mod_seqs[chrom]
    return seqs