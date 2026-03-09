from Bio import SeqIO
import os
import json

from apply_mutations import *

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
    
def wgz(floc, ob):
    with open(floc, 'w') as f:
        for i in ob:
            f.write(f'{i}\n')
            
            
info = os.path.join(snakemake.input['info'], 'information_list.json')
full_genome = snakemake.input['genome']
outdir = snakemake.output['outdir']

# Load genome and bed files
chrom_names, chroms = map(list, zip(*((record.id, bytearray(str(record.seq).upper(), 'utf-8')) for record in SeqIO.parse(full_genome, "fasta")))) # Parse sequences and names together
numchrommap = dict(zip(range(len(chrom_names)), chrom_names))
rev_numchrommap = {v: k for k, v in numchrommap.items()}

# Load infos
with open(info, 'r') as f:
    infos = {int(k): v for k, v in json.load(f).items()}
    
# Write genomes
for clone in infos.keys():
    mut_chroms = applyMutations(chroms, infos, clone)
    os.makedirs(outdir, exist_ok=True)
    wgz(os.path.join(outdir, f'{clone}.gz'), mut_chroms)