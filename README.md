# phylodnasim

**Phylogenetic DNA tumour simulator** — generates synthetic bulk and
single-cell sequencing reads (FASTQ) from a simulated tumour phylogeny,
including a comprehensive catalogue of somatic mutations:

| Mutation class | Events |
|---|---|
| Small variants | SNV, small insertion, small deletion |
| Copy number | CNV (tandem duplication), aneuploidy |
| Structural | Large deletion, inversion, reciprocal translocation |
| Complex | BFB, chromothripsis, kataegis |

---

## Installation

```bash
pip install phylodnasim
# or, from source:
git clone https://github.com/gbotta7/phylodnasim.git
cd phylodnasim
pip install .
```

---

## Quick start

### 1 — Generate mutations only

```bash
phylodnasim \
  --genome hg38.fa \
  --outdir out/ \
  --num_leaves 5 \
  --generate_mutations
```

Produces `out/information_list.yaml` with all mutation coordinates.

### 2 — Bulk WGS simulation

```bash
phylodnasim \
  --genome hg38.fa \
  --outdir out/ \
  --num_leaves 5 \
  --bulk \
  --coverage 30 \
  --paired \
  --threads 8
```

Outputs `out/bulkleft.fq.gz` (and `out/bulkright.fq.gz` if `--paired`).

### 3 — Single-cell WGS simulation

```bash
phylodnasim \
  --genome hg38.fa \
  --outdir out/ \
  --num_leaves 5 \
  --num_single_cells 50 \
  --cell_coverage 0.05 \
  --r 2 --p 0.5 \
  --threads 8
```

### 4 — Targeted (panel) sequencing

```bash
phylodnasim \
  --genome hg38.fa \
  --bed panel.bed \
  --outdir out/ \
  --num_leaves 5 \
  --targeted \
  --bulk \
  --coverage 500 \
  --paired
```

### 5 — Resume from pre-computed mutations

```bash
phylodnasim \
  --genome hg38.fa \
  --outdir out/ \
  --num_leaves 5 \
  --mut_infos out/information_list.yaml \
  --bulk \
  --coverage 30
```

---

## All options

| Flag | Default | Description |
|---|---|---|
| `--outdir` | *(required)* | Output directory |
| `--genome` | *(required)* | Reference FASTA |
| `--num_leaves` | *(required)* | Number of tumour clones (tree leaves) |
| `--bed` | — | BED panel file (required with `--targeted`) |
| `--mut_infos` | — | Pre-computed mutation YAML (skip tree/mutation step) |
| `--clone_prop` | — | Pre-computed clone proportions file |
| `--prefix` | — | Filename prefix for all outputs |
| `--threads` | 1 | Parallel workers |
| `--generate_mutations` | off | Only generate mutations, then exit |
| `--paired` | off | Paired-end reads |
| `--targeted` | off | Targeted sequencing mode |
| `--bulk` | off | Bulk simulation (default: single-cell) |
| `--num_single_cells` | — | Cells to simulate (single-cell mode) |
| `--coverage` | — | Total coverage |
| `--cell_coverage` | — | Per-cell coverage (single-cell mode) |
| `--read_len` | 150 | Read length (bp) |
| `--frag_len` | 400 | Fragment length (bp) |
| `--error_rate` | 0.001 | Sequencing error rate |
| `--prop_hc` | 0.1 | Proportion of healthy (non-tumour) cells |
| `--alpha` | 10 | Dirichlet concentration for clone proportions |
| `--r` | — | NB dispersion parameter (single-cell mode) |
| `--p` | — | NB success probability (single-cell mode) |
| `--population_size` | 8e8 | Effective population size for coalescent |

---

## Python API

```python
from phylodnasim import (
    getTree, getPaths_and_TimeMatrix, generateOrder,
    saveMutations, applyMutations, getDirichletClone,
    wgsSim_bulk_parallel,
)
```

All public functions are importable directly from `phylodnasim`.

---

## Output files

| File | Description |
|---|---|
| `information_list.yaml` | Mutation catalogue for every node |
| `clones_proportions.txt` | Sampled Dirichlet clone proportions |
| `parameter_list.yaml` | Run parameters |
| `*left.fq.gz` | R1 reads |
| `*right.fq.gz` | R2 reads (paired-end only) |
| `tree_sequence.tree` | msprime tree sequence |
