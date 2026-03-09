import os

WORKDIR = "/cluster/work/bewi/members/gbotta/phyloDNAsim"
configfile: os.path.join(WORKDIR, "workflow/config/config.yaml")

### GENERATE MUTATIONS
rule generate_mutations_simulation:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/simulation/mutations/clone_{n_clones}"))
    log:
        "logs/benchmark/targeted/phyloDNAsim/mutations/clone_{n_clones}/generate_muts.log"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --generate_mutations \
            --num_leaves {wildcards.n_clones} \
            &> {log}
        """


rule simulate_tumor_data_sc_dna_seq_phyloDNAsim:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        infos_folder=os.path.join(WORKDIR, "results/simulation/mutations/clone_{n_clones}")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/simulation/clone_{n_clones}/cell_{n_cells}"))
    threads: config["simulate_threads"]
    params:
        num_samples = 1,
        read_len = 150,
        read_len_lb = 140,
        frag_len = 400,
        prop_hc = 0.1,
        prop_hc_plb = 0.9,
        prop_hc_ref = 1,
        cell_coverage = 200,
        coverage = 2000,
        prefix_tlb = "T",
        prefix_plb = "P",
        prefix_ref = "N"
    log:
        "logs/simulation/clone_{n_clones}/cell_{n_cells}/phyloDNAsim_sc.log"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        ### SINGLE-CELL
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --threads {threads} \
            --paired \
            --targeted \
            --num_samples {params.num_samples} \
            --num_leaves {wildcards.n_clones} \
            --num_single_cells {wildcards.n_cells} \
            --read_len {params.read_len} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc} \
            --alpha 10 \
            --r 3.78 \
            --p 0.0143 \
            --cell_coverage {params.cell_coverage} \
            &> {log};
        ### TUMOR LIQUID BIOPSY
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --clone_prop {output.outdir}/clones_proportions.txt \
            --threads {threads} \
            --paired \
            --bulk \
            --targeted \
            --num_samples {params.num_samples} \
            --num_leaves {wildcards.n_clones} \
            --read_len {params.read_len_lb} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc} \
            --coverage {params.coverage} \
            --prefix {params.prefix_tlb}
            &> {log};
        ### PLASMA LIQUID BIOPSY
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --clone_prop {output.outdir}/clones_proportions.txt \
            --threads {threads} \
            --paired \
            --bulk \
            --targeted \
            --num_samples {params.num_samples} \
            --num_leaves {wildcards.n_clones} \
            --read_len {params.read_len_lb} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc_plb} \
            --coverage {params.coverage} \
            --prefix {params.prefix_plb}
            &> {log};
        ### REFERENCE
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --clone_prop {output.outdir}/clones_proportions.txt \
            --threads {threads} \
            --paired \
            --bulk \
            --targeted \
            --num_samples {params.num_samples} \
            --num_leaves {wildcards.n_clones} \
            --read_len {params.read_len_lb} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc_ref} \
            --coverage {params.coverage} \
            --prefix {params.prefix_ref}
            &> {log};
        """