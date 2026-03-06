import os

WORKDIR = "/cluster/work/bewi/members/gbotta/phyloDNAsim"
configfile: os.path.join(WORKDIR, "workflow/config/config.yaml")

rule simulate_sc_tumor_data_phyloDNAsim:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/simulation/clone_{n_clones}"))
    threads: config["simulate_threads"]
    params:
        num_samples = 1,
        num_sc_list = " ".join(map(str, [500, 1000, 2000, 4000, 8000])),
        read_len_tb = 150,
        read_len_lb = 140,
        frag_len = 400,
        prop_hc_tb = 0.1,
        prop_hc_pb = 0.9,
        cell_coverage_tb = 250,
        coverage_lb = 2000
    log:
        "logs/simulation/clone_{n_clones}/phyloDNAsim.log"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --threads {threads} \
            --paired \
            --targeted \
            --num_samples {params.num_samples} \
            --num_leaves {wildcards.n_clones} \
            --num_single_cells_list {params.num_sc_list} \
            --read_len_tb {params.read_len_tb} \
            --frag_len_tb {params.frag_len} \
            --prop_hc_tb {params.prop_hc_tb} \
            --alpha 10 \
            --r 3.78 \
            --p 0.0143 \
            --cell_coverage_tb {params.cell_coverage_tb} \
            --tumor_liquid_biopsy \
            --read_len_lb {params.read_len_lb} \
            --frag_len_lb {params.frag_len} \
            --coverage_lb {params.coverage_lb} \
            --prop_hc_tlb {params.prop_hc_tb} \
            --plasma_liquid_biopsy \
            --prop_hc_plb {params.prop_hc_pb} \
            &> {log}
        """