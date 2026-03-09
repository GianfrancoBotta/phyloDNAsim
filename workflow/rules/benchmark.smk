import os

WORKDIR = "/cluster/work/bewi/members/gbotta/phyloDNAsim"
configfile: os.path.join(WORKDIR, "workflow/config/config.yaml")


### GENERATE MUTATIONS
rule generate_mutations_targeted:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/targeted/clone_{n_clones}"))
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

rule generate_mutations_wgs:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/wgs/clone_{n_clones}"))
    log:
        "logs/benchmark/phyloDNAsim/wgs/mutations/clone_{n_clones}/generate_muts.log"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --generate_mutations \
            --num_leaves {wildcards.n_clones} \
            &> {log}
        """



### MUTATE GENOMES FOR MOSAICSIM
rule generate_MosaicSim_genomes_bulk:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        info=os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/{sim_type}/clone_{n_clones}")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/bulk/clone_{n_clones}/tumor_0"))
    priority: 1
    log:
        "logs/benchmark/MosaicSim/{sim_type}/bulk/genomes/clone_{n_clones}/generate_genomes.log"
    benchmark:
        "benchmarks/MosaicSim/{sim_type}/bulk/clone_{n_clones}/generate_genomes.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    script:
        "../scripts/py/generate_MosaicSim_genome.py"

rule generate_MosaicSim_genomes_sc:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        info=os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/{sim_type}/clone_{n_clones}")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/sc/clone_{n_clones}/tumor_0"))
    priority: 1
    log:
        "logs/benchmark/MosaicSim/{sim_type}/sc/genomes/clone_{n_clones}/generate_genomes.log"
    benchmark:
        "benchmarks/MosaicSim/{sim_type}/sc/clone_{n_clones}/generate_genomes.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    script:
        "../scripts/py/generate_MosaicSim_genome.py"



### MOSAICSIM SIMULATE
rule benchmark_MosaicSim_bulk:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        signatures=os.path.join(CLONEDIR, "data/signatures.txt"),
        genomedir=os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/bulk/clone_{n_clones}/tumor_0")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/bulk/clone_{n_clones}/reference"))
    priority: 1
    params:
        targeted = lambda wc: True if wc.sim_type == "targeted" else False,
        n_cells = 0,
        n_clones = lambda wc: int(wc.n_clones),
        coverage = lambda wc: 10000 if wc.sim_type == "targeted" else 10
    log:
        "logs/benchmark/MosaicSim/{sim_type}/bulk/clone_{n_clones}/MosaicSim.log"
    benchmark:
        "benchmarks/MosaicSim/{sim_type}/bulk/clone_{n_clones}/MosaicSim.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    script:
        "../MosaicSim_revised.py"

rule benchmark_MosaicSim_sc:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        signatures=os.path.join(CLONEDIR, "data/signatures.txt"),
        genomedir=os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/sc/clone_{n_clones}/tumor_0")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/sc/clone_{n_clones}/reference"))
    priority: 1
    params:
        tool_dir=CLONEDIR,
        targeted = lambda wc: True if wc.sim_type == "targeted" else False,
        n_cells = 100,
        n_clones = lambda wc: int(wc.n_clones),
        coverage = lambda wc: 100 if wc.sim_type == "targeted" else 0.1
    log:
        "logs/benchmark/MosaicSim/{sim_type}/sc/clone_{n_clones}/MosaicSim.log"
    benchmark:
        "benchmarks/MosaicSim/{sim_type}/sc/clone_{n_clones}/MosaicSim.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    script:
        "../MosaicSim_revised.py"



### PHYLODNASIM SIMULATE
rule benchmark_phyloDNAsim_sc_targeted:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        infos_folder=os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/targeted/clone_{n_clones}")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/core_{n_cores}/targeted/sc/clone_{n_clones}"))
    params:
        n_cells = config["benchmark_cells"],
        read_len = 150,
        frag_len = 400,
        prop_hc = 0.1,
        cell_coverage = 100
    threads:
        lambda wc: int(wc.n_cores)
    log:
        "logs/benchmark/phyloDNAsim/core_{n_cores}/targeted/sc/clone_{n_clones}/phyloDNAsim.log"
    benchmark:
        "benchmarks/phyloDNAsim/core_{n_cores}/targeted/sc/clone_{n_clones}/phyloDNAsim.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --threads {threads} \
            --paired \
            --targeted \
            --num_leaves {wildcards.n_clones} \
            --num_single_cells {params.n_cells} \
            --read_len {params.read_len} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc} \
            --alpha 10 \
            --r 3.78 \
            --p 0.0143 \
            --cell_coverage {params.cell_coverage_tb} \
            &> {log}
        """

rule benchmark_phyloDNAsim_bulk_targeted:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        infos_folder=os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/targeted/clone_{n_clones}")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/core_{n_cores}/targeted/bulk/clone_{n_clones}"))
    params:
        read_len = 150,
        frag_len = 400,
        prop_hc = 0.1,
        coverage = 10000
    threads:
        lambda wc: int(wc.n_cores)
    log:
        "logs/benchmark/phyloDNAsim/core_{n_cores}/targeted/bulk/clone_{n_clones}/phyloDNAsim.log"
    benchmark:
        "benchmarks/phyloDNAsim/core_{n_cores}/targeted/bulk/clone_{n_clones}/phyloDNAsim.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --threads {threads} \
            --bulk \
            --paired \
            --targeted \
            --num_leaves {wildcards.n_clones} \
            --read_len {params.read_len} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc} \
            --alpha 10 \
            --coverage {params.coverage} \
            &> {log}
        """

rule benchmark_phyloDNAsim_sc_wgs:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        infos_folder=os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/wgs/clone_{n_clones}")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/core_{n_cores}/wgs/sc/clone_{n_clones}"))
    params:
        n_cells = config["benchmark_cells"],
        read_len = 150,
        frag_len = 400,
        prop_hc = 0.1,
        cell_coverage = 0.1
    threads:
        lambda wc: int(wc.n_cores)
    log:
        "logs/benchmark/phyloDNAsim/core_{n_cores}/wgs/sc/clone_{n_clones}/phyloDNAsim.log"
    benchmark:
        "benchmarks/phyloDNAsim/core_{n_cores}/wgs/sc/clone_{n_clones}/phyloDNAsim.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --threads {threads} \
            --paired \
            --num_leaves {wildcards.n_clones} \
            --num_single_cells {params.n_cells} \
            --read_len {params.read_len} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc} \
            --alpha 10 \
            --r 3.78 \
            --p 0.0143 \
            --cell_coverage {params.cell_coverage} \
            &> {log}
        """

rule benchmark_phyloDNAsim_bulk_wgs:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        infos_folder=os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/mutations/wgs/clone_{n_clones}")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/core_{n_cores}/wgs/bulk/clone_{n_clones}"))
    params:
        read_len = 150,
        frag_len = 400,
        prop_hc = 0.1,
        coverage = 10
    threads:
        lambda wc: int(wc.n_cores)
    log:
        "logs/benchmark/phyloDNAsim/core_{n_cores}/wgs/bulk/clone_{n_clones}/phyloDNAsim.log"
    benchmark:
        "benchmarks/phyloDNAsim/core_{n_cores}/wgs/bulk/clone_{n_clones}/phyloDNAsim.txt"
    conda:
        "../envs/phyloDNAsim.yaml"
    shell:
        """
        python ../src/phyloDNAsim.py \
            --outdir {output.outdir} \
            --genome {input.genome} \
            --bed {input.bed} \
            --mut_infos {input.infos_folder}/information_list.json \
            --threads {threads} \
            --bulk \
            --paired \
            --num_leaves {wildcards.n_clones} \
            --read_len {params.read_len} \
            --frag_len {params.frag_len} \
            --prop_hc {params.prop_hc} \
            --alpha 10 \
            --coverage {params.coverage} \
            &> {log}
        """



### GET STATS
rule get_sim_stats:
    input:
        MosaicSim_benchs=os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/{experiment}/clone_{n_clones}"),
        MosaicSim_reads=os.path.join(WORKDIR, "results/benchmark/MosaicSim/{sim_type}/{experiment}/clone_{n_clones}"),
        phyloDNAsim=expand(
            os.path.join(WORKDIR, "results/benchmark/phyloDNAsim/core_{n_cores}/{sim_type}/{experiment}/clone_{n_clones}"),
            n_cores=config["benchmark_core_list"],
            sim_type=["{sim_type}"],
            experiment=["{experiment}"])
    output:
        outdir=directory(os.path.join(WORKDIR, "results/benchmark/stats"))
    params:
        num_clones_list=config["benchmark_clone_list"],
        num_cores_list=config["benchmark_core_list"],
        num_cells=config[""]
    log:
        "logs/benchmark/stats/get_stats.log"
    conda:
        "../envs/phyloDNAsim.yaml"
    script:
        "../get_stats.py"