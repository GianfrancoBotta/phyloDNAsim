import os

WORKDIR = "/cluster/work/bewi/members/gbotta/phyloDNAsim"
configfile: os.path.join(WORKDIR, "workflow/config/config.yaml")

# Get signatures for single-cell tumor data simulation
rule download_MosaicSim_signatures:
    output:
        signatures=os.path.join(WORKDIR, "results/prepare/signatures/MosaicSim_signatures.txt")
    params:
        signatures_url=config["signatures_url"]
    log:
        "logs/prepare/signatures/MosaicSim_signatures.log"
    shell:
        r"""
        mkdir -p $(dirname {output.signatures});
        wget -O {output.signatures} {params.signatures_url} &> {log}
        """

rule simulate_sc_tumor_data_phyloDNAsim:
    input:
        genome=os.path.join(WORKDIR, "data/genome.fa"),
        bed=os.path.join(WORKDIR, "data/panel.bed"),
        signatures=os.path.join(WORKDIR, "results/prepare/signatures/MosaicSim_signatures.txt")
    output:
        outdir=directory(os.path.join(WORKDIR, "results/simulation"))
    threads: config["simulate_threads"]
    params:
        yaml_config=os.path.join(WORKDIR, "workflow/config/simulate.yaml"),
    log:
        "logs/simulation/phyloDNAsim.log"
    conda:
        "../envs/phyloDNAsim.yaml"
    script:
        "../../src/phyloDNAsim.py"