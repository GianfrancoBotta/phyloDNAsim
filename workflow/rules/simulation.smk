configfile: "/cluster/work/bewi/members/gbotta/MosaicSim_speedup/config/config.yaml"

rule simulate_sc_tumor_data_MosaicSim:
    input:
        genome="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/fasta_files/ucsc_hg19.fa",
        binned_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{region}/bed_files/{region}.binned.bed",
        signatures="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/signatures/MosaicSim_signatures.txt"
    output:
        outdir=directory("/cluster/work/bewi/members/gbotta/MosaicSim_speedup/results/simulation/{region}")
    threads: 10
    params:
        yaml_config="/cluster/work/bewi/members/gbotta/MosaicSim_speedup/config/simulate.yaml",
        region=lambda wc: wc.region
    log:
        "logs/simulation/MosaicSim_{region}.log"
    conda:
        "../envs/mosaic_sim.yaml"
    script:
        "../scripts/ms_sim.py"