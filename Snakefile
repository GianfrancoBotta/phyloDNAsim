from snakemake.utils import min_version
import os

configfile: "/cluster/work/bewi/members/gbotta/MosaicSim_speedup/config/config.yaml"

min_version(config["snakemake_min_version"])

container: f"docker://condaforge/mambaforge:{config['mambaforge_version']}" # to allow reproducibility

include: "workflow/rules/simulation.smk"

rule all:
    input:
        expand(
            "/cluster/work/bewi/members/gbotta/MosaicSim_speedup/results/simulation/{region}",
            region=config["chr_regions"]
        )
