#!/bin/bash
#SBATCH --job-name=snakemake_qc
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=40G     
#SBATCH --time=300:00:00      
#SBATCH --output snakemake_qc.log
#SBATCH --mail-type=END
#SBATCH --mail-user=gbotta@ethz.ch

source ~/.bashrc
conda activate snakemake
module load eth_proxy

LOCKDIR=".snakemake/locks"
if [ -d "$LOCKDIR" ]; then
    echo "[INFO] Detected Snakemake lock directory at $LOCKDIR."
    echo "[INFO] Unlocking workflow..."
    snakemake --unlock
else
    echo "[INFO] No lock detected. Proceeding normally."
fi

echo "[INFO] Starting Snakemake run..."
snakemake --use-conda --use-singularity --cores 20 --singularity-args '-B /scratch -B /cluster/work/bewi/members/gbotta:/cluster/work/bewi/members/gbotta:rw' --rerun-incomplete --keep-incomplete