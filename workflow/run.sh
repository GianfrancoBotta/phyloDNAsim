#!/bin/bash
#SBATCH --job-name=run_simulation
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=16G     
#SBATCH --time=300:00:00      
#SBATCH --output simulation.log
#SBATCH --mail-type=END
#SBATCH --mail-user=gbotta@ethz.ch

source /cluster/home/gbotta/miniforge3/etc/profile.d/conda.sh
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
snakemake --use-conda --use-singularity --cores 50 --singularity-args '-B /scratch -B /cluster/work/bewi/members/gbotta:/cluster/work/bewi/members/gbotta:rw' --rerun-incomplete --keep-incomplete