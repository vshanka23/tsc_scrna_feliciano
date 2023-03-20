#!/bin/bash
#
#SBATCH --job-name=TSCscRNA
#SBATCH --ntasks=1   
#SBATCH --partition=bigmem
#SBATCH --time=30-00:00:00
#SBATCH --mem=2gb
#SBATCH --output=/data/gsl/bsl/vshanka/vshanka_tsc_scrna/Counts/log/test_output_%j.txt
#SBATCH --error=/data/gsl/bsl/vshanka/vshanka_tsc_scrna/Counts/log/test_error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

cd /data/gsl/bsl/vshanka/vshanka_tsc_scrna/Counts
#mkdir -p ./{log,logs_slurm}

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

#--dag | display | dot
#-p -n \

snakemake \
-s Snakefile \
--profile slurm \
--configfile cellranger_counts.yaml \
--latency-wait 120
