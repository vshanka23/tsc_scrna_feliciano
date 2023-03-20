#!/bin/bash
#
#SBATCH --job-name=TSCscRNA
#SBATCH --ntasks=1   
#SBATCH --partition=compute
#SBATCH --time=30-00:00:00
#SBATCH --mem=2gb
#SBATCH --output=/data/gsl/bsl/vshanka/vshanka_tsc_scrna/output_%j.txt
#SBATCH --error=/data/gsl/bsl/vshanka/vshanka_tsc_scrna/error_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=vshanka@clemson.edu

cd /data/gsl/bsl/vshanka/vshanka_tsc_scrna/FASTQS

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate gcc9_libstdc6
ml cellranger/7.0.1

cellranger mkfastq --id=TSCscRNA \
--run=/data/gsl/bsl/vshanka/vshanka_tsc_scrna/BCL/Files \
--samplesheet=/data/gsl/bsl/vshanka/vshanka_tsc_scrna/BCL/Files/E6440_NovaSeq_SampleSheet_DF_sc_09_26_2022.csv \
--jobmode=slurm
