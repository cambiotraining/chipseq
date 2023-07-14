#!/bin/bash
#SBATCH -A HELARIUTTA-SL2-CPU
#SBATCH -D /home/hm533/rds/rds-bioinfo-training-Wvut5whigiI/chipseq/utilities/prep_data
#SBATCH -o logs/00-fastqdump_%a.log
#SBATCH -p icelake
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH -a 19-20

source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra

# fetch current SRA id
sra=$(cat nagarajan2014_sra_accessions.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# prefetch
prefetch ${sra}

# validate
vdb-validate ${sra}

# convert
fasterq-dump --outdir  reads/ ${sra}

# gzip
cd reads/
gzip ${sra}.fastq
