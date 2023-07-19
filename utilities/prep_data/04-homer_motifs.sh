#!/bin/bash
#SBATCH --account HELARIUTTA-SL2-CPU
#SBATCH --partition icelake-himem
#SBATCH --chdir /rds/project/rds-Wvut5whigiI/chipseq/utilities/prep_data
#SBATCH --output logs/04-homer_%a.log
#SBATCH --time 06:00:00
#SBATCH --cpus-per-task 6 
#SBATCH --array 1-2

set -e # exit upon any error

source $(conda info --base)/etc/profile.d/conda.sh
conda activate chipseq

# select one of the antibodies
antibody=$(echo "BRD4,H2Bub1" | cut -d "," -f $SLURM_ARRAY_TASK_ID)

# create output dir
mkdir -p results/homer/$antibody

# run HOMER
# Usage: findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]
findMotifsGenome.pl \
  results/nf-chipseq/bwa/mergedLibrary/macs2/broadPeak/consensus/$antibody/$antibody.consensus_peaks.bed \
  results/nf-chipseq/genome/GRCh38.109.fasta \
  results/homer/$antibody/ \
  -mask -preparse -preparsedDir results/homer/$antibody/preparse \
  -size 200

# remove preparse dir
rm -r results/homer/$antibody/preparse