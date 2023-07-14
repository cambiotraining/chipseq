#!/bin/bash

# load module
module load singularity/current

# load conda environment - has the correct java version for Nextflow
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nextflow

# run from prep_data directory
cd utilities/prep_data/

# launch workflow
nextflow run nf-core/chipseq -profile singularity \
  --input samplesheet.csv \
  --outdir results/nf-chipseq \
  --fasta $(pwd)/resources/GRCh38.109.fasta.gz \
  --gff $(pwd)/resources/GRCh38.109.gff.gz \
  --blacklist $(pwd)/resources/ENCFF356LFX_exclusion_lists.bed.gz \
  --save_reference \
  --macs_gsize 2700000000 \
  --min_reps_consensus 2
