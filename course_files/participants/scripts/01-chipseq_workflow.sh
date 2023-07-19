#!/bin/bash

# run nextflow pipeline - fix the input and outdir first!
nextflow run nf-core/chipseq -profile singularity \
  --max_memory 23.GB --max_cpus 8 \
  --input FIXME \
  --outdir FIXME \
  --fasta $(pwd)/resources/GRCh38.109.chr21.fasta.gz \
  --gff $(pwd)/resources/GRCh38.109.chr21.gff.gz \
  --blacklist $(pwd)/resources/ENCFF356LFX_exclusion_lists.chr21.bed.gz \
  --macs_gsize 40000000 \
  --min_reps_consensus 2
