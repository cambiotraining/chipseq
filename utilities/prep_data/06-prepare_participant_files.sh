#!/bin/bash

########
# Reads

# subsampled reads
mkdir -p course_files/participants/data/
cp -r utilities/prep_data/data/reads_sampled course_files/participants/data/reads


###############
# Preprocessed

# copy all results
cp -r utilities/prep_data/results course_files/participants/preprocessed

# space-saving: remove BAM files and replace them with empty mock files
find course_files/participants/preprocessed -type f -name "*.bam" -delete -exec touch {} \;

# space-saving: remove genome directory
rm -r course_files/participants/preprocessed/nf-chipseq/genome


############
# Resources 

# reference genomes, etc
mkdir -p course_files/participants/resources
#cp -r utilities/prep_data/resources course_files/participants/

# subset reference to include chromosome 21 only
# awk tip from here: https://onestopdataanalysis.com/get-sequence-fasta/
zcat utilities/prep_data/resources/GRCh38.109.fasta.gz |\
  awk -v seq="21" -v RS='>' '$1 == seq {print RS $0}' |\
  gzip > course_files/participants/resources/GRCh38.109.chr21.fasta.gz

# subset annotation to include chromosome 21 only
zcat utilities/prep_data/resources/GRCh38.109.gtf.gz |\
  grep -P "^(#!|21)\b" |\
  gzip > course_files/participants/resources/GRCh38.109.chr21.gtf.gz

# full annotation as well for annotating peaks
cp utilities/prep_data/resources/GRCh38.109.gtf.gz course_files/participants/resources/GRCh38.109.gtf.gz

# subset exclusion lists to include chrom 21 only
zcat utilities/prep_data/resources/ENCFF356LFX_exclusion_lists.bed.gz |\
  grep -P "^21\s+" |\
  gzip > course_files/participants/resources/ENCFF356LFX_exclusion_lists.chr21.bed.gz

# chromosome sizes (useful for plyranges analysis)
cp utilities/prep_data/resources/GRCh38.109.chrom_sizes.tsv course_files/participants/resources/


##############
# Samplesheet

# nf-core samplesheet
# note the sed in the second command is to include some "FIXME" for exercise (see https://unix.stackexchange.com/a/486896)
head -n 1 utilities/prep_data/samplesheet.csv > course_files/participants/samplesheet.csv
cat utilities/prep_data/samplesheet.csv | grep -e "^brd4_.*_rep[1,2]" | sed '/^brd4/s/mcf7_input_e2/FIXME/g' >> course_files/participants/samplesheet.csv
cat utilities/prep_data/samplesheet.csv | grep -e "^mcf7_input" >> course_files/participants/samplesheet.csv

# copy the original samplesheet to show as example
cp utilities/prep_data/samplesheet.csv course_files/participants/preprocessed/samplesheet.csv


###########
# DiffBind

# diffbind samplesheet and RDS files
sed 's/results/preprocessed/g' utilities/prep_data/diffbind_samplesheet.csv > course_files/participants/diffbind_samplesheet.csv
cp -r utilities/prep_data/results/diffbind course_files/participants/preprocessed/


##########
# Scripts

mkdir -p course_files/participants/scripts

# nf-core/chipseq
echo '#!/bin/bash

# run nextflow pipeline - fix the input and outdir first!
nextflow run nf-core/chipseq -profile singularity \
  --max_memory 23.GB --max_cpus 8 \
  --input FIXME \
  --outdir FIXME \
  --fasta $(pwd)/resources/GRCh38.109.chr21.fasta.gz \
  --gff $(pwd)/resources/GRCh38.109.chr21.gff.gz \
  --blacklist $(pwd)/resources/ENCFF356LFX_exclusion_lists.chr21.bed.gz \
  --macs_gsize 40000000 \
  --min_reps_consensus 2' > course_files/participants/scripts/01-chipseq_workflow.sh