#!/bin/bash

# subsampled reads
mkdir -p course_files/participants/data/
cp -r utilities/prep_data/data/reads_sampled course_files/participants/data/reads

# preprocessed data
cp -r utilities/prep_data/results course_files/participants/preprocessed

# space-saving: remove BAM files and replace them with empty mock files
find course_files/participants/preprocessed -type f -name "*.bam" -delete -exec touch {} \;

# space-saving: remove genome directory
rm -r course_files/participants/preprocessed/nf-chipseq/genome

# reference genomes, etc
cp -r utilities/prep_data/resources course_files/participants/

# subset reference to include chromosome 21 only
# awk tip from here: https://onestopdataanalysis.com/get-sequence-fasta/
zcat course_files/participants/resources/GRCh38.109.fasta.gz |\
  awk -v seq="21" -v RS='>' '$1 == seq {print RS $0}' |\
  gzip > course_files/participants/resources/GRCh38.109.chr21.fasta.gz

# space-saving: remove full genome
rm course_files/participants/resources/GRCh38.109.fasta.gz

# nf-core samplesheet
# note the sed in the second command is to include some "FIXME" for exercise (see https://unix.stackexchange.com/a/486896)
head -n 1 utilities/prep_data/samplesheet.csv > course_files/participants/samplesheet.csv
cat utilities/prep_data/samplesheet.csv | grep -e "^brd4_.*_rep[1,2]" | sed '/^brd4/s/mcf7_input_e2/FIXME/g' >> course_files/participants/samplesheet.csv
cat utilities/prep_data/samplesheet.csv | grep -e "^mcf7_input" >> course_files/participants/samplesheet.csv

# copy the original samplesheet to show as example
cp utilities/prep_data/samplesheet.csv course_files/participants/preprocessed/samplesheet.csv

# diffbind samplesheet and RDS files
sed 's/results/preprocessed/' utilities/prep_data/diffbind_samplesheet.csv > course_files/participants/diffbind_samplesheet.csv
cp -r utilities/prep_data/results/diffbind course_files/participants/preprocessed/