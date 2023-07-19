#!/bin/bash

cd utilities/prep_data/
mkdir resources

# download latest genome version from ENSEMBL
# see https://nf-co.re/usage/reference_genomes#custom-genomes
version=109
wget -L -O resources/GRCh38.${version}.fasta.gz ftp://ftp.ensembl.org/pub/release-${version}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L -O resources/GRCh38.${version}.gtf.gz ftp://ftp.ensembl.org/pub/release-${version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${version}.gtf.gz
wget -L -O resources/GRCh38.${version}.gff.gz ftp://ftp.ensembl.org/pub/release-${version}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${version}.gff3.gz

# get chromosome lengths
# this approach misses contigs
# zcat resources/GRCh38.${version}.gff.gz | awk -v OFS='\t' '{ if ($3 == "chromosome") { print $1, $5 - $4 + 1}}' > resources/GRCh38.${version}.chrom_sizes.tsv

# get chromosome lengths from FASTA
# https://www.danielecook.com/generate-fasta-sequence-lengths/
# note the sed is to remove extra information in the sequence names
zcat resources/GRCh38.${version}.fasta.gz | sed '/>/ { s/ .*// }' | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > resources/GRCh38.${version}.chrom_sizes.tsv

# download human exclusion list regions from Encode
# this is equivalent to v3 of the nf-core/chipseq pipeline: https://github.com/nf-core/chipseq/blob/master/assets/blacklists/
wget -O temp.bed.gz https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
# remove "chr" from chromosome names to match ENSEMBL naming
zcat temp.bed.gz | sed 's/chr//' | gzip > resources/ENCFF356LFX_exclusion_lists.bed.gz
rm temp.bed.gz
