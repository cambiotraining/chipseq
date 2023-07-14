#!/bin/bash
#SBATCH -A HELARIUTTA-SL2-CPU
#SBATCH -D /home/hm533/rds/rds-bioinfo-training-Wvut5whigiI/chipseq
#SBATCH -o logs/03-sample_fastq.log
#SBATCH -p icelake-himem
#SBATCH -c 4
#SBATCH -t 02:00:00

# exit when any command fails
set -e 

# All filepaths relative to parent directory

# load conda environment - with seqtk
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ngs

# moved raw reads to prep_data folder ahead of time
# mkdir utilities/prep_data/reads
# mv course_files/participants/data/reads/* utilities/prep_data/reads/

# loop through fastq files
for fastq in $(ls utilities/prep_data/reads/)
do 
  # sample reads
  seqtk sample -s1 utilities/prep_data/reads/$fastq 1000000 |\
    gzip > course_files/participants/data/reads/$fastq
  
  echo "Finished $fastq"
done
