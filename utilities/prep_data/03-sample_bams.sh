#!/bin/bash
#SBATCH -A HELARIUTTA-SL2-CPU
#SBATCH -D /home/hm533/rds/rds-bioinfo-training-Wvut5whigiI/chipseq/utilities/prep_data
#SBATCH -o logs/03-sample_bams.log
#SBATCH -p icelake-himem
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH -a 2,3,5,6,20,21  # only a few samples for the demonstration (saves time)

set -e # exit upon any error

# load conda environment - with samtools and seqtk
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnaseq

# obtain sample information from samplesheet
info=$(cat samplesheet.csv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

# split relevant parts
sample=$(echo $info | cut -d "," -f 1)
fastq1=$(echo $info | cut -d "," -f 2)

# ncpus for samtools (advised to be NCPUs - 1)
ncpus=$(($SLURM_CPUS_PER_TASK - 1))

# chromosomes to sample, separated by space
chroms="21"

# output directory
outdir="data/reads_sampled/"
mkdir -p $outdir

# BAM file path - assuming nf-core/chipseq has been run already on the full data
bam="results/nf-chipseq/bwa/mergedLibrary/${sample}.mLb.clN.sorted.bam"

# get mapped reads
samtools view -F 4 -h -@ $ncpus $bam $chroms |\
  samtools sort -@ $ncpus -u |\
  samtools fastq -@ $ncpus |\
  sed -n '1~4p' | sed 's/^@//' > ${outdir}/${sample}_aligned_read_names.txt

# sample reads
seqtk subseq $fastq1 ${outdir}/${sample}_aligned_read_names.txt |\
  # seqtk sample -s1 - 200000 |\
  gzip > ${outdir}/$(basename $fastq1)

# clean things
rm ${outdir}/${sample}_aligned_read_names.txt
echo "Finished ${sample} (array # ${SLURM_ARRAY_TASK_ID})"