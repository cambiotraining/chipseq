#!/bin/bash
#SBATCH --account HELARIUTTA-SL2-CPU
#SBATCH --partition icelake-himem
#SBATCH --chdir /rds/project/rds-Wvut5whigiI/chipseq/utilities/prep_data
#SBATCH --output logs/05-diffbind_%j.log
#SBATCH --time 06:00:00
#SBATCH --cpus-per-task 4

set -e # exit upon any error

# load DiffBind environment, created with:
# mamba create --name diffbind -c bioconda bioconductor-diffbind
source $(conda info --base)/etc/profile.d/conda.sh
conda activate diffbind

# create output dir
mkdir -p results/diffbind/

# create temporary script to run analysis
echo '
library(DiffBind)
samplesheet <- read.csv("diffbind_samplesheet.csv")
dba(sampleSheet = samplesheet[samplesheet$Antibody == "H2Bub1", ]) |> 
  dba.count(minOverlap = 2) |> 
  saveRDS("results/diffbind/H2Bub1_dba.count.rds")

dba(sampleSheet = samplesheet[samplesheet$Antibody == "BRD4", ]) |> 
  dba.count(minOverlap = 2) |> 
  saveRDS("results/diffbind/BRD4_dba.count.rds")
' > temp-diffbind.R

# run analysis
Rscript temp-diffbind.R

# remove temporary script
rm temp-diffbind.R