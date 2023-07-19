# Packages ----

# load packages
library(rtracklayer) # for importing BED/GFF/etc.
library(plyranges)   # for working with GenomicRanges
library(ChIPseeker)  # to annotate peaks
library(profileplyr) # for profile heatmaps
library(ggplot2)

# change the default ggplot theme
theme_set(theme_classic(base_size = 14))

# read consensus GRanges
brd4_consensus <- readRDS("preprocessed/r_objects/brd4_consensus_granges.rds")


# Profile heatmaps ----

# import the pre-computed matrix from deeptools
brd4_prof <- import_deepToolsMat("preprocessed/deeptools/brd4.log2.chr21.mat.gz")

brd4_prof

# ranges information
rowRanges(brd4_prof)

# sample information
sampleData(brd4_prof)

# subset the object for:
# first 10 peaks
# first 5 bins
# first 3 samples
brd4_prof[1:10, 1:5, 1:3]

# sample 200 peaks randomly
random_peaks <- sample(1:nrow(brd4_prof), 200)

# Plot random peaks to give an idea of occupancy
generateEnrichedHeatmap(brd4_prof[random_peaks, ], 
                        include_group_annotation = FALSE)

## # save as PNG - adjust width, height and resolution to fit your needs
## png("results/brd4_heatmap.png", width = 3000, height = 1500, res = 300)
## generateEnrichedHeatmap(brd4_prof, include_group_annotation = FALSE)
## dev.off()

# read DEGs from Nagarajan 2017
degs <- read.csv("resources/degs_nagarajan2017.csv")

# subset annotated intervals 
brd4_consensus_degs <- brd4_consensus |> 
  filter(geneId %in% degs$ensembl_gene_id)

# subset our profile object to include only the filtered intervals
brd4_prof_degs <- brd4_prof |> 
  subsetByOverlaps(brd4_consensus_degs)
  
# visualise them
generateEnrichedHeatmap(brd4_prof_degs, include_group_annotation = FALSE)


# Exercise ----

# Generate a profile heatmap for the H2Bub1 ChIP

# !!!FIX!!! import the deeptools matrix preprocessed/deeptools/h2bub1.log2.chr21.mat.gz
h2bub1_prof <- FIXME

# Sample 200 random peaks
random_peaks <- sample(1:nrow(h2bub1_prof), 200)

# !!!FIX!!! Visualise them
FIXME