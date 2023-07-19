# Packages ----

# load packages
library(rtracklayer) # for importing BED/GFF/etc.
library(plyranges)   # for working with GenomicRanges
library(ChIPseeker)  # to annotate peaks
library(profileplyr) # for profile heatmaps
library(DiffBind)    # for ChIP peak analysis
library(ggplot2)

# change the default ggplot theme
theme_set(theme_classic(base_size = 14))

# read DiffBind samplesheet
samplesheet <- read.csv("diffbind_samplesheet.csv")

head(samplesheet)

# create DBA object for BRD4 antibody
brd4_dba <- dba(sampleSheet = samplesheet[samplesheet$Antibody == "BRD4", ])

brd4_dba

# correlation plot using caller score
plot(brd4_dba)

# count reads overlapping peaks 
# this takes a long time to run! So we load pre-computed one
# brd4_dba <- dba.count(brd4_dba)
brd4_dba <- readRDS("preprocessed/r_objects/brd4_dba.count.rds")

# correlation plot based on raw counts
plot(brd4_dba)

# normalise counts by library size (default)
brd4_dba <- dba.normalize(brd4_dba, normalize = DBA_NORM_LIB)

# set contrast 
brd4_dba <- dba.contrast(brd4_dba, 
                         design = ~ Treatment,
                         reorderMeta = list(Treatment = "veh"),
                         minMembers = 3)

# run the analysis
brd4_dba <- dba.analyze(brd4_dba, 
                        method = DBA_DESEQ2,
                        bBlacklist = FALSE, 
                        bGreylist = FALSE)

# extract diffbound sites
# keep all peaks, even those that are non-significant
brd4_diffbound <- dba.report(brd4_dba, th = 1)

brd4_diffbound

# count how many up or down
brd4_diffbound |> 
  filter(FDR < 0.05) |> 
  summarise(up = sum(Fold > 0), down = sum(Fold < 0))


# Visualisation ----

# PCA plot
dba.plotPCA(brd4_dba, label = DBA_REPLICATE)

# MA plot
dba.plotMA(brd4_dba)

# can also do it with ggplot2
brd4_diffbound |>
  as.data.frame() |>
  mutate(sig = ifelse(FDR < 0.05, Fold, NA)) |>
  ggplot(aes(Conc, Fold)) +
  geom_point(colour = "grey") +
  geom_point(aes(y = sig), colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed")

# volcano plot
dba.plotVolcano(brd4_dba)

# can also do it with ggplot2
brd4_diffbound |>
  as.data.frame() |>
  mutate(sig = ifelse(FDR < 0.05, Fold, NA)) |>
  ggplot(aes(Fold, -log10(FDR))) +
  geom_point(colour = "grey") +
  geom_point(aes(x = sig), colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed")

# boxplot of normalised counts
dba.plotBox(brd4_dba)

# heatmap of DB peaks
dba.plotHeatmap(brd4_dba, contrast = 1, 
                correlations = FALSE, scale = "row")


# Pipes ----

# # full pipeline - do not run, it will take too long!
# brd4_dba <- dba(sampleSheet = samplesheet[samplesheet$Antibody == "BRD4", ]) |>
#   dba.count() |>
#   dba.normalize(normalize = DBA_NORM_LIB) |>
#   dba.contrast(reorderMeta = list(Treatment = "veh")) |>
#   dba.analyze(method = DBA_DESEQ2,
#               bBlacklist = FALSE,
#               bGreylist = FALSE)
# 
# # extract diffbound sites
# brd4_diffbound <- dba.report(brd4_dba, th = 1)


# Exercise ----

# !!!FIX!!! run diffbind analysis
h2bub1_dba <- readRDS("preprocessed/r_objects/h2bub1_dba.count.rds") |>
  FIXME

# !!!FIX!!! extract results as a GRanges object
brd4_diffbound <- FIXME

# !!!FIX!!! correlation heatmap for the samples
FIXME

# !!!FIX!!! PCA plot
FIXME

# !!!FIX!!! MA plot
FIXME

# !!!FIX!!! boxplot
FIXME

# !!!FIX!!! heatmap of differentially bound peaks
FIXME
