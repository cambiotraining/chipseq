# Packages ----

# load packages
library(rtracklayer) # for importing BED/GFF/etc.
library(plyranges)   # for working with GenomicRanges
library(ChIPseeker)  # to annotate peaks
library(profileplyr) # for profile heatmaps
library(ggplot2)

# change the default ggplot theme
theme_set(theme_classic(base_size = 14))


# Chromosome info ----

# read chromosome sizes (for GRanges annotation)
chroms <- read.table("resources/GRCh38.109.chrom_sizes.tsv", 
                     col.names = c("seqnames", "seqlengths"))

# order chromosomes in a more intuititve manner
# and retain only autosomes (no contigs, no MT)
chroms <- chroms[match(c(1:22, "X", "Y"), chroms$seqnames), ]

# if you had MT, you can use this code to set it as circular
chroms$is_circular <- chroms$seqnames == "MT"

# view table
chroms


# Import peaks ----

# list peak files
brd4_files <- list.files(path = "preprocessed/nf-chipseq", 
                         pattern = "brd4_.*broadPeak", 
                         recursive = TRUE,
                         full.names = TRUE)
names(brd4_files) <- gsub("_peaks.broadPeak", "", 
                          basename(brd4_files))

brd4_files

# take the peak files, and then...
brd4_ranges <- brd4_files |> 
  # ... loop through and import them, and then...
  lapply(import, 
         format = "broadPeak") |> 
  # ... bind them all together
  bind_ranges(.id = "sample")

brd4_ranges

# subset ranges to contain only main chromosomes
brd4_ranges <- brd4_ranges[seqnames(brd4_ranges) %in% chroms$seqnames, ]
seqlevels(brd4_ranges) <- chroms$seqnames
brd4_ranges <- set_genome_info(brd4_ranges, 
                               genome = "GRCh38",
                               seqnames = chroms$seqnames,
                               seqlengths = chroms$seqlengths,
                               is_circular = chroms$is_circular)

# add treatment variable
brd4_ranges <- brd4_ranges |> 
  mutate(treatment = ifelse(grepl("_e2_", sample), "e2", "veh"))


# Coverage ranges ----

# calculate coverage across genome
brd4_coverage <- brd4_ranges |> 
  compute_coverage()

# visualise occupancy rates
brd4_coverage |> 
  # remove intervals with no coverage at all
  filter(score > 0) |> 
  # convert to data.frame
  as.data.frame() |> 
  # barplot of counts for each coverage score
  ggplot(aes(x = score)) + 
  geom_bar() +
  scale_x_continuous(breaks = 1:6) +
  labs(x = "# overlaps")

# get intervals with coverage >= 2
brd4_coverage2 <- brd4_coverage |> 
  filter(score >= 2)

# create consensus peak intervals
brd4_consensus <- brd4_ranges |> 
  # filter to retain ranges with enough coverage
  filter_by_overlaps(brd4_coverage2) |> 
  # merge ranges within 1kb of each other
  reduce(min.gapwidth = 1e3)


# Annotate peaks ----

# import GTF as a TxDb object
genes <- GenomicFeatures::makeTxDbFromGFF("resources/GRCh38.109.gtf.gz")

# we use ChIPseeker to annotate the peaks
brd4_consensus <- brd4_consensus |> 
  annotatePeak(tssRegion = c(-3e3, 3e3),
               TxDb = genes) |> 
  # convert back to GRanges
  as.GRanges()

brd4_consensus

# barplot of annotations
brd4_consensus |> 
  # remove gene IDs from exon/intro annotations for cleaner plot
  mutate(annotation = gsub("Exon .*", "Exon", annotation)) |> 
  mutate(annotation = gsub("Intron .*", "Intron", annotation)) |> 
  # make plot
  as.data.frame() |> 
  ggplot(aes(annotation)) +
  geom_bar() +
  coord_flip()


# Exercise ----

# !!!FIX!!! list peak files 
h2bub1_files <- list.files(path = "preprocessed/nf-chipseq",
                           pattern = "FIXME",
                           recursive = TRUE,
                           full.names = TRUE)
names(h2bub1_files) <- gsub("_peaks.broadPeak", "",
                            basename(h2bub1_files))

# take the peak files, and then...
h2bub1_ranges <- h2bub1_files |>
  # ... loop through and import them, and then...
  lapply(import,
         format = "broadPeak") |>
  # ... bind them all together
  bind_ranges(.id = "sample")

# subset ranges to contain only main chromosomes
h2bub1_ranges <- h2bub1_ranges[seqnames(h2bub1_ranges) %in% chroms$seqnames, ]
seqlevels(h2bub1_ranges) <- chroms$seqnames
h2bub1_ranges <- set_genome_info(h2bub1_ranges,
                                 genome = "GRCh38",
                                 seqnames = chroms$seqnames,
                                 seqlengths = chroms$seqlengths,
                                 is_circular = chroms$is_circular)

# add treatment variable
h2bub1_ranges <- h2bub1_ranges |>
  mutate(treatment = ifelse(grepl("_e2_", sample), "e2", "veh"))

# !!!FIX!!! calculate coverage across genome
h2bub1_coverage <- FIXME

# occupancy rates
h2bub1_coverage |>
  # remove intervals with no coverage at all
  filter(score > 0) |>
  # convert to data.frame
  as.data.frame() |>
  # barplot of counts for each coverage score
  ggplot(aes(x = score)) +
  geom_bar() +
  scale_x_continuous(breaks = 1:6) +
  labs(x = "# overlaps")

# !!!FIX!!! get intervals with coverage >= 2
h2bub1_coverage2 <- FIXME

# !!!FIX!!! create consensus peak intervals
h2bub1_consensus <- h2bub1_ranges |>
  # filter to retain ranges with enough coverage
  FIXME |>
  # merge ranges within 1kb of each other
  reduce(min.gapwidth = 1e3)

# !!!FIX!!! use ChIPseeker to annotate the peaks
h2bub1_consensus <- FIXME

# barplot of annotations
h2bub1_consensus |>
  # remove gene IDs from exon/intro annotations for cleaner plot
  mutate(annotation = gsub("Exon .*", "Exon", annotation)) |>
  mutate(annotation = gsub("Intron .*", "Intron", annotation)) |>
  # make plot
  as.data.frame() |>
  ggplot(aes(annotation)) +
  geom_bar() +
  coord_flip()



# Subset peaks ----

# read DEGs from Nagarajan 2017
degs <- read.csv("resources/degs_nagarajan2017.csv")

# subset annotated intervals
brd4_consensus_degs <- brd4_consensus |> 
  filter(geneId %in% degs$ensembl_gene_id)
