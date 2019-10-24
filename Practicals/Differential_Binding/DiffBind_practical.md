# Differential binding analysis practical
(based on the practical by Bori Mifsud)

## Resources used

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[R](https://cran.r-project.org/)

[DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)


## Introduction

One of the common downstream analyses of ChIP-seq data is comparing the binding
profile of a transcription factor or histone modification in different
conditions (e.g. healty vs. disease, wild-type vs. mutant, treatment vs.
non-treatment, two different cell types ...). In this practical we will do
differential binding analysis on a H3K27ac (active enhancer mark) dataset
comparing two mutant Drosophila cell lines (gd7 and tl10b), which correspond to
mesodermal and ectodermal precursor cells (GEO: GSE68983) from Koenecke et al.,
Genome Biology (2016).

## Differential binding

Finding differentially bound regions in the genome is analogous to identifying
differentially expressed genes in RNA-seq data. In both cases we are dealing
with count data summarised over features (genes/transcripts in the case of
RNA-seq and peaks in the case of ChIP-seq). In both cases the biological
replicates show larger variability than technical replicates, and the negative
binomial model is suitable to compare binding affinities accross samples. 

There are a number of differential expression packages in R that use the
negative binomial model e.g. DESeq2 and edgeR. These methods are wrapped in the
DiffBind package that is geared towards analysing differential binding in
ChIP-seq data and provides a number of analytical plots as well. In this
practical we will use the DiffBind R package to identify differentially
acetylated regions in the two Drosophila cell lines.

## Quality control and alignment

Reads were downloaded and converted into fastq format, read quality was also
investigated with fastqc and quality looked ok. Reads were mapped onto the
Drosophila genome (dm3). In order to keep the dataset small, we used chr2L only
from the reference genome. For mapping we used the aligner bowtie and we used
samtools to process alignments and sort the bam files according to genome
coordinates. PCR duplicates can be filtered from the sorted bam files using the
samtools markdup function.

All the steps above were already performed for you so you do not need to run
them again on your computer.

## Peak calling

Files were renamed to reflect cell type and sample information: 

* **gd7** - gd7_K27ac_rep1, gd7_K27ac_rep2, gd7_input_rep1, gd7_input_rep2 
* **tl10b** - tl10b_K27ac_rep1, tl10b_K27ac_rep2, tl10b_input_rep1,
  tl10b_input_rep2

**Note**: The tl10b_input_rep1 had very low number of reads therefore we use the second
replicate as input for both tl10b_K27ac_rep samples. 

The bam files are located in ~/Desktop/Course_Materials/Practicals/DiffBind/bam.

|CellLine |SampleType |Replicate |Deduplicated |Filename                         |
|:--------|:----------|:---------|:------------|:--------------------------------|
|gd7      |input1     |          |No           |gd7_input1_chr2L.bam             |
|gd7      |input1     |          |Yes          |gd7_input1_chr2L.dedup.bam       |
|gd7      |input2     |          |No           |gd7_input2_chr2L.bam             |
|gd7      |input2     |          |Yes          |gd7_input2_chr2L.dedup.bam       |
|gd7      |H3K27ac    |rep1      |No           |gd7_K27ac_rep1_chr2L.bam         |
|gd7      |H3K27ac    |rep1      |Yes          |gd7_K27ac_rep1_chr2L.dedup.bam   |
|gd7      |H3K27ac    |rep2      |No           |gd7_K27ac_rep2_chr2L.bam         |
|gd7      |H3K27ac    |rep2      |Yes          |gd7_K27ac_rep2_chr2L.dedup.bam   |
|tl10b    |input      |          |No           |tl10b_input_chr2L.bam            |
|tl10b    |input      |          |Yes          |tl10b_input_chr2L.dedup.bam      |
|tl10b    |H3K27ac    |rep1      |No           |tl10b_K27ac_rep1_chr2L.bam       |
|tl10b    |H3K27ac    |rep1      |Yes          |tl10b_K27ac_rep1_chr2L.dedup.bam |
|tl10b    |H3K27ac    |rep2      |No           |tl10b_K27ac_rep2_chr2L.bam       |
|tl10b    |H3K27ac    |rep2      |Yes          |tl10b_K27ac_rep2_chr2L.dedup.bam |

The first step in the differential binding analysis is identifying the regions
to be tested. It could be regions that you are interested in (e.g. promoters,
enhancers etc.), but in most cases we test differential binding at the peaks
called in each dataset to be compared. In order to make sure that the peaks are
not PCR artefacts, we use the deduplicated files for peak calling.

First, go to the folder where the bam files are stored and create the output folder.

```
cd ~/Desktop/Course_Materials/Practicals/Differential_Binding/bam
mkdir dedup_peaks
```

Call peaks for each sample using MACS2. H3K27ac marks are relatively sharp,
therefore we use the standard macs2 callpeak function as you did yesterday. As
a bonus exercise you could try the --broad option later, which is specifically
designed for broad histone modifications and joins smaller peaks into wider
peaks. For effective genome size we use 90% of chr2L, which is 21Mb.


Run the following to call peaks for the **gd7** cells:

```
# gd7 cells
macs2 callpeak -t gd7_K27ac_rep1_chr2L.dedup.bam \
    -c gd7_input1_chr2L.dedup.bam \
    --format BAM \
    --name dedup_peaks/gd7_H3K27ac_rep1 \
    --gsize 21000000 \
    --qvalue 0.01 \
    --call-summits

macs2 callpeak -t gd7_K27ac_rep2_chr2L.dedup.bam \
    -c gd7_input2_chr2L.dedup.bam \
    --format BAM \
    --name dedup_peaks/gd7_H3K27ac_rep2 \
    --gsize 21000000 \
    --qvalue 0.01 \
    --call-summits
```

Run the following to call peaks for replicate 1 of the **tl10b** cells:

```
macs2 callpeak -t tl10b_K27ac_rep1_chr2L.dedup.bam \
    -c tl10b_input_chr2L.dedup.bam \
    --format BAM \
    --name dedup_peaks/tl10b_H3K27ac_rep1 \
    --gsize 21000000 \
    --qvalue 0.01 \
    --call-summits
```

Now modify the last bit of code to the call peaks for replicate 2 of the
**tl10b** cells. 

The next steps will be performed within RStudio, please open it and set the
working directory to `~/Desktop/Course_Materials/Practicals/DiffBind/bam`. For
the DiffBind analysis we need to prepare a comma separated file containing
information about the samples to compare. We have prepared this file for you
and you can view it by loading it in R.

```
setwd("~/Desktop/Course_Materials/Practicals/Differential_Binding/bam") 
read_csv("H3K27ac_diffbind2.csv")
```

### Questions

1. Back in the Terminal, have a look at the output files obtained using macs2,
   how many peaks per sample did you obtain?  
2. In RStudio, which column in the file H3K27ac_diffbind2.csv we use to
   contrast the two cell types?

## Differential binding analysis using DiffBind

As the immunoprecipitation largely enriches a small fraction of the genome, it
is expected that in the ChIP samples we see exact duplicates, and removing
those can reduce your power to identify significant binding differences, as it
caps the dynamic range. In order to see the real level of change, we use the
original bam files including duplicates (before deduplication). 

**Note**: a higher than expected level of duplication indicates that there are
extensive PCR amplification artefacts. In this case it is advisable to remove
duplicates.

In RStudio, first load the DiffBind library and create a DiffBind object by reading the
H3K27ac_diffbind2.csv description file.

```
library(DiffBind)

H3K27ac <- dba(sampleSheet = "H3K27ac_diffbind2.csv")
```

Then, check the correlation between the peak sets (occupancy analysis).

```
plot(H3K27ac)
```

This plot only considers the read counts globally, and not the read counts at
the sample peaks. In the next step, count the reads falling in peak regions
using the dba.count function. With the minOverlap option we can adjust the
number samples where the analysed peaks should be present. By setting it to 2
we consider peaks that are shared in two biological replicates.

```
H3K27ac <- dba.count(H3K27ac, minOverlap=2)
print(H3K27ac)
```

The last step might take a few seconds. Next, we specify the factor that
separates our samples into groups for comparison. The option minMembers sets
the number of replicates per group required in the analysis.

```
H3K27ac <- dba.contrast(H3K27ac, categories=DBA_CONDITION, minMembers=2)
print(H3K27ac)
```

We now use the DESeq2 method to identify the differentially acetylated regions.

```
H3K27ac <- dba.analyze(H3K27ac, method=c(DBA_DESEQ2))
```

Then extract the differential regions:

```
H3K27ac.DB <- dba.report(H3K27ac)
H3K27ac.DB
```

Finally we create plots depicting the differences in the differentially
acetylated regions.

```
dba.plotMA(H3K27ac)
```

The MA plot shows the average normalised read counts on the X-axis and the
log-fold change on the Y-axis. This plot can be used to assess the effect of
the normalisation on the data as well as highlighting the significant
differentially bound regions (in red).

The PCA (principle component analysis) plot shows how the analysed samples
cluster according to the normalised read counts at the significantly
differential regions (FDR < 0.05):

```
dba.plotPCA(H3K27ac, contrast=1, label=DBA_CONDITION)
```

Boxplots are useful to view how read distributions differ between classes of
regions.

```
dba.plotBox(H3K27ac)
```

The first two boxes show normalised read counts in all regions in gd7 and
tl10b. The following two pairs of boxes show regions of significant increased
affinity in gd7 or tl10b respectively. 

**Note**: if there are differences in all regions (first two boxes), then
either the normalisation is not good or the regions selected are specific to
one of the conditions only.

Another plot that gives an idea about how different the identified regions are
is the heatmap, which shows significantly differential regions as rows and the
samples as columns. The samples are clustered according to the normalised read
counts in these regions:

```
dba.plotHeatmap(H3K27ac, contrast=1, correlations=FALSE)
```

### Questions: understanding DiffBind

1. How many differential regions you obtained?  
2. Is the experimental design appropriate for differential binding analysis?
3. Do the significantly differential regions show consistent differences
   between the two cell types?

CONGRATULATIONS! :) You’ve made it to the end of the practical.

We hope you enjoyed it! Don’t hesitate to ask any questions and feel free to
contact us any time.

## Bonus Exercise 1

Within DiffBind there is the possibility to use edgeR as the method for
differential binding analysis (DBA_EDGER). Repeat the analysis EdgeR instead of
DESeq2 in the analysis.

### Questions

1. How do the results from DESeq2 and edgeR compare? Which identified more
   regions? Which set of identified differential regions segregate the two cell
   types better in PCA?  
2. Do both methods use acceptable normalisation?

## Bonus Exercise 2

Given the high rate of PCR replicates in some of our samples it would be
advisable to remove duplicates from the BAM files used in the differential
analysis. Run the analysis without duplicates. Hint: We have already removed
the duplicates for the MACS2 peak calling, therefore you only need to change
the name of the path to the bam files in the DiffBind samplesheet.

### Questions

1. Do the samples cluster better with or without duplicates?  
2. How does removing the duplicates affect the observed level of change?

## Bonus Exercise 3

As H3K27ac peaks are broader than most transcription factor peaks, it would be
worthwhile to compare the results we obtained with using standard callpeaks in
MACS2 with those you would get if the peak regions were called with the –broad
option. Call broad peaks and use those as the regions tested for differential
acetylation. Note: broadPeaks file format should still be indicated as
narrowPeaks in the sample sheet, as the columns are the same, and there is no
separate option within DiffBind for broadPeaks.

### Question

1. How do the results compare in terms of number of differential regions and
   clustering of biological replicates?
