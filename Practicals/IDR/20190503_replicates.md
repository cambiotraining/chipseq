
# Biological replicates practical

Sandra Cortijo sandra.cortijo@slcu.cam.ac.uk
Sergio Martinez Cuesta sermarcue@gmail.com
Ashley Sawle Ashley.Sawle@cruk.cam.ac.uk
Denis Seyres Denis.seyres@bioresource.nihr.ac.uk

(based on the practical by Myrto Kostadima)


## Introduction

In this practical, we will use ChIP-seq data for the PAX5 transcription factor, generated on the GM12878 lymphoblastoid cell line. BAM files for the two replicates and a matched control Input file have been downloaded from http://www.encodeproject.org, named PAX5_BR1TR1_ENCFF614SHU.bam, PAX5_BR2TR1_ENCFF475VQL.bam and Control.bam, respectively.



## Genome-wide correlation

We will first calculate genome-wide correlation between the read coverage of the two PAX5 ChIP-seq replicates, a basic measure of reproducibility across experiments.

Open the Terminal and go to the folder containing the data.

```bash
cd ~/Desktop/Course_Materials/Practicals/idr
```

The first step towards calculating pairwise correlation between the two biological replicates is to compute read coverages over genomic regions using the BAM files. The analysis can be performed for the entire genome by running the tool `multiBamSummary` in 'bins' mode. By default the coverage calculation is done for consecutive bins of 10 kilobases long.

```bash
multiBamSummary bins --bamfiles processed_data/bam/PAX5_BR1TR1_ENCFF614SHU.bam processed_data/bam/PAX5_BR2TR1_ENCFF475VQL.bam --labels PAX5_BR1TR1 PAX5_BR2TR1 -out processed_data/coverage/PAX5_read_counts.npz --outRawCounts processed_data/coverage/PAX5_read_counts.tab
```

Finally we compute the Pearson correlation of these two PAX5 ChIP-seq profiles, using the command `plotCorrelation`.

```bash
plotCorrelation -in processed_data/coverage/PAX5_read_counts.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Average Reads per Bin" --whatToPlot scatterplot -o processed_data/coverage/PAX5_pearson_corr.png --outFileCorMatrix processed_data/coverage/PAX5_pearson_corr.tab --removeOutliers
```


### Questions

- What is the correlation between these two samples?
- Is this higher or lower than you expected?



## Working with two biological replicates - CHANCE statistics

CHANCE (CHIP-seq ANalytics and Confidence Estimation) is a standalone package for ChIP-seq quality control and protocol optimization first published in [Diaz et al., Genome Biology (2012)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r98). The authors provide a user-friendly graphical software that among other functionalities it estimates the strength and quality of immunoprecipitations (CHANCE is available at https://github.com/songlab/chance). For the purposes of our practical we will calculate the enrichment of our two PAX5 biological replicates over our Control sample, not using the CHANCE tool itself, but using [deepTools](https://deeptools.readthedocs.io/en/latest/).

First, we need to make sure that all our BAM files, for both replicates and Control samples are indexed. Have a look at the `processed_data/bam` folder to see if all BAM files have their respective `.bam.bai` file. If that's not the case, create indexes for all the bam files using `samtools index`, e.g. `samtools index PAX5_BR1TR1_ENCFF614SHU.bam`.

Once all the index files have been created, run the following command to calculate the enrichment of the sample PAX5_BR1TR1_ENCFF614SHU over the Control using bin sizes of 1000bp.

```bash
plotFingerprint -b processed_data/bam/PAX5_BR1TR1_ENCFF614SHU.bam processed_data/bam/Control.bam --binSize 1000 --labels PAX5_BR1TR1 Control --JSDsample processed_data/bam/Control.bam --outQualityMetrics processed_data/bam/PAX5_BR1TR1_ENCFF614SHU_fingerprint.txt -plot processed_data/bam/PAX5_BR1TR1_ENCFF614SHU_fingerprint.png
```

For a detailed explanation of the output *_fingerprint.txt* file format, read [this](https://deeptools.readthedocs.io/en/latest/content/feature/plotFingerprint_QC_metrics.html). Now run the above command for the other replicate as well. Once the tool has finished, open both the PNG files created.


### Question

- Is the enrichment similar for both PAX5 biological replicates? If not, what differences do you observe?



## Working with two biological replicates - Irreproducible Discovery Rate (IDR) analysis

Peak calling using MACS2 was performed for both replicates individually and after merging the alignments from both replicates, e.g. using `samtools merge`, to create a pooled peak set. These datasets are already provided for you in the `processed_data/macs2` directory, they are labelled `PAX5_BR1TR1_peaks.bed`, `PAX5_BR2TR1_peaks.bed` and `PAX5_pooled_peaks.bed`, respectively.

The narrowPeak format in which the peaks are reported is an ENCODE format, using an extension of the BED format for providing peaks and associated scores and p-values. See [this](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) for full information.

The ChIP-seq experiment and alignment may be affected by read mapping artefacts caused by biases in chromatin accessibility and ambiguous alignments. Some of these spurious regions can be removed by filtering out any peaks that overlap a blacklist, believed to contain experiment and cell type independent areas of high artefactual signal. The ENCODE blacklist has been generated by combining regions of known repeats and manually curated genomic regions of ubiquitous open chromatin and input sequence signal. The file can be downloaded from [here](https://www.encodeproject.org/annotations/ENCSR636HFF/).

We can now remove any peaks which overlap with the ENCODE blacklist:

```bash
bedtools intersect -a processed_data/macs2/PAX5_BR1TR1_peaks.bed -b annotation/GRCh38_blacklisted_regions.bed -v > processed_data/macs2/PAX5_BR1TR1_peaks_filtered.bed
bedtools intersect -a processed_data/macs2/PAX5_BR2TR1_peaks.bed -b annotation/GRCh38_blacklisted_regions.bed -v > processed_data/macs2/PAX5_BR2TR1_peaks_filtered.bed
bedtools intersect -a processed_data/macs2/PAX5_pooled_peaks.bed -b annotation/GRCh38_blacklisted_regions.bed -v > processed_data/macs2/PAX5_pooled_peaks_filtered.bed
```

Next, we need to sort peaks in the files by significance (p-value). For `MACS2`, the p-value works best for ranking using the IDR procedure. See [this](https://sites.google.com/site/ anshulkundaje/projects/idr#TOC-Peak-callers-tested-with-IDR) for the best ranking measures to use with other peak callers.

Here, we sort by p-value and then use the 100,000 most significant peaks:

```bash
sort -k 8fr processed_data/macs2/PAX5_BR1TR1_peaks_filtered.bed | head -n 100000 > processed_data/macs2/PAX5_BR1TR1_top100000_peaks.bed
sort -k 8fr processed_data/macs2/PAX5_BR2TR1_peaks_filtered.bed | head -n 100000 > processed_data/macs2/PAX5_BR2TR1_top100000_peaks.bed
sort -k 8fr processed_data/macs2/PAX5_pooled_peaks_filtered.bed | head -n 100000 > processed_data/macs2/PAX5_pooled_top100000_peaks.bed
```

Finally, we can perform the IDR analysis on the peaks called in the two technical replicates, using the peak list from the merged replicates

```bash
idr --samples processed_data/macs2/PAX5_BR1TR1_top100000_peaks.bed processed_data/macs2/PAX5_BR2TR1_top100000_peaks.bed --peak-list processed_data/macs2/PAX5_pooled_top100000_peaks.bed --idr-threshold 0.05 --output-file idr/PAX5_replicates_idr.txt --plot
```


### Questions: understanding the IDR output plot

- How reproducible are the peaks called in these two biological replicates for PAX5 binding in the GM12878 cell line?
- What factors could lower the reproducibility between the two ChIP-seq experiments?

CONGRATULATIONS! :) You’ve made it to the end of the practical.



## Bonus Exercise 1

The IDR statistics can also be used to flag data sets with low reproducibility. This may be due to one of the two replicates being of lower ChIP enrichment, hence having a high signal-noise ratio. In this case, the standard IDR protocol would record few reproducible peaks, despite one replicate having high information content.

ENCODE has developed a rescue strategy in this case by using pseudo-replicates. These pseudo-replicates are generated by pooling all the reads, and then randomly splitting them into two files. These pseudo-replicates do not represent true biological or technical replicates, but attempt to model the stochastic noise in the sampling of sequenced reads from a population of fragments.

The pseudo-replicates analysis uses a lower IDR threshold than with biological/technical replicates, typically 0.0025, due to the reduced noise. Check the `idr` help in the Terminal and modifiy the command above to run the IDR analysis for pseudo-replicates of the GM12878 PAX5 ChIP-seq data.


### Question

- Does the IDR method select more peaks using the original technical replicates or the pseudo-replicates approach?
