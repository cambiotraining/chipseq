---
title: Peak Calling
---

::: {.callout-tip}
#### Learning Objectives

- Bulleted list of learning objectives
:::


## Peak Calling Workflow

One of the main steps in analysing ChIP-seq data is to identify regions of the genome enriched for sequencing reads, usually referred to as **peak calling**. 
These peaks are indicative of regions of the genome where our protein of interest binds to the DNA. 

There are several steps involved before we do the actual peak calling, namely filtering the sequencing reads for quality and aligning them to the reference genome. 
Because many of these steps are relatively standard, the bioinformatic community has built pipelines that can be used to automate these data processing steps in a scalable manner. 

We will use the Nextflow pipeline developed by the nf-core project, which we will refer as `nf-core/chipseq`. 
These pipelines are [very well documented](https://nf-co.re/chipseq/), with many options available to customise our analysis. 
The main advantage of these workflows is that they chain together dozens of tools, can process an arbitrary number of samples and can be run both locally and on HPC clusters. 

As part of their output they also provide an interactive quality control report (HTML file), which is very useful to identify any issues with our samples early on. 


## Running `nf-core/chipseq`

A typical command to run this workflow is given here: 

```bash
nextflow run nf-core/chipseq \
  -profile singularity \
  --input samplesheet.csv \
  --outdir path/to/results \
  --fasta path/to/genome.fasta.gz \
  --gff path/to/annotation.gff.gz \
  --blacklist path/to/exclusion_lists.bed.gz \
  --macs_gsize 2700000000 \
  --min_reps_consensus 2
```

Where: 

- `-profile singularity` is the mode we want to use for software management. 
  Singularity is recommended when running analysis on HPC clusters. Other alternatives include `conda` and `docker`. 
- `--input` is a CSV file containing information about our samples names and the directory paths to their respective FASTQ files. 
  The format of this file is fully detailed in the [documentation](https://nf-co.re/chipseq/2.0.0/docs/usage#samplesheet-input). 
- `--outdir` is the output directory where we want our results to be saved. This directory will be created if it does not already exist.
- `--fasta` and `--gff` are the directory paths to the reference genome and gene annotation, respectively. 
  These can be typically be downloaded from public repositories such as ENSEMBL.
- `--blacklist` is the directory path to a BED file containing regions of the genome to exclude from the analysis (regions identified as problematic for peak calling). 
- `--macs_gsize` is the estimated _mappable_ genome size to be used as input to the MACS software. 
  The MACS documentation recommends using 2.7e9 for the human genome and 1.87e9 for mouse. 
- `--min_reps_consensus` is the minimum number of replicates where a peak should be observed to be included in a "consensus" peak set.

When you start running the pipeline, you will get a progress log printed on the terminal and a message will be printed when it completes. 

By default the pipeline runs MACS in **broad peak** mode. 
If you want to run it in **narrow peak** mode (e.g. if you're studying transcription factors or narrow histone marks) you have to add the option `--narrow_peak` to your command. 


## Pipeline outputs

The `nf-core/chipseq` workflow generates many output files, which are explained in detail in the [workflow documentation](https://nf-co.re/chipseq/2.0.0/docs/output). 
We highlight some of the more relevant files we will use throughout our analysis: 

- **MultiQC report**: this is the quality control report, which is usually the first thing to look at once we run our pipeline. 
  This is located in `multiqc/<PEAK TYPE>/multiqc_report.html`.
- **Bigwig tracks**: these files provide a compact representation of the genome coverage of each sample and can be used in genome browsers such as IGV (more on this below). 
  You can find these files in `bwa/mergedLibrary/bigwig/*.bigWig`, with one file per sample.
  The coverage values in these files are normalised to 1M mapped reads. 
- **Called peaks**: files containing the peak intervals identified by the MACS software. 
  These files are in standard BED format (more details in the [MACS documentation page](https://macs3-project.github.io/MACS/docs/callpeak.html)). 
  These files are located in `bwa/mergedLibrary/macs2/<PEAK_TYPE>/*.broadPeak`.
- **Consensus peaks:** these are the BED files containing the consensus peak intervals identified in a minimum number of replicates set by the pipeline option mentioned above. 
  A consensus file is generated _per antibody_ and these can also be loaded into IGV. 
  These files are located in `bwa/mergedLibrary/macs2/<PEAK_TYPE>/consensus/<ANTIBODY>/*.bed`.

There are many more output files, which can be useful depending on the analysis you want to do. 
But these are the main ones we will focus on for now. 


## Quality report

The _MultiQC_ report generated by `nf-core/chipseq` contains a range of visualisations and metrics that can be used to assess the quality of our samples, from the sequencing reads, to the mapping all the way down to the peak calling. 

### Read quality

The first few sections of the report show:

- The results of the FastQC software (both before and after quality filtering), reporting average read quality, adapter contamination, GC content of our reads, amongst others. 
  See the [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) for full details. 
- A summary of the results from running the `cutadapt` software, which is used for quality-filtering of our reads and trimming adapter contaminations. 

Nowadays sequences tend to be very high quality, so it is unusual to have serious problems at this stage. 
However, if a very high fraction of your reads was filtered out, you should revisit what may have happened during library preparation and sequencing. 


### Mapping quality

TODO

One of the main considerations here is to check how many reads are mapped to the genome _after de-duplication_. 
In humans we usually want to have 10-20M reads for narrow peaks and 40-50M reads for broader peaks. 


### Peak Quality 

This is the more specific part of the QC report for ChIP-seq analysis. 
There are several quality metrics that can be used to assess the quality of the called peaks as well as the ChIP-seq experiment overal: 

- [Complexity curve](https://www.nature.com/articles/nmeth.2375)
- [Fingerprint plot](https://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html#background)
- [Associated fingerprint metrics](https://deeptools.readthedocs.io/en/latest/content/feature/plotFingerprint_QC_metrics.html)
- [FRiP score](https://www.encodeproject.org/data-standards/terms/)
- [phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools) - strand-shift correlation, NSC coefficient, RSC coefficient
- Correlation between replicates and PCA - in the pipeline this is done with DESeq2 (but `deeptools` can do this also)


<!-- https://www.encodeproject.org/data-standards/terms/#library -->


<!--
## Consensus Calling

The way the `nf-core/chipseq` pipeline seems to do the consensus calling is the following: 

- Identify peaks in at least `--min_reps_consensus` replicates per group (a group is a combination of antibody + control, I think)
- Then take the peaks across all the replicates for that antibody and take the _union_ of peaks

An alternative to do consensus calling more "manually" is to call peaks only if they are covered by at least X replicates.
In this case, the simplest way to do it is to find peaks across all the replicates of an antibody. 
This results in  more peaks than `nf-core/chipseq` gives, because it calls peaks if there was 1 rep in vehicle and 1 rep in e2 with a peak. 
Whereas `nf-core/chipseq` doesn't seem to consider that one as a peak (I checked this by loading the BED consensus from the two approaches on IGV). 
-->

## Visualising peaks in IGV

After thoroughly analysing our _MultiQC_ report, it is also advisable to visually inspect the coverage tracks for our samples (including input controls). 
This can give us a good indication of whether our peaks are sharp or broad (or somewhat in-between), whether they tend to occur over gene bodies, promoters, or spread across the genome. 

Although we can obtain this information from our _MultiQC_ report, it is still advisable to perform this visual inspection of the data. 
For example, you can look at genes/regions where you expect binding to occur for your protein of interest, or where you expect differences due to the treatment/conditions you used (e.g. from previous experiments, literature, qPCR, ...).

![TODO: IGV snapshot]()


## Exercise

:::{.callout-exercise}
# Running `nf-core/chipseq`

- Fix the `samplesheet.csv` file where the word "FIXME" appears.
- Fix the script `scripts/01-chipseq_workflow.sh` where the word "FIXME" appears. Output the results to a directory called `results/nf-chipseq`.
- Run the script to execute the workflow. Hint: to run a shell script you use the program `bash name_of_script.sh`
  - Check that the workflow starts running successfully. You should get some progress of the different steps printed on the terminal.
  - The workflow will take some time to run (~20 minutes). You don't need to wait for it to finish before moving on to the next exercise.

:::{.callout-answer collapse=true}


:::

:::


:::{.callout-exercise}
# MultiQC report

We've already processed the full dataset through the `nf-core/chipseq` workflow. 
The full data includes 14 samples: 

- 6 samples pulled down with BRD4 antibody: 3 control (vehicle, 'veh') and 3 estrogen-treated ('e2');
- 6 samples pulled down with H2Bub1 antibody: 3 control (vehicle, 'veh') and 3 estrogen-treated ('e2');
- 2 input controls for 'veh' and 'e2' conditions.

Use the file browser to open the MultiQC report generated by the pipeline, which you can find in `preprocessed/nf-chipseq/multiqc/broadPeak/multiqc_report.html`.
Answer the following questions:

- Was the quality of the raw FASTQ files acceptable? Was there any evidence for Illumina adapter contamination? Was this solved after trimming?
- According to the Encode guidelines do you think there were enough mapped reads for our analyses?
- Was there a high percentage of duplicate reads?
- Looking at the fingerprint plot, what can you tell about the differences between the two antibodies? Do you expect different types of peaks?
- From the read distribution profiles, do these two proteins bind to similar regions of annotated genes?
  - Can also look at the annotation from HOMER.
- Something about FRiP score
- Something about peak count
- From PCA do you think treatment had a significant effect on the binding profiles?

Check [ENCODE standards](https://www.encodeproject.org/chip-seq/histone-encode4/#standards).

:::{.callout-answer collapse=true}

For the PCA it's quite interesting to see a correlation between FRiP score and PC1 for BRD4.

:::
:::

:::{.callout-exercise}
#### Visualising results in IGV

- Choose right genome
- Change track height, colour, etc.

:::{.callout-answer collapse=true}

TODO

:::
:::

## Summary

::: {.callout-tip}
#### Key Points

- Last section of the page is a bulleted summary of the key points
:::
