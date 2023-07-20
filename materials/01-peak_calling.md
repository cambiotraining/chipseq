---
title: Peak Calling
---

::: {.callout-tip}
#### Learning Objectives

- Apply the `nf-core/chipseq` pipeline for pre-processing, mapping, peak calling and QC of ChIP-seq data.
- Evaluate the quality of the ChIP from the MultiQC report.

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
- **Coverage tracks**: these files are in "bigwig" format, a compact representation of the genome coverage of each sample and can be used in genome browsers such as IGV (more on this below). 
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

In terms of mapping quality, we want to look out for: 

- What fraction of reads mapped to the genome?
- How many reads are mapped to the genome _after de-duplication_?
  ENCODE's guidelines suggest having 10-20M reads for narrow peaks and 40-50M reads for broader peaks. 
- What is the % of duplicated reads?
- What was our library complexity ([complexity curve](https://www.nature.com/articles/nmeth.2375))?


### ChIP Enrichment

This part of the QC report is more specific for ChIP-seq analysis. 
There are several quality metrics that can be used to assess the quality of the ChIP enrichment: 

- Strand-shift correlation and the derived metrics NSC and RSC
- [FRiP score](https://www.encodeproject.org/data-standards/terms/)
- [Fingerprint plot](https://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html#background)
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

![Screenshot of the IGV program, showing the coverage (normalised per 1M reads) for 3 replicate samples. In this region the peaks seem consistent across replicates, although the bottom one seems to have overal lower coverage. This is something we could investigate further.](images/igv_screenshot.png)


## Exercises

:::{.callout-note}
#### Running `nf-core/chipseq`

In the course data, you will find a shell script named `scripts/01-chipseq_workflow.sh`. 
This contains some code to run the Nextflow pipeline on our samples. 
However, there's a few things that need fixing first: 

- Fix the `samplesheet.csv` file where the word "FIXME" appears.
  You can open this file in a spreadsheet program to fix the issue. 
- Open the script `scripts/01-chipseq_workflow.sh` using a text editor (you can use `nano` from the command line, or double-click the file to open it with a graphical text editor). 
  Fix the code where the word "FIXME" appears. 
  Output the results to a directory called `results/nf-chipseq`.
- Run the script to execute the workflow. Hint: to run a shell script you use the program `bash name_of_script.sh`
  - Check that the workflow starts running successfully. You should get some progress of the different steps printed on the terminal.
  - The workflow will take some time to run (~15 minutes). You don't need to wait for it to finish before moving on to the next exercise.

:::{.callout-answer collapse=true}

To fix the samplesheet we needed to add the correct input sample name to E2-treated samples, which in this case was called "mcf7_input_e2". 

The fixed code in the script is: 

```bash
nextflow run nf-core/chipseq -profile singularity \
  --max_memory 23.GB --max_cpus 8 \
  --input samplesheet.csv \
  --outdir results/nf-chip \
  --fasta $(pwd)/resources/GRCh38.109.chr21.fasta.gz \
  --gtf $(pwd)/resources/GRCh38.109.chr21.gtf.gz \
  --blacklist $(pwd)/resources/ENCFF356LFX_exclusion_lists.chr21.bed.gz \
  --macs_gsize 40000000 \
  --min_reps_consensus 2
```

After saving the fixed script, we ran with with:

```bash
bash scripts/01-chip_workflow.sh
```

The script started running successfully, with the progress being printed on the screen. 

:::

:::


:::{.callout-exercise}
#### MultiQC report

In the previous exercise we've only processed a small number of samples, with fewer reads and one chromosome, due to time constraints. 
However, the full data includes 14 samples: 

- 6 samples pulled down with BRD4 antibody: 3 control (vehicle, 'veh') and 3 estrogen-treated ('e2');
- 6 samples pulled down with H2Bub1 antibody: 3 control (vehicle, 'veh') and 3 estrogen-treated ('e2');
- 2 input controls for 'veh' and 'e2' conditions.

We've processed the full dataset through the exact same `nf-core/chipseq` workflow and the output is in the `preprocessed/nf-chip` folder. 

Use the file browser to open the MultiQC report generated by the pipeline, which you can find in `preprocessed/nf-chipseq/multiqc/broadPeak/multiqc_report.html`.
Answer the following questions:

1. Was the quality of the raw FASTQ files acceptable? Was there any evidence for Illumina adapter contamination? Was this solved after trimming? (section "LIB: FastQC (raw)")
1. What was the percentage of mapped reads? (section "MERGED LIB: SAMTools (unfiltered)")
1. Were there any samples with >20% duplicate reads? (section "MERGED LIB: Picard (unfiltered)")
1. Looking at the fingerprint plot, what can you tell about the differences between the two antibodies? Do you expect different types of peaks? (section "MERGED LIB: deepTools")
1. From the read distribution profiles, do these two proteins bind to similar regions of annotated genes? (section "MERGED LIB: deepTools")
1. Was the peak count similar between the two antibodies? How about between replicates and treatments? (section "MERGED LIB: MACS2 peak count")
1. The ENCODE project uses the criteria below to evaluate if samples are good quality. Do our samples pass these thresholds? (section "MERGED LIB: MACS2 FRiP score" onwards)
  - FRiP > 1%
  - NSC > 1.05
  - RSC > 0.8

:::{.callout-answer collapse=true}

1. Yes, the quality was very high, as is often the case with current Illumina chemistry.
2. Generally the % of mapped reads was above 90%. The mapping rates seemed to be slightly lower for H2Bub1, but it doesn't seem to be a reason for concern. 
3. No, the lowest sample had 87% unique reads ("h2bub1_veh_rep3"), suggesting no reason for concern. 
4. The fingerprint plot suggests that H2Bub1 has more consistently localised peaks compared to BRD4, as it shows "sharper" fingerprint curves.
5. No, it seems like BRD4 binds primary upstream of the transcription start site (TSS), whereas H2Bub1 seems to occupy gene bodies (with some bias towards the 5' of the genes). 
6. No, H2Bub1 had substantially more peaks than BRD4. This histone modification marks transcribed genes, so perhaps it is no surprise that so many peaks are found for it.  
7. None of our samples have FRiP > 1, but this is perhaps to be expected as these do not seem to occur as narrow, sharp peaks, but rather more broad. The same applies for the other criteria, as they are generally more suited to use with sharp peaks. 

:::
:::

:::{.callout-exercise}
#### Visualising results in IGV

To visually inspect your results: 

- Open IGV (from the toolbar on the left)
- Choose the "Human (GRCh39/h19)" genome on the top-right drop-down menu
- Go to "File" -> "Load from file..." and then navigate to `preprocessed/nf-chipseq/bwa/mergedLibrary/bigwig/*.bigWig`
- Choose all the "bigWig" files in this folder (you can use the <kbd>Shift</kbd> key to select multiple files) and click "Open"
- On the search bar at the top search for "ACTB" to look for enrichment in this region. 

IGV is reasonably user-friendly and has many options to customise the display. 
By right-clicking on the track names on the left you can see many options. 

:::

