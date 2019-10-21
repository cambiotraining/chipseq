# ChIPseq practical

##### 3rd May 2019, Cambridge, UK

## Resources used

Bowtie: [http://bowtie-bio.sourceforge.net/index.shtml](http://bowtie-bio.sourceforge.net/index.shtml)

Samtools: [http://samtools.sourceforge.net/](http://samtools.sourceforge.net/)

Bedtools: [http://bedtools.readthedocs.io/en/latest/](http://bedtools.readthedocs.io/en/latest/)

UCSC tools: [http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)

IGV genome browser: [http://www.broadinstitute.org/igv/](http://www.broadinstitute.org/igv/)

MACS2: [https://github.com/taoliu/MACS](https://github.com/taoliu/MACS)

ChIPseeker: [https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)

MEME: [http://meme.sdsc.edu/meme/cgi-bin/meme.cgi](http://meme.sdsc.edu/meme/cgi-bin/meme.cgi)

TOMTOM: [http://meme.sdsc.edu/meme/cgi-bin/tomtom.cgi](http://meme.sdsc.edu/meme/cgi-bin/tomtom.cgi)

DeepTools: [http://deeptools.readthedocs.io/en/latest/index.html](http://deeptools.readthedocs.io/en/latest/index.html)

## Additional resources:

Ensembl: [http://www.ensembl.org](http://www.ensembl.org)

Original Data from: 
[http://www.ebi.ac.uk/ena/data/view/PRJNA242533](http://www.ebi.ac.uk/ena/data/view/PRJNA242533)

These data are reported in  Buecker, C. et al. Reorganization of Enhancer
Patterns in Transition from Naive to Primed Pluripotency. Cell Stem Cell 14,
838–853 (2014).

# Introduction

The goal of this hands-on session is to perform some basic tasks in the
analysis of ChIP-seq data. The first step includes an short read alignment for
a small subset of raw reads. We will align raw sequencing data to the mouse
genome using Bowtie and then we will manipulate the SAM output in order to
visualise the alignment on the *IGV* browser. Then based on these aligned reads
we will find enriched regions using the peak caller *MACS*. We will then
perform functional annotation and motif analysis on the predicted binding
regions.

# Prepare the environment

We will use one data set in this practical, which can be found in the ChIP-seq
directory on your desktop. Throughout this practical we will try to identify
potential transcription factor binding sites of Oct4 in mouse embryonic stem
cells.

Open the Terminal.

First, go to the right folder, where the data are stored.

```bash
cd ~/Course_Materials/Practicals/ChIPseq
```

The .fastq file that we will align is called Oct4.fastq. This file is based on
Oct4 ChIP-seq data published by Buecker et al. (2014). We will align these
reads to the mouse chromosome 1.

# Alignment 

There are a number of competing tools for short read alignment, each with its
own set of strengths, weaknesses, and caveats. Here we will try *Bowtie*, a
widely used ultrafast, memory efficient short read aligner.


*Bowtie* has a number of parameters in order to perform the alignment.

To view them all type

```
bowtie --help
```

Bowtie uses indexed genome for the alignment in order to keep its memory
footprint small. Because of time constraints we will build the index only for
one chromosome of the mouse genome. For this we need the chromosome sequence in
fasta format. This is stored in a file named mm10, under the subdirectory
`bowtie_index`.

The indexed chromosome is generated using the command (you need to uncompressed
the fasta file first):

```
gunzip bowtie_index/mm10.fa.gz
bowtie-build bowtie_index/mm10.fa bowtie_index/mm10
```

This command will output 6 files that constitute the index. These files that
have the prefix mm10 are stored in the bowtie_index subdirectory. 

To view if these files have been successfully created type:

```
ls -l bowtie_index
```

Now that the genome is indexed we can move on to the actual alignment. The
first two arguments make sure that the output is in SAM format using the `-S`
parameter and that Bowtie reports only uniquely mapped reads using the `-m 1`
option. The following argument for bowtie is the basename of the index for the
genome to be searched; in our case is `mm10`. The last argument is the name of
the fastq file.


Align the Oct4 reads using *Bowtie*:

```
gunzip fastq/Oct4.fastq.gz
bowtie -m 1 -S bowtie_index/mm10 fastq/Oct4.fastq > Oct4.sam
```

The above command outputs the alignment in SAM format and stores them in the
file `Oct4.sam`.

In general before you run *Bowtie*, you have to know which fastq format you
have. The available fastq formats in bowtie are:

```
--phred33-quals input quals are Phred+33 (default, same as --solexa-quals)
--phred64-quals input quals are Phred+64 (same as solexa1.3-quals)
--solexa-quals  input quals are from GA Pipeline version < 1.3
--solexa1.3-quals input quals are from GA Pipeline version >= 1.3
--integer-quals qualities are given as space-separated integers (not ASCII)
```

The fastq files we are working on is of Sanger format (Phred+33), which is the
default for *Bowtie*.

*Bowtie* will take 2-3 minutes to align the file. This is fast compared to
other aligners that sacrifice some speed to obtain higher sensitivity.

Look at the SAM format by typing:

```
head -n 10 Oct4.sam
```

Look at https://samtools.github.io/hts-specs/SAMv1.pdf to get more imformation
about SAM format.

## Questions

1. Can you distinguish between the header of the SAM format and the actual alignments?

2. What kind of information does the header provide you with?

3. To which chromosome are the reads mapped?

# Manipulate SAM output

SAM files are rather big and when dealing with a high volume of NGS data,
storage space can become an issue. We can convert SAM to BAM files (their
binary equivalent files that are not human readable) that occupy much less
space.


Convert SAM to BAM using samtools and store the output in the file `Oct4.bam`.
You have to instruct samtools that the input is in SAM format (`-S`{.ngs}), the
output should be in BAM format (`-b`{.ngs}) and that you want the output to be
stored in the file specified by the `-o`{.ngs} option:

```
samtools view -bSo Oct4.bam Oct4.sam
```

# Visualise alignments in IGV



It is often instructive to look at your data in a genome browser. Here, we use
IGV, a stand-alone browser, which has the advantage of being installed locally
and providing fast access. Web-based genome browsers, like Ensembl or the UCSC
browser, are slower, but provide more functionality. They do not only allow for
more polished and flexible visualisation, but also provide easy access to a
wealth of annotations and external data sources. This makes it straightforward
to relate your data with information about repeat regions, known genes,
epigenetic features or areas of cross-species conservation, to name just a few.
As such, they are useful tools for exploratory analysis.

Visualisation will allow you to get a 'feel' for the data, as well as detecting
abnormalities and problems. Also, exploring the data in such a way may give you
ideas for further analyses.

IGV is a stand-alone genome browser. Please check their website
[http://www.broadinstitute.org/igv/](http://www.broadinstitute.org/igv/) for
all the formats that IGV can display. For our visualisation purposes we will
use the BAM and bigWig formats.


When uploading a BAM file into the genome browser, the browser will look for
the index of the BAM file in the same folder where the BAM files is. The index
file should have the same name as the BAM file and the suffix `.bai`. Finally,
to create the index of a BAM file you need to make sure that the file is sorted
according to chromosomal coordinates. 


Sort alignments according to chromosome position and store the result in the
file called `Oct4.sorted.bam`:

```
samtools sort -o Oct4.sorted.bam Oct4.bam
```

Index the sorted file.

```
samtools index Oct4.sorted.bam
```

The indexing will create a file called `Oct4.sorted.bam.bai`. Note that you do
no??t have to specify the name of the index file when running samtools.


Another way to visualize the alignments is to convert the BAM file into a
bigWig file. The bigWig format is for display of dense, continuous data and the
data will be displayed as a graph. The resulting bigWig files are in an indexed
binary format.



The BAM to bigWig conversion takes place in two steps. Firstly, we convert the
BAM file into a bedgraph, called `Oct4.bedgraph`, using the tool
genomeCoverageBed from BEDTools:

``` 
genomeCoverageBed -bg -ibam Oct4.sorted.bam -g bowtie_index/mouse.mm10.genome >
Oct4.bedgraph
```

Then we convert the bedgraph into a binary graph, called `Oct4.bw`, using the
tool bedGraphToBigWig from the UCSC tools:

```
bedGraphToBigWig Oct4.bedgraph bowtie_index/mouse.mm10.genome Oct4.bw
```

Both of the commands above take as input a file called `mouse.mm10.genome` that
is stored under the subdirectory bowtie_index. These genome files are
tab-delimited and describe the size of the chromosomes for the organism of
interest. When using the UCSC Genome Browser, Ensembl, or Galaxy, you typically
indicate which species/genome build you are working. The way you do this for
BEDTools is to create a "genome" file, which simply lists the names of the
chromosomes (or scaffolds, etc.) and their size (in basepairs).

BEDTools includes pre-defined genome files for human and mouse in the
`/genomes` directory included in the BEDTools distribution.


Now we will load the data into the IGV browser for visualization. In order to
launch IGV type the following on your terminal:

```
igv &
```

On the top left of your screen choose from the drop down menu **Mus musculus
(mm10)**. Then in order to load the desire files go to:

**File --> Load from File**

On the pop up window navigate to **Desktop > ChIP-seq** folder and select the
file `Oct4.sorted.bam`.

Repeat these steps in order to load `Oct4.bw` as well.

Select chr1 from the drop down menu on the top left. Right click on the name of
Oct4.bw and choose Maximum under the Windowing Function. Right click again and
select Autoscale.

In order to see the aligned reads of the BAM file, you need to zoom in to a
specific region.


## Questions



1. Look for gene Lemd1 in the search box. Can you see an Oct4 binding site in
   the Lemd1 gene?

2. Using the '+' button on the top right zoom in more to see the details of the
   alignment. What is the main difference between the visualization of BAM and

3. What do you think the different colors mean?

# Alignment of control .fastq file

In the `ChIP-seq` folder you will find another .fastq file called
`Input.fastq`. Follow the steps described above for this dataset in order to
align the control reads to the mouse genome as well.

**This step is essential to be able to perform the rest of the analyses.**

# Finding enriched areas using MACS


MACS stands for Model based analysis of ChIP-seq. It was designed for
identifying transcription factor binding sites. MACS captures the influence of
genome complexity to evaluate the significance of enriched ChIP regions, and
improves the spatial resolution of binding sites through combining the
information of both sequencing tag position and orientation. MACS can be easily
used for ChIP-Seq data alone, or with a control sample to increase specificity.

Consult the MACS help file to see the options and parameters.

``` 
macs2 --help macs2 callpeak --help 
```

The input for MACS can be in ELAND, BED, SAM, BAM or BOWTIE formats (you just
have to set the --format flag). Options that you will have to use include:

``` 
-t to indicate the input ChIP file  
-c to indicate the name of the control file  
--format  the tag file format. If this option is not set MACS automatically  
detects which format the file is.   
--name to set the name of the output files  
--gsize This is the mappable genome size. With the read length we have, 70% of  
the genome is a fair estimation. Since in this analysis we include only reads  
from chromosome 1, we will use as gsize 70% of the length of chromosome 1 (197  
Mb). MACS also offers  shortcuts for human, 'mm' for mouse 'ce' for C. elegans  
and 'dm' for fruitfly.   
--call-summits when this option is set MACS detects all subpeaks in each
enriched region and returns their summits
--pvalue the P-value cutoff for peak detection. 
```


Now run macs using the following command:

``` 
macs2 callpeak -t Oct4.bam  -c Input.bam --format BAM --name Oct4 \
    --gsize 138000000 --pvalue 1e-3 --call-summits
```

Open the Excel peak file and view the peak details. Note that the number of
tags (column 6 `pileup`) refers to the number of reads at the summit (i.e.
summit height).

MACS generates its peak files in a file format called .narrowPeak file. This is
a simple text format containing genomic locations, specified by chromosome,
begin and end positions, and information on the statistical significance of
each called peak.

See

[http://genome.ucsc.edu/FAQ/FAQformat.html#format12](http://genome.ucsc.edu/FAQ/FAQformat.html#format12)

for details.

NarrowPeak files can also be uploaded to IGV or other genome browsers.


Try uploading the peak file generated by MACS to IGV. Find the first peak in
the file (use the head command to view the beginning of the `.narrowPeak`
file).



## Questions
1. Is the first peak that was called convincing to you?

2. Compare the Oct4 peak profile with that of H3K4me3 histone modification. To
   do so, upload the `H3K4me3.bw` file, which is located within the ChIP-seq
folder.

3. In IGV jump to location 1:35,987,885-36,101,395 for a sample peak. Do you
   think H3K4me3 peaks regions contain one or more modification sites? What
about the Oct4 peak located within the Lemd1 gene?

# Annotation: From peaks to biological interpretation


In order to biologically interpret the results of ChIP-seq experiments, it is
usually recommended to look at the genes and other annotated elements that are
located in proximity to the identified enriched regions. We will use R package
ChIPseeker
(https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html).
We need to modify the peak files into BED format, which is a tab-delimited file
that contains information on chromosome, start and end position for each
region. See http://genome.ucsc.edu/FAQ/FAQformat.html#format1 for details.



To convert `Oct4_peaks.narrowPeak` to a BED file,  type:

```
cut -f 1-4 Oct4_peaks.narrowPeak > Oct4_peaks.bed
```

Consult cut help file to see the options and parameters.

```
cut --help
```

Then start R in Terminal 

```
R
```

``` 
#Load requied R libraries 
library(dplyr)
library(ChIPseeker) 
library(TxDb.Mmusculus.UCSC.mm10.ensGene) 
library(clusterProfiler)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene 

#Read bed file 
peak <- readPeakFile("Oct4_peaks.bed")

#change chromosome annotation 
seqlevelsStyle(peak) <- "UCSC" 
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb)

#Peak localisation plot
plotAnnoPie(peakAnno)

#Annotate peaks 
peakAnno <- annotatePeak(peak, 
                         tssRegion = c(-1000, 1000),
                         flankDistance = 3000, 
                         TxDb=txdb)
peakTab <- as.data.frame(peakAnno)

#Add Gene Symbols 
geneSymbols <- data.frame(select(org.Mm.eg.db,
                                 keys = peakTab$geneId, 
                                 columns = c('ENSEMBL', 'SYMBOL'), 
                                 keytype = 'ENSEMBL'))
# This returns multiple results per an ENSEMBL ID. Let's just keep the first
# results for each.
geneSymbols <- geneSymbols[!duplicated(geneSymbols$SYMBOL),]

#Add to the peak table 
peakTab <- left_join(peakTab, geneSymbols, by=c(geneId="ENSEMBL"))

#Export the table 
write.table(peaktab, 
            file = "Oct4_annotated_peaks.txt", 
            row.names=FALSE, 
            sep="\t", 
            na="NA")

#Export a list of the closest genes

genes <- na.omit(unique(peakTab$SYMBOL)) 
write(genes, file = "Oct4_peaks_to_genes.txt")
```

Finally, quit R 
```
quit()
```

This list of closest downstream genes found under the link **The Full
Annotation File** can be the basis of further analysis. For instance, you could
look at the Gene Ontology terms associated with these genes to get an idea of
the biological processes that may be affected. Web-based tools like DAVID
(http://david.abcc.ncifcrf.gov) or GOstat (http://gostat.wehi.edu.au) take a
list of genes and return the enriched GO categories.

# Motif Analysis

It is often interesting to find out whether we can associate identified binding
sites with a sequence pattern or motif. We will use MEME for motif analysis.
The input for MEME should be a file in fasta format containing the sequences of
interest. In our case, these are the sequences of the identified peaks that
probably contain Oct4 binding sites.

Since many peak-finding tools merge overlapping areas of enrichment, the
resulting peaks tend to be much wider than the actual binding sites.
Sub-dividing the enriched areas by accurately partitioning enriched loci into a
finer-resolution set of individual binding sites, and fetching sequences from
the summit region where binding motifs are most likely to appear enhances the
quality of the motif analysis. Sub- peak summit sequences have already been
called by MACS with the `--call-summits` option.

De Novo motif finding programs take as input a set of sequences in which to
search for repeated short sequences. Since motif discovery is computationally
heavy, we will restrict our search for the Oct4 motif to the genome regions
around the summits of the 300 most significant Oct4 subpeaks on Chromosome 1.


Sort the Oct4 peaks by the height of the summit (the maximum number of
overlapping reads).

```
sort -k5 -nr Oct4_summits.bed > Oct4_summits.sorted.bed
```

Using the sorted file, select the top 300 peaks. 

```
head -n 300 Oct4_summits.sorted.bed > Oct4_top300_summits.sorted.bed
```

To create a BED file for the regions of 61 base pairs centred around the peak
summit, first download a file containing information about the chromosomes size
for the mouse genome. This file is required by 'bedtools slop' to add 30bp
around the summit peak.

```
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \ 
    -e "select chrom,size from mm10.chromInfo" > mm10.genome

bedtools slop -b 30 -i Oct4_top300_summits.sorted.bed -g mm10.genome \
    > Oct4_top300_summits.bed
```

Consult bedtools slop help file to see the options and parameters.

```
bedtools slop --help
```

The mouse genome sequence is available in FASTA format in the `bowtie_index`
directory. You can now use `bedtools` to extract the sequences around the Oct4
peak summits in FASTA format, which we save in a file named
`Oct4_top300_summits.fa`.

```
bedtools getfasta -fi bowtie_index/mm10.fa -bed Oct4_top300_summits.bed \
    -fo Oct4_top300_summits.fa
```

We are now ready to perform de novo motif discovery, for which we will use the
tool MEME.


Open a web bowser, go to the MEME website at
[http://meme-suite.org/](http://meme-suite.org/), and choose the 'MEME' tool.
Fill in the necessary details, such as:

 * the sub-peaks fasta file `Oct4_top300_summits.fa` (will need uploading), or
   just paste in the sequences.
 * the number of motifs we expect to find (1 per sequence)
 * the width of the desired motif (between 6 to 20) (under `Advanced options`)
 * the maximum number of motifs to find (3 by default). For Oct4 one classical
   motif is known.

Start Search. The results page will refresh automatically and once the tool has
finished running, please follow the link **"MEME html output"**

Scroll down until you see the first motif logo. We would like to know if this
motif is similar to any other known motif. We will use TOMTOM for this. Click
under the option **Submit/Download** and choose the TOMTOM button to compare to
known motifs in motif databases, and on the new page choose to compare your
motif to those in the JASPAR CORE and UniPROBE Mouse database.

## Questions

1. Which motif was found to be the most similar to your motif?

# Analysis centered around annotations 

You have seen until now how to perform peak-centered analyses of ChIP-seq data.
However, sometimes it is also of interest to analyse ChIP-seq profiles around
annotations (Transcription start sites, genes, transposable elemenst etc...),
without performing peak calling. This is often the case when analysing histone
modifications, in order to investigate the chromatin environment around gene
features and regulatory elements. The goal of this hands-on session is to plot
where the different histone modifications occur in relation to genes. 

We will be using ChIP-seq data for H3K27ac and H3K4me1 from Buecker C. et al.
(2014) and H3K4me3 from Marks H. et al. (2012). Although the former experiments
had 2 technical replicates, we will only be using the first replicate for
simplicity. For time purposes we have aligned the data for you and saved it in
the folder `Histone_modifications` of your `ChIP-seq` directory.

We will be using **deepTools** (http://deeptools.readthedocs.io/en/latest/), a
suite of tools for visualization, quality control and normalization of ChIP-seq
data.  We will start by looking at the genome-wide profiles of the selected
histone marks around TSS. A bed file of the TSS can be created using the UCSC
Table Browser at http://genome.ucsc.edu/cgi-bin/hgTables. We have prepared the
file for the mouse genome mm10 for you as TSS_mm10.bed in the ChIP-seq folder.

First, to go in the correct directory, type:

```
cd ./Histone_modifications
```

We will need to index the bam files we want to use:

```
samtools index H3K4me3.bam samtools index Input.bam
```

To analyse the ChIP-seq profile, the first step is to create a bigwig file from
the BAM file that is going to be used for plotting, normalised by the Input,
using bamCompare:

```
bamCompare --help

bamCompare -b1 H3K4me3.bam -b2 Input.bam -o Log2Ratio_H3K4me3_INPUT.bw
```

Then create a matrix, which is required for plotting the ChIP profiles, with
computeMatrix, using the mode to center profiles around TSS or TTS:

```
computeMatrix --help

computeMatrix reference-point -S Log2Ratio_H3K4me3_INPUT.bw -R ../TSS_mm10.bed \
     -a 500 -b 500 -out Matrix_log2_H3K4me3_allTSS.mat.gz
```

To plot profiles of ChIP-seq, use plotProfile:

```
plotProfile --help plotProfile -m Matrix_log2_H3K4me3_allTSS.mat.gz -out
Profile_log2_H3K4me3_allTSS.png
```

To plot hetmaps of the ChIP-seq signal at each TSS, use plotHeatmap.  Find how
to use plotHeatmap using the help function and generate a heatmap for H3K4me3

```
plotHeatmap --help
```

You can find a list of all deepTools and how to use them on
http://deeptools.readthedocs.io/en/latest/content/list_of_tools.html

**Repeat** these plots for another histone modification, H3K4me1.  For this use
the following functions: 'bamCompare', 'computeMatrix', 'plotProfile' and
'plotHeatmap'


## Questions

1. You can now view the H3K4me3 and H3K4me1 aggregation plots and heatmaps in
your folder. Where does the majority of the H3K4me3 and H3K4me1 signal appear
in relation to TSS?  

2. What proportion of genes display this pattern?  

3. Are the patterns you observe consistent with the putative functions of the
modifications in the figure and table given below? 

![histone_mods](../../Lectures/ChIPseq_Introduction/images/histone_mods.jpg)
![histone-mods-table](../../Lectures/ChIPseq_Introduction/images/histone-mods-table.jpg)

CONGRATULATIONS! You’ve made it to the end of the practical. 

Hope you enjoyed it!

Don’t hesitate to ask any questions and feel free to contact us any time (email
addresses on the front page).






# Bonus Exercise I: How to check if the ChIP worked





Before to explore any results form ChIP-seq data, it is important to verify if
the ChIP worked. In other words,  did the antibody-treatment enrich
sufficiently so that the ChIP signal can be separated from the background
signal?  To test this, we are using the 'plotFingerprint' function from
deepTools
(http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html).
This tool plots a profile of cumulative read coverages for the INPUT and ChIP.
It is used to determine how well the signal in the ChIP-seq sample can be
differentiated from the background distribution of reads in the control sample.
Also different type of profiles are expected for ChIPs resulting in narrow or
broad regions.

First, you should go back to the correct directory:

```
cd ~/Course_Materials/Practicals/ChIPseq
```




To test if the Oct4 ChIP worked, type:

```
plotFingerprint -b Oct4.sorted.bam Input.sorted.bam -plot
Oct4_fingerprint.png
```

## Questions



1. Compare 'Oct4_fingerprint.png' with the examples below. Can you conclude if
the ChIP worked?  \hrulefill\par \hrulefill\par

![fingerprint_examples](../../Lectures/ChIPseq_Introduction/images/fingerprint_examples.png)


# Bonus Exercise II: How to remove duplicates by Using Picard


Duplicate reads are the ones having the same start and end coordinates. This may be the result of technical duplication (too many PCR cycles), or over-sequencing (very high fold coverage). It is very important to put the duplication level in context of your experiment. For example, duplication level in targeted or re-sequencing projects may mean different than in RNA-seq experiments. In RNA-seq experiments over-sequencing is usually necessary when detecting the low expressed transcripts.
The duplication level computed by **FastQC** is based on sequence identity at the end of reads. Another tool, **Picard**, determines duplicates based on identical start and end positions.


Picard is a suite of tools for performing many common tasks with SAM/BAM format files. For more information see the Picard website and information about the various command-line tools available:

[http://picard.sourceforge.net](http://picard.sourceforge.net)

[http://picard.sourceforge.net/command-line-overview.shtml](http://picard.sourceforge.net/command-line-overview.shtml)

One of the Picard tools (MarkDuplicates) can be used to analyse and remove duplicates from the raw sequence data. The input for Picard is a sorted alignment file in `.bam` format. Short read aligners such as, bowtie, BWA, tophat etc. can be used to align fastq files against a reference genome to generate SAM/BAM alignment format.
However interested users can use the following general command to run the MarkDuplicates tool at their leisure and only need to provide a BAM file for the INPUT argument.



To remove duplicates reads for the Input, Type:

```
java -jar /applications/local/picard-tools/picard-tools/picard.jar MarkDuplicates \
  INPUT=Input.sorted.bam  VALIDATION_STRINGENCY=LENIENT OUTPUT=Input.sorted.picard.dup  \
  METRICS_FILE=Input.sorted.picard.metric  ASSUME_SORTED=true REMOVE_DUPLICATES=true
```

Changing the code, remove duplicates reads for the Oct4 ChIP.



## Questions
1. What is the percentage of duplication for the INPUT and Oct4 ChIP-seq data? You can find this information in 'Input.sorted.picard.matric' and 'Oct4.sorted.picard.matric' files.
\hrulefill\par


CONGRATULATIONS! You’ve made it to the end of the bonus exercises. 

Hope you enjoyed it!

Don’t hesitate to ask any questions and feel free to contact us any time (email addresses on the front page).


