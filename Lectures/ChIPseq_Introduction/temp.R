

### ANALYSIS OVERVIEW

```
Enriched
areas
```
```
Sequencing
```
```
Control
sample
```
```
Image
analysis Alignment Peak calling
```
```
Data
integration
```
```
Other
```
```
SubPeaks
```
```
Peak
annotation
```
```
Genome
browser
```
```
Motif
analysis
```
```
50 - 100bp
reads
```
```
Sequencing
Bioinformatics analysis
Downstream analysis
```

###### Reference	Sequence GCTGATGTGCCGCCTCACTTCGGTGG

###### CTGATGTGCCGCCTCACTTCGGTGGT

###### TGATGTGCCGCCTCACTACGGTGGTG

###### GATGTGCCGCCTCACTTCGGTGGTGA

###### GCTGATGTGCCGCCTCACTACGGTG

###### GCTGATGTGCCGCCTCACTACGGTG

**Short-reads**

### (Short) Read Alignment

###### GOAL: Given a reference sequence and a set of short

###### reads, align each read to the reference sequence


- Not all of the genome is ‘available’ for mapping
- Align your reads to the unmasked genome
- For ChIP-seq, usually short reads are used (50/100bp)
- Limited gain in using longer reads (again, unless you

###### have a specific interest in repeat regions)

### MAPPABILITY

```
*Calculated based on 30nt sequence tags
Rozowsky, (2009)
```

###### RNA

###### DNA

###### Methylatio

###### n

###### microRNA

[http://www.ebi.ac.uk/~nf/hts_mappers/](http://www.ebi.ac.uk/~nf/hts_mappers/)


### Reads can map in multiple locations

- Some parts of the genome will not be unique:
- Common, repeated motifs (proteins domains)
- Repeat regions

This can not always be resolved

##### >CHROMOSOME_1

##### GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAA

##### TCCATTTGTTCAACTCACATTAAATAGTCGATCAAAT

##### AGTTGGTTCAAAGCAGTCCATT

###### TTCAAAGC

###### ?


### And what do I do with duplicates?

- Reads that align in exactly the same place (same start +
    
    ###### same CIGAR string)
    
    - Duplicates can occur from:
    - Artefacts from sequencing (PCR artefacts)
  - Real biological signal
  - We cannot tell apart which one, unless we use barcodes.
  
  ```
  Parekh S et al., SciRep.2016 May 9;6:25533.
  ```
  
  ##### LET’S TALK ABOUT QUALITY CONTROL
  
  
  #### Read alignment statistics
  
  ```
  Sample1 Sample2 Sample3 Sample4 Sample5
  ```
  
  ### Did my ChIP experiment work?
  
  ```
  60% of the
  genome didn’t
  have sufficient
  coverage in the
  IP
  ```
  ```
  IP worked
  ```
  ```
  over 60% of the
  reads map to a
  small percentage
  of the genome,
  indicating
  amplification
  bias
  ```
  ```
  weak IP: IP and
  Input curves are
  not well
  separated
  ```
  Diaz et al. (2012)
  https://github.com/songlab/chance https://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html
  
  
  ### ANALYSIS OVERVIEW
  
  ```
  Enriched
  areas
  ```
  ```
  Sequencing
  ```
  ```
  Control
  sample
  ```
  ```
  Image
  analysis Alignment Peak calling
  ```
  ```
  Data
  integration
  ```
  ```
  Other
  ```
  ```
  SubPeaks
  ```
  ```
  Peak
  annotation
  ```
  ```
  Genome
  browser
  ```
  ```
  Motif
  analysis
  ```
  ```
  50 - 100bp
  reads
  ```
  ```
  Sequencing
  Bioinformatics analysis
  Downstream analysis
  ```
  
  ```
  Kharchenko(2008) Nature Biotechnology
  ```
  ### STRAND SPECIFIC PROFILE
  
  
  ### PEAK CALLING
  
  - Basic -regions are scored by the
  
  ```
  number of tags in a window of a
  given size. Then assess by
  enrichment over control and
  minimum tag density.
  ```
  - Advanced -take advantage of
  
  ```
  the directionality of the reads.
  Kharchenko(2008) Nature Biotechnology
  ```
  
  - Adjust for sequence alignability-regions that contain repetitive elements
  have different expected tag count
  - Different ChIP-seqapplications produce different type of peaks. Most current
  tools have been designed to detect sharp peaks (TF binding, histone
                                                  modifications at regulatory elements)
  - Alternative tools exist for broader peaks (histone modifications that mark
                                               domains -transcribed or repressed), e.g. SICER
  
  ### PEAK CALLING-CHALLENGES
  
  
  Park J, Nature Reviews Genetics, 2009
  
  
  ### MACS TOOL
  
  - Model the shift size between
  +/-strand tags
  - Scan the genome to find
  regions with tags more than
  mfold(between 10-30)
  enriched relative to random
  tag distribution
  - Randomly sample 1000 of
  these (high quality peaks)
  and calculate the distance
  between the modes of their
  +/-peaks
  - Shift all the tags by d/2
  toward the 3’ end.
  
  
  - x = observed read count
  - λ=expected read number
  - P = Probability to find a peak higher than x
  
  ## Is a peak greater than expected by chance?
  
  ## €
  
  ## P = 1 −
  
  ## e
  
  −λ
  
  ## λ
  
  _k_
  
  _k_ = (^0) _k_!
    _x_ − 1
  
  # ∑
  
  ###### €
  
  ###### λ =
  
  ###### ( read _ size )*( mapped _ reads )
  
  ###### Alignable _ genome _ size
  
  ## • Tag distribution along the genome can be modelled by
  
  ###### a Poisson distribution
  
  ```
  Slide adapted from http://www.slideshare.net/lucacozzuto/macs-course
  ```
  
  ```
  tag count = 2
  total reads = 30,000,000
  read length = 36
  mappable human genome = 2,700,000,000
  ```
  ## e.g. non –significant peak
  
  ###### €
  
  ###### λ BG =
  
  ###### ( 36 )*( 30 , 000 , 000 )
  
  ###### 2 , 700 , 000 , 000
  
  ###### = 4
  
  ## €
  
  ## P = 1 −
  
  ## e
  
  − 4
  
  ## * 4
  
  _k_
  
  _k_ = (^0) _k_!
    1
  
  # ∑ = 0.9
  
  ```
  Slide adapted from http://www.slideshare.net/lucacozzuto/macs-course
  ```
  
  ```
  tag count = 10
  total reads = 30,000,000
  read length = 36
  mappable human genome = 2, 700,000,000
  ```
  ## e.g. significant peak
  
  ###### €
  
  ###### λ BG =
  
  ###### ( 36 )*( 30 , 000 , 000 )
  
  ###### 2 , 700 , 000 , 000
  
  ###### = 4
  
  ## €
  
  ## P = 1 −
  
  ## e
  
  − 4
  
  ## * 4
  
  _k_
  
  _k_ = (^0) _k_!
    9
  
  # ∑ = 0.008
  
  ```
  Slide adapted from http://www.slideshare.net/lucacozzuto/macs-course
  ```
  
  - Candidate peaks are evaluated via
  comparing them against a “local”
  distribution.
  - Fold enrichment = Enrichment over the
  λlocal
  - False Discovery Rate (FDR) is
  calculated (#peaks in control) / (#peaks
    in IP) Control peaks are calculated by
  swapping control and sample.
  
  #### The benefit of having a control ChIP experiment
  
  ```
  Slide adapted from http://www.slideshare.net/lucacozzuto/macs-course
  ```
  
  - Remove duplicate tags (in excess of what can be expected by
                           
                           chance)
  
  - Slide window across the genome to find candidate peaks with a
  
  ```
  significant tag enrichment (Poisson distribution, global
                              background, p-value 10e-5)
  ```
  - Also looks at local background levels and eliminates peaks that
  
  are not significant with respect to local background
  
  - Uses the control sample to eliminate peaks that are also present
  
  there
  
  ### MACS PEAK DETECTION
  
  
  - To ensure that experiments are reproducible at least 2
  
  replicates must be performed.
  
  - To access replicate agreement, Irreproducible Discovery
  
  Rate analysis can be used:
    
    - If replicates measure the same biology high scores (e.g. low p-
                                                            values) represent strong evidence of being genuine signals and
  are ranked high and also more consistently on the replicates
  than those of spurious signals.
  - Broad histone modifications are hard to quantify for
  
  reproducibility
  
  ### Do we need replicates?
  
  
  - Plotting the consistency between a pair of rank lists (that contains both significant
                                                           and insignificant peaks) will indicate how many peaks have been reliably detected
  as the point showing change in consistency
  - IDR score for each signal: reflects the probability for the signal to belong to the
  irreproducible group.
  
  ### IDR ANALYSIS
  
  ###### reasonable
  
  ###### consistency
  
  ###### poor
  
  ###### consistency
  
  
  ### ANALYSIS OVERVIEW
  
  ```
  Enriched
  areas
  ```
  ```
  Sequencing
  ```
  ```
  Control
  sample
  ```
  ```
  Image
  analysis Alignment Peak calling
  ```
  ```
  Data
  integration
  ```
  ```
  Other
  ```
  ```
  SubPeaks
  ```
  ```
  Peak
  annotation
  ```
  ```
  Genome
  browser
  ```
  ```
  Motif
  analysis
  ```
  ```
  50 - 100bp
  reads
  ```
  ```
  Sequencing
  Bioinformatics analysis
  Downstream analysis
  ```
  
  #### ANALYSIS DOWNSTREAM TO PEAK CALLING
  
  - **Visualisation** - genome browser: Ensembl, UCSC, IGV
  
  
  #### ANALYSIS DOWNSTREAM TO PEAK CALLING
  
  - **Visualization** - genome browser: Ensembl, UCSC, IGV
  - **Peak Annotation** - finding interesting features surrounding peak regions:
    - PeakAnalyzer
  - ChIPpeakAnno(R package)
  - GREAT
  - bedtools
  - PAVIS
  
  
  #### ANALYSIS DOWNSTREAM TO PEAK CALLING
  
  - **Visualization** - genome browser: Ensembl, UCSC, IGV
  - **Peak Annotation** - finding interesting features surrounding peak regions
  - Correlation with **expression data**
    
    
    #### ANALYSIS DOWNSTREAM TO PEAK CALLING
    
    - **Visualization** - genome browser: Ensembl, UCSC, IGV
  - **Peak Annotation** - finding interesting features surrounding peak regions:
    - Correlation with expression data
  - Discovery of **binding sequence motifs**
    
    
    #### ANALYSIS DOWNSTREAM TO PEAK CALLING
    
    - **Visualization** - genome browser: Ensembl, UCSC, IGV
  - **Peak Annotation** - finding interesting features surrounding peak regions:
    - Correlation with expression data
  - **Discovery of binding sequence motifs**
    - **Gene Ontology analysis** on genes that bind the same factor or have the
  same modification
  
  
  #### ANALYSIS DOWNSTREAM TO PEAK CALLING
  
  - **Visualization** - genome browser: Ensembl, UCSC, IGV
  - **Peak Annotation** - finding interesting features surrounding peak regions:
    - Correlation with expression data
  - **Discovery of binding sequence motifs**
    - Gene Ontology analysis on genes that bind the same factor or have the
  same modification
  - Correlation with SNP data to find **allele-specific binding**
    
    
    ### NON-PEAK BASED ANALYSIS
    
    - Can be more appropriate for non-specific binding sites
  - Does not involve calling significant peaks, and
  
  ###### discarding the rest of the signal as noise
  
  - Instead, the whole signal is used for analysis
  - Global analysis, e.g. looking for enrichment at TSS
  - Tools include CEAS and EpiChip and deepTools
  
  
  ### NON-PEAK BASED ANALYSIS
  
  (^49) • Example: Profile and Heatmap for H3K4me3
  around the TSS
  
  
  #### “It’s an absolute myth that you can send an algorithm
  
  #### over raw data and have insights pop up”
  
  - Jeffrey Heer (University of Washington)
  
  
                                                                                    