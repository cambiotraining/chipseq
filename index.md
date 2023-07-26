---
title: "Analysis of ChIP-seq Data"
author: "Hugo Tavares"
date: today
number-sections: false
---

## Overview 

Chromatin immunoprecipitation followed by sequencing (ChIP-seq) is a method used to identify binding sites for transcription factors, histone modifications and other DNA-binding proteins across the genome. 
These materials cover the fundamentals of ChIP-seq data analysis, from raw data processing to downstream applications.  
We will start with an introduction to ChIP-seq methods, including important considerations when designing your experiments. 
We will cover the bioinformatic steps in a standard ChIP-seq analysis workflow, covering raw data quality control, trimming/filtering, mapping, duplicate removal, post-mapping quality control, peak calling and peak annotation. 
We will discuss metrics used for quality assessment of the called peaks when multiple replicates are available, as well as the analysis of differential binding across sample groups. 
Finally, we will also cover tools and packages that can be used for visualising and exploring your results. 

::: {.callout-tip}
### Learning Objectives

- Describe how ChIP-seq data is generated and what information it provides about the (epi)genome
- Recall the experimental design considerations that are needed when performing ChIP-seq experiments
- Understand the bioinformatic steps involved in processing ChIP-seq data
- Interpret and assess the quality of your data and results
- Perform differential binding analysis to compare different groups of samples 
:::


### Target Audience

This course is aimed at researchers with **no prior experience in the analysis of ChIP-seq data**, who would like to get started in processing their data using a standardised pipeline and perform downstream analysis and visualisation of their results.


### Prerequisites


- Basic understanding of high-throughput sequencing technologies.
  - Watch [this iBiology video](https://youtu.be/mI0Fo9kaWqo) for an excellent overview. 
- A working knowledge of the UNIX command line ([course registration page](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-unix2)).
  - If you are not able to attend this prerequisite course, please work through our [Unix command line materials](https://cambiotraining.github.io/unix-shell/) ahead of the course (up to section 7). 
- A working knowledge of R ([course registration page](https://training.csx.cam.ac.uk/bioinformatics/course/bioinfo-introRbio)).
  - If you are not able to attend this prerequisite course, please work through [our R materials](https://cambiotraining.github.io/intro-r/) ahead of the course.


<!-- Training Developer note: comment the following section out if you did not assign levels to your exercises
### Exercises

Exercises in these materials are labelled according to their level of difficulty:

| Level | Description |
| ----: | :---------- |
| {{< fa solid star >}} {{< fa regular star >}} {{< fa regular star >}} | Exercises in level 1 are simpler and designed to get you familiar with the concepts and syntax covered in the course. |
| {{< fa solid star >}} {{< fa solid star >}} {{< fa regular star >}} | Exercises in level 2 combine different concepts together and apply it to a given task. |
| {{< fa solid star >}} {{< fa solid star >}} {{< fa solid star >}} | Exercises in level 3 require going beyond the concepts and syntax introduced to solve new problems. |
-->

## Authors
<!-- 
The listing below shows an example of how you can give more details about yourself.
These examples include icons with links to GitHub and Orcid. 
-->

About the authors (alphabetical by surname):

- **Sandra Cortijo**
  <a href="https://orcid.org/0000-0003-3291-6729" target="_blank"><i class="fa-brands fa-orcid" style="color:#a6ce39"></i></a> 
  <a href="https://github.com/scortijo" target="_blank"><i class="fa-brands fa-github" style="color:#4078c0"></i></a>  
  _Affiliation_: Centre National de la Recherche Scientifique: Montpellier  
  _Roles_: writing; conceptualisation; coding
- **Sergio Martinez Cuesta**
  <a href="https://orcid.org/0000-0001-9806-2805" target="_blank"><i class="fa-brands fa-orcid" style="color:#a6ce39"></i></a> 
  <a href="https://github.com/semacu" target="_blank"><i class="fa-brands fa-github" style="color:#4078c0"></i></a>  
  _Affiliation_: AstraZeneca, Cambridge  
  _Roles_: writing; conceptualisation; coding
- **Sankari Nagarajan**
  <a href="https://orcid.org/0000-0001-8748-6223" target="_blank"><i class="fa-brands fa-orcid" style="color:#a6ce39"></i></a> 
  _Affiliation_: University of Manchester  
  _Roles_: writing; conceptualisation
- **Ashley Sawle**
  <a href="https://orcid.org/0000-0002-2985-5059" target="_blank"><i class="fa-brands fa-orcid" style="color:#a6ce39"></i></a> 
  <a href="https://github.com/AshKernow" target="_blank"><i class="fa-brands fa-github" style="color:#4078c0"></i></a>  
  _Affiliation_: Cancer Research UK, Cambridge Institute  
  _Roles_: writing; conceptualisation; coding  
- **Denis Seyres**
  <a href="https://orcid.org/0000-0002-2066-6980" target="_blank"><i class="fa-brands fa-orcid" style="color:#a6ce39"></i></a> 
  _Affiliation_: Universitätsspital Basel: Basel  
  _Roles_: writing; conceptualisation; coding
- **Hugo Tavares**
  <a href="https://orcid.org/0000-0001-9373-2726" target="_blank"><i class="fa-brands fa-orcid" style="color:#a6ce39"></i></a> 
  <a href="https://github.com/tavareshugo" target="_blank"><i class="fa-brands fa-github" style="color:#4078c0"></i></a>  
  _Affiliation_: Bioinformatics Training Facility, University of Cambridge  
  _Roles_: writing; conceptualisation; coding



## Citation

<!-- We can do this at the end -->

Please cite these materials if:

- You adapted or used any of them in your own teaching.
- These materials were useful for your research work. For example, you can cite us in the methods section of your paper: "We carried our analyses based on the recommendations in _Cortijo S et al. (2023)_.".

You can cite these materials as:

> Cortijo S, Martinez Cuesta S, Nagarajan S, Sawle A, Seyres D, Tavares H (2023) "cambiotraining/chipseq: Analysis of ChIP-seq Data", https://cambiotraining.github.io/chipseq/

Or in BibTeX format:

```
@Misc{,
  author = {Cortijo, Sandra AND Martinez Cuesta, Sergio AND Nagarajan, Sankari AND Sawle, Ashley AND Seyres, Denis AND Tavares, Hugo},
  title = {cambiotraining/chipseq: Analysis of ChIP-seq Data},
  month = {July},
  year = {2023},
  url = {https://cambiotraining.github.io/chipseq/}
}
```


## Acknowledgements

<!-- if there are no acknowledgements we can delete this section -->

There are many online resources that inspired our own materials (e.g. package vignettes) and we cite them where relevant. 

We also recommend the following training materials:

- [Understanding chromatin biology using high throughput sequencing](https://hbctraining.github.io/Intro-to-ChIPseq/schedule/2-day.html) from the Harvard Chan Bioinformatics Core
- [Introduction to ChIPseq using HPC](https://hbctraining.github.io/Intro-to-ChIPseq/schedule/2-day.html) from the Harvard Chan Bioinformatics Core
- [ChIP-seq analysis](https://www.bioinformatics.babraham.ac.uk/training.html#chip) from the Babraham Institute