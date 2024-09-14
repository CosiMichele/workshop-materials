# An Introduction to Building a Transcriptomics Pipeline: RNA-Seq Alignment and Analysis

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/wiki/griffithlab/rnaseq_tutorial/Images/RNA-seq_Flowchart2.png" width="650">
</p>
<br>

Image credits: [Griffith Lab](https://griffithlab.org/) ([GitHub](https://github.com/griffithlab), [source](https://github.com/griffithlab/rnaseq_tutorial/wiki/Indexing)).

---

>[!important]
> :clock1: **Schedule**
> - 3:00pm-3:05pm: Welcome and introdution to topic
> - 3:05pm-3:25pm: Alignment and Data Preparation
> - 3:25pm-3:45pm: Transcript Assembly and Quantification
> - 3:45pm-4:00pm: Exploring Differential Gene Expression and closing remarks

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Hands-on Experience with RNA-seq Data Processing
> - Understanding Transcript Assembly and Quantification
> - Understand basic Differential Gene Expression Analysis

<br>

---
---

## Topic overview

This workshop introduces participants to the **fundamentals of RNA-seq data analysis** by building a streamlined pipeline for transcriptomics. Using *E. coli* as a model organism, attendees will work hands-on with real RNA-seq data, **learning to align reads** to a reference genome, **assemble transcripts**, and **perform basic differential gene expression analysis**. By the end of the session, participants will have a working knowledge of key bioinformatics tools— HISAT2, SAMtools, StringTie, and Ballgown—and how these tools fit together in a transcriptomics workflow.

> [!IMPORTANT]
> This workshop uses data obtained by the [National Center of Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/).  
> We are going to be using *E. coli* as the example organism throughout this workshop. All files are made available on CyVerse [here](https://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/datalab/biosciences) or if you are using the command line here: `/iplant/home/shared/cyverse_training/datalab/biosciences`.
>
> - The *E. coli* genome used will be the [ASM584v2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/) (4.6MB).
> - Curated RefSeq annotations (gtf) for the genome can be found in the same location as above.
> - Raw reads used are from [here](https://www.ncbi.nlm.nih.gov/sra/SRX26020219[accn]) (302.8MB) and [here](https://www.ncbi.nlm.nih.gov/sra/?term=SRR30597479) (Note: these sequences come from the [same project](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1158507) and have been processed with an Illumina NovaSeq 6000 machine).
>   - These reads have been fetched with the `sra-tools` (`prefetch <SRA sample>`) and extracted (`fastq-dump <SRA sample>.sra`).

> [!NOTE]
> **Overview of Pipeline**
> As this can be classified as an analysis pipeline, it is important for us to highlight what the major steps are:
> 
> [Raw RNA-seq Data (SRA)]                   
>         ↓                                 
> [Convert to FASTQ]                        
>         ↓
> [Align with Hisat2] 
>         ↓
> [Convert & Sort BAM with SAMtools]
>         ↓
> [Assemble & Quantify with StringTie] 
>         ↓
> [Analyze with Ballgown]

---
---

## Alignment and Data Preparation

<p align="center">
    <img src="https://www.petagene.com/wp-content/uploads/2021/01/HISAT2-performance-with-PetaGene-e1611335725386.png" width="450">
</p>

Image source: [Petagene](https://www.petagene.com/a-practical-example-of-hisat2-smaller-files-same-tools-faster-analysis/). Note, the image above uses the human genome as reference ([GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)) and specific instant specs. We will be using *E. coli* and different resources instead.

In this section, we will learn how to process and prepare RNA-seq data for analysis. FASTQ formatted data has been made available on CyVerse (users should first move this data to a local directory). We will then cover how to align these reads to the *E. coli* reference genome using HISAT2, a tool designed for efficient and accurate RNA-seq alignment. Following alignment, we'll use SAMtools to convert the alignment results from SAM to BAM format and sort them, preparing the data for downstream analysis.

- [Alignment via HISAT2](#alignment-via-hisat2) 
- [Data Formatting Using SAMtools ](#data-formatting-using-samtools)

### Alignment via HISAT2

---

### Data Formatting Using SAMtools 

---
---

## Transcript Assembly and Quantification

<p align="center">
    <img src="https://ccb.jhu.edu/software/stringtie/DE_pipeline_refonly.png" width="450">
</p>

Image source: [StringTie manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual).

This section focuses on assembling and quantifying transcripts from the aligned RNA-seq data. Learners will use StringTie to assemble the RNA-seq reads into transcript structures and quantify gene expression levels. We’ll discuss how to use the GTF annotation file to guide the assembly process and generate a GTF file representing the assembled transcripts. This step provides crucial insights into gene expression and transcript levels.

- [StringTie: Assembling RNA-seq Reads Into Transcripts](#stringtie-assembling-rna-seq-reads-into-transcripts)

> [!NOTE]
> **File formats: GTF vs GFF**

### StringTie: Assembling RNA-seq Reads Into Transcripts

[**StringTie**](https://ccb.jhu.edu/software/stringtie/) is a tool for transcript assembly and quantification from RNA-seq data. It reconstructs transcript structures and estimates their expression levels.

**Key Features**:
- Transcript Assembly: Assembles transcript models from aligned RNA-seq reads.
- Expression Quantification: Estimates transcript abundance in terms of FPKM (Fragments Per Kilobase of transcript per Million mapped reads) or TPM (Transcripts Per Million).
- Output: Provides GTF/GFF files with transcript annotations and expression levels.

---
---

## Exploring Differential Gene Expression

<p align="center">
    <img src="https://raw.githubusercontent.com/alyssafrazee/ballgown/master/figure/plotMeans-1.png" width="450">
</p>

Image source: [Ballgown GitHub repository](https://github.com/alyssafrazee/ballgown). Pictured is a comparison of mean abundances between two groups.

In the final section, participants will be introduced to the basics of differential gene expression analysis using Ballgown. We will cover how to load and interpret the output from StringTie in Ballgown, focusing on identifying and exploring differential expression patterns between samples. This section will give participants a foundational understanding of how to analyze and visualize gene expression data to uncover meaningful biological insights.

- [Visualizing Differential Expression with Ballgown](#visualizing-differential-expression-with-ballgown) 

> [!NOTE]
> We will be needing the use of the statistical software **R**, available in the [CyVerse Bioscience App](https://de.cyverse.org/apps/de/cc046834-5907-11ef-bcd7-008cfa5ae621).

### Visualizing Differential Expression with Ballgown

---
---