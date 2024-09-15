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
> ```
> [Raw RNA-seq Data (SRA)]                  # 
>         ↓                                 # These have been pre-made in preparation for the workshop.
> [Convert to FASTQ]                        #
>         ↓
> [Align with Hisat2] 
>         ↓
> [Convert & Sort BAM with SAMtools]
>         ↓
> [Assemble & Quantify with StringTie] 
>         ↓
> [Analyze with Ballgown]
> ```

---
---

## Getting Things Started

In order to get things started:

1. Execute the [JupyterLab Bioscience CyVerse App](https://de.cyverse.org/apps/de/cc046834-5907-11ef-bcd7-008cfa5ae621) and open the Terminal.
2. Initiate GoCommands, the CyVerse in-house transfer tool using `gocmd init`.
    - Press enter (don't input anything) when propted for **iRODS Host**, **Port**, **Zone**; put your CyVerse username when asked for **iRODS Username** and your CyVerse password when asked for **iRODS Password** (**Note**: you will not see the password being typed as per standard Ubuntu security).
3. Copy the genome files using the following command:

    ```
    gocmd get --progress /iplant/home/shared/cyverse_training/datalab/biosciences/e_coli_gen_dir
    ```

    This will copy:
    
    - The *E. coli* reference genome (`ASM584v2.fna `)
    - The reference transcript annotations (`ASM584v2_genomic.gtf`)
4. Copy the SRA files with:
    ```
    gocmd get --progress /iplant/home/shared/cyverse_training/datalab/biosciences/e_coli_sra
    ```
5. Optional (recommended): execute the following command to move files to a single folder
    ```
    mkdir week_3 && \ 
    mv e_coli_gen_dir/ASM584v2.fna e_coli_gen_dir/ASM584v2_genomic.gtf week_3/ && \
    mv e_coli_sra/SRR30597*/*.fastq week_3 && \
    cd week_3
    ```

You are now ready for the workshop!

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

Prior to alignment with HISAT2, we require to index the genome. This can be carried out with the following command synthax. Note, we have added the `-p <threads>` flag for threads (decreasing the indexing time for larger genomes):

```
hisat2-build -p <threads> <reference-genome.fasta> <reference_index>
```

Thus, do:

```
hisat2-build -p 8 ASM584v2.fna ecoli_index
```

This command will build the *E. coli* database using 8 threads (time estimation: ~10s).

> [!NOTE]
> If you ever need to index larger genomes, do also include the `--large-idex` flag.
> Other flags that are worth using:
> - `--ss <splice-sites.txt>`: where you point HISAT2 to a txt file where known splice sites are reported.
> - `--exon <exons.txt>:` similar to the splice site option, but for exon information. 

Now, you can align your reads. The baseline HISAT2 command for alignment is the following:

```
hisat2 -x <reference_index> -1 <sample_R1.fastq> -2 <sample_R2.fastq> -S <output.sam>
```

Since we are using single reads and not pair ends, we should use the following command:

```
hisat2 -x <index_basename> -U <reads.fastq> -S <output.sam>
```

We will also use the `<read1.fastq>,<read2.fastq>` synthax to include both the reads available to us, as well as using the `-p 8` to run 8 threads. Execute:

```
hisat2 -p 8 -x ecoli_index -U SRR30597479.fastq,SRR30597506.fastq -S output.sam
```

The command will run for ~1-2 minutes, creating the `output.sam` file (~5.1GB)

You will notice that this will wild a really poor alignmet:

```
8343493 reads; of these:
  8343493 (100.00%) were unpaired; of these:
    8343461 (100.00%) aligned 0 times
    32 (0.00%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
0.00% overall alignment rate
```

**Only 32 aligned!!!** ... so what now?

> [!NOTE]
> **Exercise: Downlaod more read sequences (not required but recommended)**
> 1. Go to the NCBI SRA repository: https://www.ncbi.nlm.nih.gov/sra/
> 2. In the search field, search for *E. coli*. Select any of the experiments, *or* if you know any experiment you know of, search for that!
> 3. In the terminal, download the reads:
>    - This uses the `sra-toolkit` which is a powerful tool that allows you to download any experimental data hosted by NCBI.
>    - Use `prefetch <sra run #>`, such as `prefetch SRR30660439`. You can get this number once you open one of the experiments (close to the bottom of the page).
> 4. extract the reads using `fastq-dump <sra run #>` such as `fastq-dump SRR30660439`. This will generate the fastq file you can use for alignment.
> You will notice that the alignment number will increase (perhaps, not by much as these are from different experiments).
> 
> Why do you think the alignment is so low?

As this workshop is concentrated around the techniques, we will continue with the poor alignment. *All is good!*

---

### Data Formatting Using SAMtools 

Using SAMtools, we are going to do a few things:

1. Convert SAM to BAM. This allows for the downstream analyses to be quicker.
2. Sort the reads by genomic coordinates.
3. Index the output BAM, allowing fast access to specific regions in the BAM file.
4. Check stats.

To convert SAM to BAM, one must do:

```
samtools view -S -b <output.sam> > <output.bam>
```

The command `samtools view` converts SAM to BAM. The `-S` flag specifies the input SAM, and -b outputs the BAM file. Notice how the command generates a much smaller file (~5.1GB vs ~1.1GB).

The next command is the `samtools sort` command, with `-@ <threads>` to specify how many threads you want to run it with, `-m <memory per thread followed by size (G for GB)>`. It's good practice to find ways to multithread your analyses.

```
samtools sort -@ 8 -m 4G output.bam -o output_sorted.bam
```

Sorting refers to organizing the alignments in a BAM (or SAM) file based on the genomic coordinates where the reads align to the reference genome. It makes a lot of downstream analyses quicker and more accurate.

Before sorting, the BAM file might look like this:

|Read ID | Chromosome |	Position |
| --- | --- | ---|
| R1 |chr1 | 500 |
| R2 | chr2 | 1000 |
| R3 | chr1 | 1500 |

Notice how read ID R2 is on a different chromosome. Post sorting, read ID R3 should be moved closer to R1, as both are on the same chromosome (quicker access downstream).

|Read ID | Chromosome |	Position |
| --- | --- | ---|
| R1 |chr1 | 500 |
| R3 | chr1 | 1500 |
| R2 | chr2 | 1000 |

Now you can index your sorted output file; Use `-@ <threads>` to speed things up and `-b` to point to the file.

```
samtools index -@ 8 -b output_sorted.bam
```

This creates a `.bai` file, which is your BAM index file.

> [!NOTE]
> Check Stats (Optional)
> You can check the stats of your indexed BAM file by doing `samtools flasgstat -@ 8 output_sorted.bam

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