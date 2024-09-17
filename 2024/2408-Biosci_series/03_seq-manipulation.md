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

This workshop introduces participants to the **fundamentals of RNA-seq data analysis** by building a streamlined pipeline for transcriptomics. Using *E. coli* as a model organism, attendees will work hands-on with real RNA-seq data, **learning to align reads** to a reference genome, **quantify gene expression**, and **perform differential gene expression analysis**. By the end of the session, participants will have a working knowledge of key bioinformatics tools— HISAT2, SAMtools, featureCounts, and DESeq2—and understand how these tools fit together in a typical RNA-seq workflow.

> [!IMPORTANT]
> This workshop uses data obtained by the [National Center of Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/).  
> We are going to be using *E. coli* as the example organism throughout this workshop. All files are made available on CyVerse [here](https://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/datalab/biosciences) or if you are using the command line here: `/iplant/home/shared/cyverse_training/datalab/biosciences`.
>
> - The *E. coli* genome used will be the [ASM75055v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000750555.1/) (4.6MB).
> - Curated RefSeq annotations (gtf) for the genome can be found in the same location as above.
> - Raw reads are taken from a recent study, [Transcriptome RNA Sequencing Data Set of Differential Gene Expression in *Escherichia coli* BW25113 Wild-Type and *slyA* Mutant Strains](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8142576/) by Frank J. Stewart (Microbiology Resoururce Announcements, 2021). 
>   - Reads used are [SRX10348166](https://www.ncbi.nlm.nih.gov/sra/SRX10348166) (run SRR13970433), [SRX10348167](https://www.ncbi.nlm.nih.gov/sra/SRX10348167) (run SRR13970434), [SRX10348168](https://www.ncbi.nlm.nih.gov/sra/SRX10348168) (run SRR13970435), and [SRX10348169](https://www.ncbi.nlm.nih.gov/sra/SRX10348169) (run SRR13970436). Reads have been processed using the Illumina NextSeq 500 machine.
>   - These reads have been fetched with the `sra-tools` (`prefetch <SRA sample>`) and extracted (`fastq-dump <SRA sample>.sra`).
>   - Note: 

> [!NOTE]
> **Overview of Pipeline**
> As this can be classified as an analysis pipeline, it is important for us to highlight what the major steps are:
> ```
> [Raw RNA-seq Data (SRA)]                  # 
>         ↓                                 # Pre-made data in preparation for the workshop.
> [Convert to FASTQ]                        #
>         ↓
> [Align with HISAT2] 
>         ↓
> [Convert & Sort BAM with SAMtools]
>         ↓
> [Quantify Reads with featureCounts] 
>         ↓
> [Perform Differential Expression Analysis with DESeq2]
> ```

---
---

## Getting Things Started

In order to get things started:

1. Execute the [JupyterLab Bioscience CyVerse App](https://de.cyverse.org/apps/de/cc046834-5907-11ef-bcd7-008cfa5ae621) and open the Terminal with 16 CPU cores (min and max) and 16GB RAM memory.
2. Initiate GoCommands, the CyVerse in-house transfer tool using `gocmd init`.
    - Press enter (don't input anything) when propted for **iRODS Host**, **Port**, **Zone**; put your CyVerse username when asked for **iRODS Username** and your CyVerse password when asked for **iRODS Password** (**Note**: you will not see the password being typed as per standard Ubuntu security).
3. Copy the genome files using the following command:

    ```
    gocmd get --progress /iplant/home/shared/cyverse_training/datalab/biosciences/e_coli_genomic/
    ```

    This will copy:
    
    - The *E. coli* reference genome (`ASM75055v1_genomic.fna`)
    - The reference transcript annotations (`ASM75055v1_genomic.gtf`)
    - The processed fastq files:
        - SRR13970433.fastq  
        - SRR13970434.fastq  
        - SRR13970435.fastq  
        - SRR13970436.fastq

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
hisat2-build -p 16 ASM75055v1_genomic.fna ecoli_index
```

This command will build the *E. coli* database using 16 threads (time estimation: ~<10s).

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
hisat2 -p 16 -x ecoli_index -U SRR13970433.fastq,SRR13970434.fastq,SRR13970435.fastq,SRR13970436.fastq -S output.sam
```

The command will run for ~1-2 minutes, creating the `output.sam` file (~3.8GB, time estimate ~10-15m)

You will notice that this will wild a really poor alignmet:

```
11395529 reads; of these:
  11395529 (100.00%) were unpaired; of these:
    11232207 (98.57%) aligned 0 times
    148028 (1.30%) aligned exactly 1 time
    15294 (0.13%) aligned >1 times
1.43% overall alignment rate
```

**Only 1.30% aligned!!!** ... so what now?

> [!NOTE]
> **Exercise: Downlaod more read sequences (not required but recommended)**
> Base your searches on the SRA files presented in the [Transcriptome RNA Sequencing Data Set of Differential Gene Expression in *Escherichia coli* BW25113 Wild-Type and *slyA* Mutant Strains](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8142576/) publication.
> 1. Go to the NCBI SRA repository: https://www.ncbi.nlm.nih.gov/sra/
> 2. In the search field, search for *E. coli*. Select any of the experiments, *or* if you know any experiment you know of, search for that! \*wink\* \*wink\* 
> 3. In the terminal, download the reads:
>    - This uses the `sra-toolkit` which is a powerful tool that allows you to download any experimental data hosted by NCBI.
>    - Use `prefetch <sra run #>`, such as `prefetch SRR30660439`. You can get this number once you open one of the experiments (close to the bottom of the page).
> 4. extract the reads using `fastq-dump <sra run #>` such as `fastq-dump SRR30660439`. This will generate the fastq file you can use for alignment.
> You will notice that the alignment number will increase (perhaps, not by much as these are from different experiments).
> 
> Question: Why do you think the alignment is so low?

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
samtools view -@ 16 -S -b output.sam > output.bam
```

The command `samtools view` converts SAM to BAM. The `-S` flag specifies the input SAM, `-@ <threads>` to speed up processing, and `-b` outputs the BAM file. Notice how the command generates a much smaller file (~3.8GB vs ~700MB).

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
samtools index -@ 16 -b output_sorted.bam
```

This creates a `.bai` file, which is your BAM index file.

> [!NOTE]
> Check Stats (Optional)
> You can check the stats of your indexed BAM file by doing `samtools flasgstat -@ 16 output_sorted.bam

---
---

## Quantifying Gene Expression

<p align="center">
    <img src="https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/union.png" width="350">
</p>

Image source: [Harvard Chan Bioinformatics Core](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html).

In this step, learners will learn how to use featureCounts to assign aligned reads from RNA-seq data to genomic features, such as genes. This tool efficiently counts the number of reads mapping to each gene, providing a crucial input for downstream differential expression analysis. By accurately quantifying gene expression, participants will better understand how different genes are expressed in their RNA-seq dataset.

- [Quantify Reads with featureCounts](#quantify-reads-with-featurecounts)

> [!NOTE]
> **File formats: GFF vs GTF**
> Common feature annotations are found in 2 formats: GFF and GTF. Both files contain tabular information with 9 columns: chromosome, source, feature, start, end, score, strand, frame, attributes.
> GFF:
> - General Feature Format
> - Describes any genomic feature
> GTF:
> - Gene Transfer Format
> - Contains only gene annotation 

### Quantify Reads with featureCounts

 [featureCounts](https://subread.sourceforge.net/featureCounts.html) is used for counting the number of reads that map to each feature (e.g., genes, exons) in a reference annotation. It operates on aligned RNA-Seq reads (e.g., SAM/BAM files) and is focused on providing read counts for predefined genomic features.

> [!NOTE]
> There are different path one can take in order to proceed the pipeline. Another option, for example, is to use StringTe and Ballgown.
>  StringTie is primarily used for assembling transcripts and estimating their expression levels from RNA-Seq data. It can create a transcriptome assembly from aligned RNA-Seq reads and generate gene and transcript expression levels.
> Meanwhile, featureCounts is used for counting the number of reads that map to each feature (e.g., genes, exons) in a reference annotation. It operates on aligned RNA-Seq reads (e.g., SAM/BAM files) and is focused on providing read counts for predefined genomic features.
>
> **Both serve important purposes, it all depends what the end goal is.**

We run featureCounts with the following command synthax:

```
featureCounts -a <reference.gtf> -o counts.txt <output_sorted.bam>
```

There are a few flags that we want to use:

- `-a <reference.gtf>`: Specifies the annotation file in GTF format.
- `-o <output.txt>`: Specifies the output file for counts.
- `-T <threads>`: Uses 4 threads (adjust according to the number of available cores).
- `-g gene-id`: specifies the gene ID in the output.
- `--ignoreDup`: counts only non-duplicate reads (useful if your data might include duplicates).

thus we run:

```
featureCounts -T 16 -a ASM75055v1_genomic.gtf -o counts.txt -g gene_id --fraction --ignoreDup output_sorted.bam
```

We can then check the `counts.txt` output by doing `head counts.txt`. Output should look something like this:

```
# Program:featureCounts v2.0.3; Command:"featureCounts" "-T" "16" "-a" "ASM75055v1_genomic.gtf" "-o" "counts.txt" "-g" "gene_id" "--fraction" "--ignoreDup" "output_sorted.bam" 
Geneid  Chr     Start   End     Strand  Length  output_sorted.bam
BW25113_RS01010 NZ_CP009273.1   220258  221799  +       1542    10
BW25113_RS01015 NZ_CP009273.1   221868  221944  +       77      0
BW25113_RS01020 NZ_CP009273.1   221987  222062  +       76      0
BW25113_RS01025 NZ_CP009273.1   222246  225149  +       2904    51
BW25113_RS01030 NZ_CP009273.1   225245  225360  +       116     0
BW25113_RS01035 NZ_CP009273.1   225415  225491  +       77      0
BW25113_RS01085 NZ_CP009273.1   233418  233494  +       77      0
BW25113_RS01235 NZ_CP009273.1   258582  258657  +       76      0
```

This output tells us the following:

- **Geneid**: The identifier for the gene, here using the *E. coli* gene names (e.g., BW25113_RS01010). This tells you which gene each row is summarizing.
- **Chr**: The chromosome (or contig) where the gene is located. For instance, NZ_CP009273.1 refers to the name of the contig (which could be the chromosome or a plasmid).
- **Start**: The starting position of the gene on the chromosome (in base pairs). For example, the gene BW25113_RS01010 starts at position 220258.
- **End**: The ending position of the gene on the chromosome. For BW25113_RS01010, the gene ends at position 221799.
- **Strand**: The strand the gene is located on, either "+" (forward strand) or "−" (reverse strand). In this case, all genes are on the forward strand (+).
- **Length**: The length of the gene (in base pairs). For example, the length of BW25113_RS01010 is 1542 bp.
- *output_sorted*.bam: This column shows the number of reads that were assigned to each gene in your sorted BAM file (output_sorted.bam). For BW25113_RS01010, there are 10 reads assigned, while other genes like BW25113_RS01015 and BW25113_RS01020 have 0 reads assigned.

As we are going to be using DESeq2, we need to modify the formatting and keep only the features that we care about.

>[!IMPORTANT]
> ALWAYS. KEEP. YOUR DATA. Make a backup, rename files differently, but always find a way to keep your data (so that when you make mistakes, you can always go back!).

There are a few ways to modify this data; Using R is a popular option, but we can be quicker with the Command Line. We are going to use a powerful editing software called [awk](https://www.geeksforgeeks.org/awk-command-unixlinux-examples/) to remove the columns we do not need. Run the following:

```
awk 'NR > 1 {print $1, $NF}' counts.txt > cleaned_counts.txt
```

- `NR > 1` skips the first line (metadata).
- `{print $1, $NF}` prints only the first column (gene ID) and the last column (count). If we do `head cleaned_counts.txt` you will see something similar to this:

```
Geneid output_sorted.bam
BW25113_RS01010 10
BW25113_RS01015 0
BW25113_RS01020 0
BW25113_RS01025 51
BW25113_RS01030 0
BW25113_RS01035 0
BW25113_RS01085 0
BW25113_RS01235 0
BW25113_RS23035 0
```

We then need to use an editing sofware (e.g., `nano`) to edit output_sorted.bam to `sample1` (as we only have a single sample).

```
Geneid sample1
BW25113_RS01010 10
BW25113_RS01015 0
BW25113_RS01020 0
BW25113_RS01025 51
BW25113_RS01030 0
BW25113_RS01035 0
BW25113_RS01085 0
BW25113_RS01235 0
BW25113_RS23035 0
```

The file is not ready for the next step!

---

## It's *R* time

As we are going to be using R for the next step, in the JupyterLab click **Laucher :heavy_plus_sign:** and the **RStudio :arrow_upper_right:** vignette.

This will open RStudio in a new tab.

---

## Exploring Differential Gene Expression

<p align="center">
    <img src="https://hbctraining.github.io/DGE_workshop/img/normalization_methods_depth.png" width="450">
</p>

Image source: [Harvard Chan Bioinformatics Core](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html). In the figure above, each pink and green rectangle represents a read aligned to a gene. Reads connected by dashed lines connect a read spanning an intron.

In the following analysis we use R and [DESeq2](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) to perform normalization and some differential gene expression analysis. The analysis itself leads us to uncover which genes are most highly expressed in the sample, normalize the data to make it comparable, and visualize gene expression patterns. 

Code is available in the R script provided at `data/iplant/home/shared/cyverse_training/datalab/biosciences/deseq2_viz.R`; copy it to your location (`cp data/iplant/home/shared/cyverse_training/datalab/biosciences/deseq2_viz.R .`). Here's a general summary of what is covered in this section:

1. **Loading the Count Data**: We begin by reading in the raw counts matrix, which contains the number of reads assigned to each gene across different samples. This count data is necessary for downstream analysis.
2. **Creating the DESeqDataSet Object**: We create a DESeqDataSet object, which structures the raw count data and defines the experimental design (conditions or groups being compared). This object is essential for normalizing and analyzing RNA-seq data.
3. **Normalization of Count Data**: RNA-seq data can have biases due to varying sequencing depths across samples. DESeq2 performs size factor estimation to account for these differences, generating normalized counts that can be compared across genes and samples.
4. **Visualization**: With normalized counts, we generate visualizations like heatmaps using ggplot2 in R. These visualizations provide an overview of gene expression patterns across different conditions or samples.
5. **Identifying Highly Expressed Genes**: By calculating the average expression for each gene, we can rank and identify highly expressed genes. These genes may be of particular interest for further study.

By the end, you should be able to create a small table of the highest expressed genes in our alignment (which should look like the following)!

```
                           Gene Average_Expression
BW25113_RS19530 BW25113_RS19530                378
BW25113_RS16985 BW25113_RS16985                173
BW25113_RS19520 BW25113_RS19520                140
BW25113_RS13525 BW25113_RS13525                132
BW25113_RS17000 BW25113_RS17000                 96
BW25113_RS01025 BW25113_RS01025                 51
BW25113_RS13535 BW25113_RS13535                 34
BW25113_RS20600 BW25113_RS20600                 30
BW25113_RS19995 BW25113_RS19995                 28
BW25113_RS20800 BW25113_RS20800                 12
```

---
---