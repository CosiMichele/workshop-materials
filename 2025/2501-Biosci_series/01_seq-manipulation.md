# An Introduction to Sequence manipulation, Alignment, and Assessment

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
> - 2:00pm-2:15pm: Welcome and introdution to topic
> - 2:15pm-2:25pm: Quality Control for Raw Reads
> - 2:05pm-2:25pm: Alignment and Data Preparation
> - 2:25pm-2:45pm: Transcript Assembly and Quantification
> - 2:45pm-3:00pm: Exploring Differential Gene Expression and closing remarks

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

This workshop introduces participants to the **fundamentals of Sequence manipulation, Alignment, and Assessment** by building a streamlined pipeline for transcriptomics. Using *M. Musculus* (house mouse) as a model organism, attendees will work hands-on with real RNA-seq data, **learning to align reads** to a reference genome, **quantify gene expression**. By the end of the session, participants will have a working knowledge of key bioinformatics tools— HISAT2, SAMtools, featureCounts, and DESeq2—and understand how these tools fit together in a typical RNA-seq workflow.

> [!IMPORTANT]
> This workshop uses data obtained by the [National Center of Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/).  
> We are going to be using *M. Musculus* as the example organism throughout this workshop. All files are made available on CyVerse [here](https://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/datalab/biosciences/) or if you are using the command line here: `/iplant/home/shared/cyverse_training/datalab/biosciences`.
>
> - The *M. Musculus* genome used will be the [GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/).
> - Curated RefSeq annotations (gff) for the genome can be found in the same location as above.
> - RNA-seq data was obtained from the  "RNA-seq and Differential Expression" Workshop by Texas A&M HPC group ([ACES](https://hprc.tamu.edu/aces/)) for the University of Arizona, first held in April 2024.

> [!NOTE]
> **Overview of Pipeline**
> As this can be classified as an analysis pipeline, it is important for us to highlight what the major steps are:
> ```
> [Raw RNA-seq Data QC]               
>         ↓
> [Trimming]               
>         ↓
> [HISAT2 indexing and alignment] 
>         ↓
> [Convert & Sort BAM with SAMtools]
>         ↓
> [Counting Reads in Features with HTSeq] 
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
    - The compiled reads in sam format (`output.sam`)
    - The Jupyter Notebook you can use to replicate the workshop.
    - The R script we will use for the last part of the workshop.

You are now ready for the workshop!

---
---

## Handling Raw Reads with FastQC and TrimGalore 

<p align="center">
    <img src="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/per_base_quality.png" width="450">
</p>

This section covers popular methods to assess quality of data. **FastQC** assesses raw data from sequencing runs, useful before an assembly or analysis starts.

Upon obtaining raw sequencing data and before starting your downstream analysis, it is common that at the beginning of your pipeline you run **FastQC**.

FastQC will help you identify potential issues such as 
- Low quality scores at certain positions.
- Adapter contamination.
- Overrepresented sequences.
- GC content anomalies.

The typical workflow including FastQC is the following: 
- Run FastQC on the raw reads.
- Review the report for quality issues.
- Perform trimming/cleaning if needed (using tools like Trimmomatic or Cutadapt).
- Optionally, re-run FastQC to confirm the data is clean.

To run fastqc, the commands are

```
fastqc <sample_data>.fastq -o /path/to/output
```

Or, if you have multiple samples

```
fastqc *.fastq -o /path/to/output
```

### Reading the FastQC Results

Let's use the figure above to help us understand the output FastQC, called the the **Per base sequence quality** plot, and it provides the distribution of quality scores at each position in the read across all reads.
- The y-axis gives the quality scores, while the x-axis represents the position in the read. 
- The color coding of the plot denotes what are considered high, medium and low quality scores.
- The yellow box represents the 25th and 75th percentiles, with the red line as the median. The whiskers are the 10th and 90th percentiles.
- The blue line represents the average quality score for the nucleotide.

> [!IMPORTANT]
> **What's the Phred Quality Score?**
> 
> The y-axis typically represents [Phred](https://en.wikipedia.org/wiki/Phred_(software)) quality scores, which measure the probability of an incorrect base call. Quality scores are logarithmic and can theoretically go beyond 40. Scores above 40 represent very high confidence in base calls, with probabilities of incorrect base calls being very low (less than 1 in 10,000). Setting the y-axis limit at 40 helps in standardizing the plot and avoiding unnecessary scaling for the majority of sequencing data, where scores often peak around 30 to 40. Scores higher than 40 are relatively rare in most datasets.
><p align="center">
>    <img src="https://upload.wikimedia.org/wikipedia/commons/e/e7/Phred_Figure_1.jpg" width="450">
></p>

Other important results:
- **Per tile sequence quality** (only if sequencing was carried out with Illumina): The plot shows the deviation from the average quality for each tile, blue = good, red = bad (heat scale blue to red). [Here's an example of a not so great result](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/per_tile_quality.png) ([source](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/12%20Per%20Tile%20Sequence%20Quality.html)).
- **Per sequence quality scores**: the average quality score on the x-axis and the number of sequences with that average on the y-axis.
- **Per base sequence content**:  plots out the proportion of each base position in a file for which each of the four normal DNA bases has been called. [In a random library you would expect that there would be little to no difference between the different bases of a sequence run, so the lines in this plot should run parallel with each other](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html). **Note:** this will always fail with RNA-seq data.
- **Per sequence GC content**: plots the GC distribution over all sequences
- **Per base N content**: If a sequencer is unable to make a base call with sufficient confidence then it will normally substitute an N rather than a conventional base call. This module plots out the percentage of base calls at each position for which an N was called.
- **Sequence Length Distribution**: shows the distribution of sequence length.

### Raw Quality Control

Execute the command: 

```
 fastqc -t 2 -o . Control1_R1.fastq.gz Control1_R2.fastq.gz
```

The output from FastQC is:

<p align="center">
    <img src="https://raw.githubusercontent.com/CosiMichele/workshop-materials/refs/heads/main/2024/2408-Biosci_series/assets/control_r1.png" width="450">
</p>

and

<p align="center">
    <img src="https://raw.githubusercontent.com/CosiMichele/workshop-materials/refs/heads/main/2024/2408-Biosci_series/assets/control_r2.png" width="450">
</p>

How would you interpret these readings?

### Trimming

**TrimGalore**, is a tool for trimming low-quality bases and adapter sequences from high-throughput sequencing reads.

Based on the QC output above, we use TrimGalore's default settings which trims any read that is below length 20. We add the `--paired` option/flag to only return paired reads.

```
trim_galore --paired --fastqc Control1_R1.fastq.gz Control1_R2.fastq.gz
```

The `--fastqc` option/flag outputs a report we can read in FastQC.

<p align="center">
    <img src="https://raw.githubusercontent.com/CosiMichele/workshop-materials/refs/heads/main/2024/2408-Biosci_series/assets/control_r1_val.png" width="450">
</p>

and

<p align="center">
    <img src="https://raw.githubusercontent.com/CosiMichele/workshop-materials/refs/heads/main/2024/2408-Biosci_series/assets/control_r2_val.png" width="450">
</p>

Do you think TrimGalore did well? (yes).

---
---

## Alignment and Data Formatting

<p align="center">
    <img src="https://www.petagene.com/wp-content/uploads/2021/01/HISAT2-performance-with-PetaGene-e1611335725386.png" width="450">
</p>

Image source: [Petagene](https://www.petagene.com/a-practical-example-of-hisat2-smaller-files-same-tools-faster-analysis/). Note: the image above uses the human genome as reference ([GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)) and specific instant specs.

In this section, we will learn how to process and prepare RNA-seq data for analysis. FASTQ formatted data has been made available on CyVerse (users should first move this data to a local directory). We will then cover how to align these reads to the *M. Musculus* reference genome using HISAT2, a tool designed for efficient and accurate RNA-seq alignment. Following alignment, we'll use SAMtools to convert the alignment results from SAM to BAM format and sort them, preparing the data for downstream analysis.

- [Alignment via HISAT2](#alignment-via-hisat2) 
- [Data Formatting Using SAMtools ](#data-formatting-using-samtools)

---

### HISAT2: Aligning Spliced RNA-Seq data

<p align="center">
    <img src="https://daehwankimlab.github.io/hisat2/assets/img/ogp.png" width="250">
</p>

[**HISAT2**](https://daehwankimlab.github.io/hisat2/) is a tool optimized to deal with spliced junctions. It uses the similar algorithms to the [BWA aligner](https://bio-bwa.sourceforge.net/bwa.shtml), paired with [FM-index](https://en.wikipedia.org/wiki/FM-index), a compression algorithm for improved memory efficiency. HISAT2 can handle larger genomes and more complex transcriptomes.

HISAT2 requires an indexed genome, obtainable with a simple synthax:

```
hisat2-build reference_genome.fa genome_index
```

Basic synthax:

```
hisat2 -x genome_index -1 read1.fastq -2 read2.fastq -S output.sam
```

The options used here are:

- `-x`: Specifies the index of the reference genome.
- `-1`: Specifies the file containing the first reads of paired-end data.
- `-2`: Specifies the file containing the second reads of paired-end data.
- `-S`: Specifies the output file in SAM format.

Some other popular options are:

- `-U <single_end_read.fq>`: for single end reads
- `--dta`: Use this option for transcriptome alignments to improve alignments for transcript quantification.
- `--rna-strandness`: Specify RNA strandness (use `FR` for forward, `RF` for reverse, or `NONE` for no specific strandness).
- `--spliced-align`: Ensure the aligner considers splicing events.
- `--max-intron-len <preferred intron length>`: Set the maximum allowed intron length for spliced alignments.

---

### Alignment via HISAT2

Prior to alignment with HISAT2, we require to index the genome. This can be carried out with the following command synthax. Note, we have added the `-p <threads>` flag for threads (decreasing the indexing time for larger genomes).

HISAT2 is used to index the model organism of the Mouse first.

```
hisat2-build -p <threads> <reference-genome.fasta> <reference_index>
```

Thus, do:

```
hisat2-build -p 16 GCF_000001635.27_GRCm39_genomic.fna GCF_000001635.27_GRCm39_genomic
```

This command will build the *M. musculus* database using 16 threads (time estimation: ~10m).

> [!NOTE]
> If you ever need to index larger genomes, do also include the `--large-idex` flag.
> Other flags that are worth using:
> - `--ss <splice-sites.txt>`: where you point HISAT2 to a txt file where known splice sites are reported.
> - `--exon <exons.txt>:` similar to the splice site option, but for exon information. 

Now, you can align your reads. The baseline HISAT2 command for alignment is the following:

```
hisat2 -x <reference_index> -1 <sample_R1.fastq> -2 <sample_R2.fastq> -S <output.sam>
```

Then, HISAT2 is used to align our trimmed reads about outputting a SAM file.

```
 hisat2 -x GCF_000001635.27_GRCm39_genomic -p 2 -1 Control1_R1_val_1.fq.gz -2 Control1_R2_val_2.fq.gz -S Control1.sam
        .
        .
        .
236499 reads; of these:
 236499 (100.00%) were paired; of these:
30736 (13.00%) aligned concordantly 0 times
197200 (83.38%) aligned concordantly exactly 1 time
8563 (3.62%) aligned concordantly >1 times
----
30736 pairs aligned concordantly 0 times; of these:
 3583 (11.66%) aligned discordantly 1 time
----
27153 pairs aligned 0 times concordantly or discordantly; of these:
 54306 mates make up the pairs; of these:
 30660 (56.46%) aligned 0 times
 21188 (39.02%) aligned exactly 1 time
 2458 (4.53%) aligned >1 times
93.52% overall alignment rate
```

With 93.52% overall alignmet, we score well!

---

### Data Formatting Using SAMtools 

Notice how the output from our most recent alignment is in SAM format. To operate SAM/BAM files, we use **SAMtools**.

SAMtools is a suite of tools used for manipulating and analyzing files in the SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) formats. These formats are commonly used to store alignment data from various sequencing technologies. SAMtools provides a variety of functionalities to process alignment files, including sorting, indexing, and variant calling.

As mentioned before, SAM and BAM files are essentially the same, with one (SAM, Sequence Alignment/Map, larger size) being human readable whilst the other (BAM, Binary Alignment/Map, more compact) is machine readable. 

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
samtools sort -@ 16 -m 8G output.bam -o output_sorted.bam
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

It is typical to first convert SAM files to BAM as these save space:

```
samtools view -bS input.sam > output.bam
```

SAMtools is able to **sort** SAM/BAM files by the genomic coordinates of the alignments, which is necessary for many downstream analyses:

```
samtools sort output.bam -o sorted.bam
```

Typically this is followed by **indexing**, which creates an index file for BAM files to allow fast access to specific regions of the genome:

```
samtools index sorted.bam
```

**Returning to our pipeline**, we then use SAMtools to sort our alignment and convert SAM into BAM ... 

```
samtools sort -@ 16 -o Control1_sorted.bam Control1.sam
```

... and idex the sorted BAM file.


```
samtools index Control1_sorted.bam
```

> [!Note]
> Notice how we did not use `samtools view -bS input.sam > output.bam`
> This is because certain actions with SAMtools can output directly in the format that we require.

> [!IMPORTANT] 
> **There is more to SAMtools that we are not covering.** 
> Users can then create aligment statistics with
> ```
> samtools flagstat sorted.bam
> ```
> Or create a pileup file for [variant calling](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/#:~:text=What%20is%20variant%20calling%3F,creating%20BAM%20or%20CRAM%20files.): 
> ```
> samtools mpileup -f reference.fasta sorted.bam > variants.bcf
> ```
> A pileup file is a data format used to represent the sequence coverage at each position in the genome across multiple reads. It provides a compact summary of read alignments and is particularly useful for variant calling and other genomic analyses.
> 
> It is useful for the following reasons:
> 
> - **Coverage Information**: shows how many reads cover each position in the genome, helping to identify regions with high or low sequencing coverage.
> - **Can be used for Base Calls**: lists the base calls (nucleotides) observed at each position along with their quality scores.
> - **Variations**: allows for the identification of variations such as SNPs (single nucleotide polymorphisms) and indels (insertions and deletions) by comparing observed bases to the reference genome.

---

## Counting Reads in Features

HTSeq-count is a tool for counting the number of reads mapped to specific genomic features (like genes) in RNA-seq data.

We use flags:
- `-r pos`: tells HTSeq that the reads in the BAM file are sorted by position (i.e., genomic coordinate), which can be faster if the BAM file is already sorted by position.
- `-i gene`: sets the feature ID attribute to use for counting ("count the genes for me please!")

```
htseq-count -r pos -i gene Control1_sorted.bam GCF_000001635.27_GRCm39_genomic.gff > Control1_counts.txt
```

As it is a txt file, one can open the txt and try to sort through the file themselves, however this is not suggested.

We will be working on this file in the following week.

---
---