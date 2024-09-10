# Intro to Essential Tools for Every Bioinformatician

<br>
<br>
<p align="center">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/C674115.jpg/662px-C674115.jpg" width="500">
</p>
<br>

---
>[!important]
> :clock1: **Schedule**
> - 3:00pm-3:05pm: Welcome and introdution to topic
> - 3:05pm-3:10pm: The Command Line in <5 Minutes
> - 3:10pm-3:20pm: Quality Control Tools
> - 3:20pm-3:45pm: Sequence Alignment and Assembly Tools
> - 3:45pm-3:55pm: Variant Calling Tools
> - 3:55pm-4:00pm: Closing remarks

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Exposure to general bioinformatic tools
> - Understanding of genome assembly pipelines

<br>

---
---

## Topic overview

This workshop focuses on the various tools every bioinformatician should know about. This is **NOT** intended to be a deep dive into each and every tool, as each software is incredibly complex and requires precious time to truly be understood. However, this workshop tries to expose the most common tools when it comes to the following popular bioinformatic topics: *quality control (QC), sequence alignment, transcript assembly, variant calling.*

***Jack of all trades, master of none, <u> but oftentimes better than a master of one. </u>***

> [!WARNING]
> **:question: "... but my favourite tool/topic is missing!"**
>
> That's completely understandable: bioinformatics/genomics is such a wide field; This list is non-exhaustive and I will add more tools you may be interested in at the end of the workshop's page.
> 
> **:heavy_exclamation_mark: Please do notice that visualization tools and analysis software have their own sessions! Check the [Data Science Institute's DataLab website](https://datascience.arizona.edu/education/uarizona-data-lab) for more information. :heavy_exclamation_mark:**

<br>

> [!IMPORTANT]
> This workshop uses data obtained by the [National Center of Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/).  
> We are going to be using *E. coli* as the example organism throughout this workshop. All files are made available on CyVerse [here](https://datacommons.cyverse.org/browse/iplant/home/shared/cyverse_training/datalab/biosciences) or if you are using the command line here: `iplant/home/shared/cyverse_training/datalab/biosciences`.
>
> - The *E. coli* genome used will be the [ASM584v2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/) (4.6MB).
> - Curated RefSeq annotations (gff) for the genome can be found in the same location as above.
> - Raw reads used are from [here](https://www.ncbi.nlm.nih.gov/sra/SRX26020219[accn]) (302.8MB) and [here](https://www.ncbi.nlm.nih.gov/sra/?term=SRR30597479) (Note: these sequences come from the [same project](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1158507) and have been processed with an Illumina NovaSeq 6000 machine).
>   - These reads have been fetched with the `sra-tools` (`prefetch <SRA sample>`) and extracted (`fastq-dump <SRA sample>.sra`).
>
> We are also going to be using some famous genes that access:
> - [HBB (hemoglobin subunit beta)](https://www.ncbi.nlm.nih.gov/gene/3043)
> - [CDC28 (yeast cell cycle regulator)](https://www.ncbi.nlm.nih.gov/gene/852457)
> - [Pr55(Gag) (HIV-1 virus encoding gene)](https://www.ncbi.nlm.nih.gov/gene/155030)
> - [rbCl (Rubisco gene)](https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=845212)
> - [lacZ (lactose metabolizing gene)](https://www.ncbi.nlm.nih.gov/gene/945006)

---
---
## The Command Line in <5 Minutes

Most tools used in bioinformatics/genomics are command line tools (not taking into consideration Python or R tools). 

The Command Line is how you communicate with your computer: it sends specific commands that tell software to *execute* specific actions. Following is a quick example of how commands work:

```
    hostname/machine           command        input file      argument/value----------------┐  
            ↓                     ↓                ↓            ↓          ↓                ↓
jovyan@aaf56ff92:~/data-store$ blastn -query <input.fasta> -db nt -out output.txt -evalue 1e-5
  ↑                    ↑                ↑                   ↑       ↑               ↑
user           current location     flag/option            flags/options------------┘

\---------------------------/  \--------------------------------------------------------------/
        shell prompt                                command
```

> [!NOTE]
> 
> - The **shell prompt** (or simply, **prompt**) is always followed by the command.
> - The command can have a number of flag/options.
>    - Usually these can be listed by doing `<command> --help` (*"help me by showing me your options!"*)
>    - Each flag/option requires one of the following:
>        - The path to an **input file**;
>        - An argument (e.g., preferred name of output) or a value (e.g., `1e-5`)
>        - Some options do not require any added argument/value (e.g., `--help` does not require any further input)

> [!TIP]
> Always remember: **Spaces in Paths are *bad*.**
> 

## Quality Control Tools

<p align="center">
    <img src="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/per_base_quality.png" width="450">
</p>

This section covers popular methods to assess quality of data. **FastQC** assesses raw data from sequencing runs, useful before an assembly or analysis starts. **MultiQC** takes outputs of other tools (but including FastQC) to give a more in depth report across multiple samples. 

### FastQC

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

#### Reading the FastQC Results

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

#### MultiQC

MultiQC is similar to FastQC in a sense that both generate reports. However, MultiQC is able to aggregate reports from various sources and for various reasons including but not limited to:

- Sequence Quality Control (from FastQC)
- Alignment and Mapping (from STAR, HISAT2, Bowtiw2, BWA, Minimap2)
- Variant Calling (from GATK, Bcftools)
- Quantifiation (from Salmon, Kallisto)
- Differential Expression (DESeq2)

Notice how MultiQC relies on these tools to create reports first. As an aggregator, the power of MultQC lies in the ability to take all of your produced reports and have them all in a single place.

To run MultiQC, one must first designate a folder with all of the generated reports and analyses. Once the folder is created, once can generate the MultiQC report by doing 

```
multiqc /path/to/folder
```

A report is going to be generated, similar to one available at in the [data commons](https://de.cyverse.org/api/download?path=%2Fiplant%2Fhome%2Fshared%2Fcyverse_training%2Fdatalab%2Fbiosciences%2Ffastqc_multiqc_reports%2Fmultiqc_report.html&attachment=0&url=display-download) or `data/iplant/home/shared/cyverse_training/datalab/biosciences/fastqc_multiqc_reports/` if you are following the workshop. Once logged in CyVerse, [you should be able to view the report from your browser](https://de.cyverse.org/api/download?path=%2Fiplant%2Fhome%2Fshared%2Fcyverse_training%2Fdatalab%2Fbiosciences%2Ffastqc_multiqc_reports%2Fmultiqc_report.html&attachment=0&url=display-download).

This MultiQC report was created by first creating the FastQC reports for samples SRR30597479 and SRR30597506. No other tool was aggregated to the report.

---
---

## Sequence Alignment Tools

<p align="center">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/5/54/Multiple_Sequence_Alignment_Using_ClustalW.jpg/584px-Multiple_Sequence_Alignment_Using_ClustalW.jpg" width="450">
</p>

Image above was generated by the alignment software [ClustalO](https://en.wikipedia.org/wiki/Clustal), a powerful multiple sequence aligner (MSA) that allows for manual manipulation of each base.

Alignment tools are key in the realm of bioinformatics as these tools help with indefying regions of intetrest, similarites between regions, assemblies, and more for DNA, RNA or protein sequences. [There's an acual plethora of aligners you have access to.](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software). 

In this section, we are going to cover [BLAST](https://en.wikipedia.org/wiki/BLAST_(biotechnology)) (Basic Local Alignment Search Tool), [BWA](https://github.com/lh3/bwa) (Burrows-Wheeler Aligner), [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and [HiSAT2](https://daehwankimlab.github.io/hisat2/)

- [BLAST](#the-tool-every-bioinformatician-needs-to-know-blast)
- [BWA](#bwa-mapping-short-reads-to-longer-genomes)
    - [Minimap2](#minimap2-bwa-for-longer-reads)
- [Bowtie2](#bowtie2-fast-and-memory-efficient)
- [SAMTools](#samtools-manipulating-sequences)
- [STAR](#star-aligning-rna-seq-data)


### The Tool Every Bioinformatician Needs to Know: BLAST

<p align="center">
    <img src="https://chanzuckerberg.zendesk.com/hc/article_attachments/20605169683988" width="450">
</p>


One of the most well known tools in the realm of bioinformatics is BLAST. It's one of the tools that is introduced very early on when learning bioinformatics, and is a staple for quick alignment searches. [NCBI-BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is available and launchable from any browser and uses the NCBI compute in order to carry out analyses.

BLAST is able to search for both nucleotide sequences and protein sequences in a quick manner, making it an adequate tool for introducing bioinformatics to a newer audience.

On the command line, BLAST is an even more powerful tool, however it comes with a few too many options. Executing BLAST first requires a genome that has been made into a database with `makeblastdb`.

The full command used is `makeblastdb -in ASM584v2.fna -dbtype nucl -out e_coli_db`. 

For this BLAST exercise, first copy the genome, annotations database and genes to BLAST onto your current location:

```
cp data/iplant/home/shared/cyverse_training/datalab/biosciences/e_coli_gen_dir/ data/iplant/home/shared/cyverse_training/datalab/biosciences/genes .
```

Then, you're ready to BLAST! The command to use is:

```
blastn -query <gene>.fa -db e_coli_gen_dir/e_coli_db/e_coli_db -out results.out
```

> [!IMPORTANT]
> **What's going on here?** Breaking down the command
> - `blastn`: calls the nucleotide version of BLAST
> - `-query <gene.fa>` : here you point to the FASTA file to query the search
> - `-db e_coli_gen_dir/e_coli_db/e_coli_db`: this points to where the *E. coli* database lives. **Note**: inside the `e_coli_gen_dir/e_coli_db` lives a number of files which end with a number of extensions; All together these files make up the database, callable with `e_coli_db`.
> - `-out results.out`: creates an output file called `results.out`

To view the output, one can do:

```
cat results.out
```

If, for example, you blasted any other gene other than the lacZ gene, the `results.out` output will be:

```
BLASTN 2.16.0+


Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb
Miller (2000), "A greedy algorithm for aligning DNA sequences", J
Comput Biol 2000; 7(1-2):203-14.



Database: ASM584v2.fna
           4,315 sequences; 4,025,874 total letters



Query= random_gene

Length=27


***** No hits found *****



Lambda      K        H
    1.33    0.621     1.12 

Gapped
Lambda      K        H
    1.28    0.460    0.850 

Effective search space used: 39525190


  Database: ASM584v2.fna
    Posted date:  Sep 10, 2024  1:42 AM
  Number of letters in database: 4,025,874
  Number of sequences in database:  4,315



Matrix: blastn matrix 1 -2
Gap Penalties: Existence: 0, Extension: 2.5
```

However, running the lacZ gene, will result with:

```
BLASTN 2.16.0+


Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb
Miller (2000), "A greedy algorithm for aligning DNA sequences", J
Comput Biol 2000; 7(1-2):203-14.



Database: ASM584v2.fna
           4,315 sequences; 4,025,874 total letters



Query= NC_000913.3:c366305-363231 Escherichia coli str. K-12 substr.
MG1655, complete genome

Length=3075
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

lcl|NC_000913.3_cds_NP_414878.1_338 [gene=lacZ] [locus_tag=b0344]...  5679    0.0  


>lcl|NC_000913.3_cds_NP_414878.1_338 [gene=lacZ] [locus_tag=b0344] 
[db_xref=UniProtKB/Swiss-Prot:P00722] [protein=beta-galactosidase] 
[protein_id=NP_414878.1] [location=complement(363231..366305)] 
[gbkey=CDS]
Length=3075

 Score = 5679 bits (3075),  Expect = 0.0
 Identities = 3075/3075 (100%), Gaps = 0/3075 (0%)
 Strand=Plus/Plus

Query  1     ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCT  60
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  1     ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCT  60

Query  61    GGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGC  120
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  61    GGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGC  120

Query  121   GAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC  180
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  121   GAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC  180

Query  181   TTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCT  240
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  181   TTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCT  240

Query  241   GAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATC  300
             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  241   GAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATC  300

...............................................................................
....TRUNCATED TO SAVE SPACE....................................................
...............................................................................

Lambda      K        H
    1.33    0.621     1.12 

Gapped
Lambda      K        H
    1.28    0.460    0.850 

Effective search space used: 11966980014


  Database: ASM584v2.fna
    Posted date:  Sep 10, 2024  1:42 AM
  Number of letters in database: 4,025,874
  Number of sequences in database:  4,315



Matrix: blastn matrix 1 -2
Gap Penalties: Existence: 0, Extension: 2.5

```

Resulting in a successful BLAST hit!

There are more ways to customize your BLAST search. Other options include:
- `-outfmt`: which allows you to pick the output format. Ther are 13 different output formats, tabular (`-outfmt 6`) is the most popular.
- `-evalue`: allows you to pick the E-value threshold (expected value), the number of hits you would expect to see by chance. Lower the E-value the more significant the match. A standard E-value is `-evalue 1e-5`
- `-max_target_seqs`: limits the number of alignments to report per query sequence. E.g., `-max_target_seqs 10` will report a max of 10 matches per sequence.
- `num_threads`: specifies the number of threads for parallel processing. A trick: a thread per computer core.

---

> [!NOTE]
> As you may have noticed, BLAST is used to search for alignments between sequences, giving you a number of information per alignment. The following alignment tools are used for a different type of alignment: mapping reads and sequences in order to assemble genomes and RNA-seq data.

> [!Important]
> Due to time constraints, we won't be able to play around with these. From here on, the workshop is going to be more theorical. We will be able to use some of these tools in the coming workshop dates.
---

### BWA: Mapping Short Reads to Longer Genomes

**BWA** (Burrows-Wheeler Aligner) is a widely used software package for mapping sequencing reads against a reference genome. It is designed for high-throughput sequence alignment and is known for its speed and accuracy. BWA can handle both short and long reads and is commonly used in genomics for aligning reads from DNA sequencing experiments.

It has 3 algorithms:
- **BWA-MEM**: The most commonly used algorithm, especially for long reads and high-throughput sequencing. It is designed for accurate and fast alignment of reads.
- **WA-ALN**: Used for aligning short reads, particularly useful for older versions of sequencing technologies.
- **BWA-SW**: Used for aligning longer reads and for more sensitive alignment, but it is slower compared to BWA-MEM.

BWA use cases include: 
- **Variant Calling**: Aligning sequencing reads to detect genetic variations such as SNPs and indels.
- **Gene Expression Analysis**: Aligning RNA-seq reads to measure gene expression levels.
- **Structural Variation Detection**: Identifying larger genomic rearrangements and structural variants.

The implementation of BWA in a pipeline is similar to the following:

1. Index the Reference Genone

```
bwa index reference_genome.fasta
```

2. Aligning Reads with BWA-MEM and redirect the output to a SAM file. 

```
bwa mem reference_genome.fasta reads.fq > aligned_reads.sam

```
**Note**: `reads.fq` in the example above is a reads in a FASTQ format, a popular format that not only contains the sequence information, but also the quality score (similar to the Phred score discussed earlier -- the higher the better.) Quality scores are represented by a mix of symbols and letters ([reference](https://help.basespace.illumina.com/files-used-by-basespace/quality-scores)).

The output SAM file is a format used to store alignment data. SAM is human readable, whilst BAM is the computer readable counterpart (both can be outputted).


> [!IMPORTANT]
> #### Minimap2: BWA for Longer Reads
> 
> As BWA is used mainly for shorter reads, a sister tool has been created in order to align longer reads (e.g., nanopore and PacBio sequences): Minimap2.
> The main difference is in the algorithms used, as Minimap2's algorithm is more efficient and accurate with longer read sequences.
> The two operate with a similar synthax.

---

### Bowtie2: Fast and (Memory) Efficient

---

### SAMtools: Manipulating Sequences

---

### STAR: Aligning RNA-Seq data

---

### Section 1 Subsection 2

---
---

## Transcript Assembly Tools

<p align="center">
    <img src="image" width="450">
</p>

- [content 1](#section-1-subsection-1) 
- [content 2](#section-1-subsection-2)

### Section 1 Subsection 1

---

### Section 1 Subsection 2

---
---

## Variant Calling Tools

<p align="center">
    <img src="image" width="450">
</p>

- [content 1](#section-1-subsection-1) 
- [content 2](#section-1-subsection-2)

### Section 1 Subsection 1

---

### Section 1 Subsection 2

---
---