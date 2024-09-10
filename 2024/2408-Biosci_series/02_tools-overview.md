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
> - 3:00pm-2:05pm: Welcome and introdution to topic
> - 3:05pm-2:10pm: The Command Line in <5 Minutes
> - 3:10pm-2:20pm: Quality Control Tools
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
    <img src="image" width="450">
</p>

- [content 1](#section-1-subsection-1) 
- [content 2](#section-1-subsection-2)

### Section 1 Subsection 1

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