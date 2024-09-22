# Using Nextflow for Streamlining Bioscience Analysis Pipelines

<br>
<br>
<p align="center">
    <img src="https://workflows.community/images/stories/nextflow-story.png" width="650">
</p>
<br>

Image credits: [Workflows community's article by Paolo di Tommaso (NextFlow creator)](https://workflows.community/images/stories/nextflow-story.png).


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
> [Trim reads with Trimmomatic]            
>         ↓
> [Align with HISAT2] 
>         ↓
> [Convert & Sort BAM with SAMtools]
>         ↓
> [Output alignment statistics]
> ```

<br>
<br>
<p align="center">
    <img src="https://i.kym-cdn.com/photos/images/newsfeed/002/546/187/fb1.jpg" width="650">
</p>
<br>

Image credits: [It's Always Sunny in Philadelphia](https://en.wikipedia.org/wiki/It%27s_Always_Sunny_in_Philadelphia).

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
    - The compiled reads in sam format (`output.sam`)
    - The NextFlow code required for the workshop.

You are now ready for the workshop!



<details>
  <summary>Click here for the raw NextFlow code</summary>

```  
#!/usr/bin/env nextflow

params.reads = "./data/*.fastq"  // Path to input FastQ files
params.index = "./hisat2_index/genome"  // Path to HISAT2 index
params.threads = 4  // Number of threads for each process
params.outputDir = "./results"  // Output directory

// Define the input channel for FastQ files
Channel
    .fromPath(params.reads)
    .set { reads_ch }

// Step 1: Run FastQC in parallel on all input reads
process fastqc {
    tag "${sample_id}"

    input:
    path read_file from reads_ch

    output:
    path "${sample_id}_fastqc.zip" into fastqc_ch

    script:
    sample_id = read_file.baseName
    """
    fastqc ${read_file} -o ./
    """
}

// Step 2: Trim adapters and low-quality sequences using Trimmomatic
process trimReads {
    tag "${sample_id}"

    input:
    path read_file from reads_ch

    output:
    path "${sample_id}_trimmed.fastq" into trimmed_ch

    script:
    sample_id = read_file.baseName
    """
    trimmomatic SE -threads ${params.threads} ${read_file} ${sample_id}_trimmed.fastq SLIDINGWINDOW:4:20 MINLEN:36
    """
}

// Step 3: Align trimmed reads using HISAT2 in parallel
process alignReads {
    tag "${sample_id}"

    input:
    path read_file from trimmed_ch

    output:
    path "${sample_id}.sam" into aligned_ch

    script:
    sample_id = read_file.baseName
    """
    hisat2 -p ${params.threads} -x ${params.index} -U ${read_file} -S ${sample_id}.sam
    """
}

// Step 4: Convert SAM to BAM and sort using SAMtools
process sortBam {
    tag "${sample_id}"

    input:
    path sam_file from aligned_ch

    output:
    path "${sample_id}.sorted.bam" into sorted_bam_ch

    script:
    sample_id = sam_file.baseName
    """
    samtools view -bS ${sam_file} | samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam
    """
}

// Step 5: Generate alignment statistics with SAMtools
process alignmentStats {
    tag "${sample_id}"

    input:
    path bam_file from sorted_bam_ch

    output:
    path "${sample_id}.stats.txt" into stats_ch

    script:
    sample_id = bam_file.baseName
    """
    samtools flagstat ${bam_file} > ${sample_id}.stats.txt
    """
}

// Step 6: Organize final results into output directory
process organizeResults {
    input:
    path stats_files from stats_ch

    output:
    path("${params.outputDir}")

    script:
    """
    mkdir -p ${params.outputDir}
    mv *.sorted.bam *.stats.txt *.fastqc.zip ${params.outputDir}
    """
}

workflow {
    // Run FastQC in parallel, trim reads, then process with HISAT2 and SAMtools
    organizeResults(stats_ch)
}

  ```
</details>

---

