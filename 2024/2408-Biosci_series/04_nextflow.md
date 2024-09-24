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
> - 3:00pm-3:10pm: Welcome and introdution to topic (What is Nextflow?)
> - 3:10pm-3:20pm: The importance of pipeline management and the direction of bioinformatics
> - 3:20pm-3:50pm: Understanding Nextflow Syntax and pipeline execution
> - 3:50pm-4:00pm: Resources and closing remarks

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Understand the concept and significance of workflow management in bioinformatics.
> - Identify the advantages of using Nextflow for data-driven analysis.
> - Navigate the Nextflow syntax and structure, including processes, channels, and workflows.
> - Break down a practical Nextflow script and explain its components.
> - Execute a basic RNA-Seq analysis pipeline, integrating tools like FastQC, Salmon, and MultiQC.

<br>

---
---

## Topic overview

This workshop is designed to introduce participants to Nextflow, a powerful workflow management system that enables the creation, execution, and sharing of reproducible computational pipelines in bioinformatics. Participants will learn the fundamentals of Nextflow, including its importance in modern research, its unique features, and how it compares to other workflow management tools.

> [!IMPORTANT]
> This workshop uses is an adaptation of the [Sequera Nextflow tutorial](https://github.com/seqeralabs/nextflow-tutorial), shortened in order to fit time limits. **Nextflow is something that cannot be learned in 1 hour or less**, <u> *but its components are worth understanding* </u>. For this tutorial, the organism used is *G. gallus*, the red junglefoul (the o.g. chicken :chicken:). All files (transcriptome, fastq files) are made available on CyVerse. 

> [!NOTE]
> **Overview of Pipeline**
> As this can be classified as an analysis pipeline, it is important for us to highlight what the major steps are:
> ```
> [Indexing trascriptome]                   
>         ↓ 
> [Perform QC]                        
>         ↓
> [Perform quantification]            
>         ↓
> [Create MultiQC report] 
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
    gocmd get --progress /iplant/home/shared/cyverse_training/datalab/biosciences/nextflow_tutorial/
    ```

    This will copy:
    
    - The *G. gallus* transcriptome (`transcriptome.fa`)
    - Gut pair-end reads (`gut_1.fq`, `gut_2.fq`)
    - The NextFlow code required for the workshop.
    - Other read files for exercise purposes (liver, lung) 

You are now ready for the workshop!

<details>
  <summary>Click here for the raw NextFlow code</summary>

```  
nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

println """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/* 
 * create a transcriptome file object given the transcriptome string parameter
 */
transcriptome_file = file(params.transcriptome)

/* 
 * Define the `index` process that creates a binary index 
 * from the transcriptome file
 */
process index {
    conda 'bioconda::salmon'

    input:
    path transcriptome_file

    output:
    path 'index'

    script:
    """
    salmon index --threads ${task.cpus} -t ${transcriptome_file} -i index
    """
}

/* 
 * Define a channel for paired reads and set the channel name to `read_pairs_ch`
 */
Channel.fromFilePairs(params.reads)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

/*
 * Define the `quantification` process that runs Salmon quantification
 */
process quantification {
    conda 'bioconda::salmon'

    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path "${pair_id}_quant"

    script:
    """
    salmon quant --threads ${task.cpus} --libType=U -i ${index} -1 ${reads[0]} -2 ${reads[1]} -o ${pair_id}_quant
    """
}

/*
 * Define the `fastqc` process that performs quality control on reads
 */
process fastqc {
    conda 'bioconda::fastqc'
    tag "FASTQC on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

/*
 * Define the `multiqc` process that generates a report from fastqc and quant outputs
 */
process multiqc {
    publishDir params.outdir, mode: 'copy'
    conda "bioconda::multiqc"

    input:
    path fastqc_files
    path quant_files

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

/* 
 * Workflow definition that orchestrates the execution of processes 
 */
workflow {
    /*
     * Pass the transcriptome file to the `index` process
     */
    index_ch = index(transcriptome_file)

    /*
     * Run the processes in the desired order and collect outputs
     */
    fastqc_ch = fastqc(read_pairs_ch)

    quant_ch = quantification(index_ch, read_pairs_ch)

    multiqc(fastqc_ch, quant_ch)
}

  ```
</details>

---

