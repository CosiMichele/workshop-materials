# Using Snakemake for Streamlining Data Analysis Pipelines

<br>
<br>
<p align="center">
    <img src="https://genomics.ed.ac.uk/wp-content/uploads/2024/08/snakemake_logo.png" width="650">
</p>
<br>

Image credits: [Snakemake logo from the University of Edinburgh's Snakemake workshop](https://genomics.ed.ac.uk/event/snakemake/).


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

## What is Nextflow?

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/nf-core/rnaseq/3.15.1/docs/images/nf-core-rnaseq_metro_map_grey_animated.svg" width="750">
</p>
<br>

Image source: [nf-core/rnaseq](https://nf-co.re/rnaseq/3.15.1/).

Written in [Groovy](https://en.wikipedia.org/wiki/Apache_Groovy), **Nextflow** is an open-source workflow management system that allows researchers to write and execute data-driven computational workflows. It is designed to handle complex workflows in a scalable and reproducible manner.

Although Nextflow is heavily used in the bioinformatics field, its processes are used in other fields including cloud computing, data science, image analysis and machine learining.

### Why Nextflow?

- **Portability**: Nextflow can run on various platforms, including local machines, clusters, and cloud environments, without modification to the workflow script.
- **Scalability/Parallel Execution**: It automatically parallelizes tasks, optimizing resource utilization and speeding up the analysis.
- **Reproducibility**: Nextflow enables the creation of workflows that can be easily shared and rerun, ensuring that analyses can be replicated by others.

Uniquely to Nextflow is its ability to manage piplelines using **Channels**: these are processes that "activate" and remain active until the workflow is turned off. Other pipeline managers fail at this as these are reliant on outputs and inputs. *Example: when you turn the water on and off, the pipes don't disappear.*

Additionally, Nextflow has a robust community that create, share and maintain their pipelines at [nf-co.re](https://nf-co.re/) with extremely well documented work, responsive support and in-person and online events.

Nextflow has seen a surge in industry use, with[ **Nanopore**](https://nanoporetech.com/) officially creating [Nextflow pipelines through "EPI2ME Workflows"](https://labs.epi2me.io/wfindex/). A great example of one of these pipelines is [nanoseq](https://nf-co.re/nanoseq/3.1.0/): an analysis pipeline for Nanopore DNA/RNA sequencing data that can be used to perform basecalling, demultiplexing, QC, alignment, and downstream analysis.


---

## Understanding the Snakemake Components

Snakemake's main features revolve around **Rules**, with each rule having an **input**, **output**, (facultatively) a **conda** directive, and a **shell** OR a **script**. Additionally, Snakemake has a *special* rule: the **All** rule. This ensures that the workflow is done by checking on the last expected output. **Wildcards**, **DAG** (Directed Acyclic Graph) are two other important components to any Snakemake workflow.

### Rules

Rules are the building blocks of a Snakemake workflow. Each rule defines how to create an output file from input files using a shell command.

Example:

```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

###  Rule All

The `all` rule specifies the final outputs of the entire workflow. This is the "target" rule that tells Snakemake what the end products of the workflow should be. It can be placed in the beginning or the end of the workflow (in the early iterations of Snakemake this was something that was strictly needed at the beginning of the workflow).

Example:

```
rule all:
    input:
        "plots/quals.svg"
```

Notice how in the example above, the expected output is defined as `input`. This is because of how Snakemake works, my having `rule all` "check" on the final output by ingesting it (thus requiring to be defined as input).

### [Directed Acyclic Graph (DAG)](https://en.wikipedia.org/wiki/Directed_acyclic_graph)

Snakemake automatically constructs a [Directed Acyclic Graph (DAG)](https://en.wikipedia.org/wiki/Directed_acyclic_graph) from the defined rules. The DAG represents the dependencies between rules, showing which rules need to be run in sequence and which can run in parallel. The DAG ensures that the workflow executes only the necessary steps, allowing Snakemake to avoid re-executing tasks if their outputs are up-to-date.

```

```

###  

```

```

Here is the [example pipeline from the snakemake official documentation](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html): 
```
SAMPLES = ["A", "B"]


rule all:
    input:
        "plots/quals.svg"

rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"

> [!Note]
>
> **Construction Worksite Analogy**
>
> Imagine Snakemake as a construction manager overseeing the construction of a large building. Each part of the building represents a rule, and the materials, tools, and workers represent the inputs, outputs, and shell commands needed to complete the project.
> <br>
> <br>
> <p align="center">
>     <img src="https://media0.giphy.com/media/v1.Y2lkPTc5MGI3NjExNXNzbWZkd3I5bmZraHg0NjZqaDQwNTVwdXlsdHIwYnpidWphamQ0biZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/rKmwe2hfFIcdq/200.webp" width="350">
> </p>
> <br>
> 
> **Rules are Construction Tasks**
> 
> - Each rule in Snakemake is like a specific construction task. It defines the materials (inputs) needed to complete the task, the finished part of the building (output), and the work instructions (shell command) to get the job done.
> - Example: Pouring the foundation requires cement, water, and gravel, which results in a finished concrete base.
>
> **DAG is the Construction Blueprint**
> 
> - The DAG (Directed Acyclic Graph) is like the blueprint that the construction manager follows. It maps out the entire project, showing how each task (rule) fits into the larger plan, ensuring nothing is missed and that tasks are completed in the right order.
> - Example: The blueprint shows that electrical wiring should be installed after the walls are built but before the drywall goes up.
> 
> **Dependencies are Construction Phases**
>
> - Just like you can’t install the roof before building the walls, some tasks depend on others being completed first. Snakemake ensures that each task happens in the right order.
> - Example: before building walls, it is important to have the foundations finished.
> 
> **The "All" Rule is the Completed Building**
> 
> Just as the manager ensures that every part of the building is finished before calling it complete, the "all" rule in Snakemake guarantees the entire workflow is executed and all final outputs are produced.
> Example: the completed building!
>
> ---
>
> There are other processes that are required in Snakemake that are missing from this larger picture:
>
> - **Parallelization**: multiple teams working on different aspect of a task (one team pours concrete, the other lays it).
> - **Resources Allocation (& Conda Environments)**: specialized tools and workers that are accessed when and where needed (an electrician is potentially not required when first levelling ground).

---

## Running the Example Pipeline

<br>
<br>
<p align="center">
    <img src="https://combine-lab.github.io/salmon/images/SalmonLogo.png" width="450">
</p>
<br>

Image source: [salmon aligner](https://combine-lab.github.io/salmon/about/) (no joke, this is their logo).

In this tutorial we are running a curated Nextflow workflow. The workflow will index the provided genome, align sequences (1/3), run QC and output a summary. 

The alignment is carried out using [**Salmon**](https://combine-lab.github.io/salmon/about/), used for quantifying the expression levels of transcripts from RNA-Seq data. Unlike traditional aligners that map reads to a reference genome, Salmon directly estimates transcript abundance by utilizing the concept of pseudo-alignment (first mapping reads to a database of known short sequences, and then using these mapped positions to infer the position of each read on the target genome/transcriptome).

FastQC is used to do QC and MultiQC to summarize FastQC and Salmon outputs.

To execute the provided example, change directory to the `nextflow_tutorial` folder (`cd nextflow_tutorial`) and execute:

```
nextflow transcriptome_qc.nf
```

The output should be similar to this:

```
 N E X T F L O W   ~  version 24.04.4

Launching `transcriptome_qc.nf` [admiring_knuth] DSL2 - revision: 3004fdc561

R N A S E Q - N F   P I P E L I N E    
===================================
transcriptome: /home/jovyan/data-store/nextflow_tutorial/data/ggal/transcriptome.fa
reads        : /home/jovyan/data-store/nextflow_tutorial/data/ggal/gut_{1,2}.fq
outdir       : results

executor >  local (4)
[e1/6abbeb] process > index                  [100%] 1 of 1 ✔
[f4/d72072] process > fastqc (FASTQC on gut) [100%] 1 of 1 ✔
[24/30d512] process > quantification (1)     [100%] 1 of 1 ✔
[98/1c93ee] process > multiqc (1)            [100%] 1 of 1 ✔
```

### Breaking Things Down

```
nextflow.enable.dsl=2
```

- **Enables DSL2**: we are enabling Nextflow's second domain-specific language (DSL2), which allows for more structured and modular workflows.

---

```
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"
```

- These lines define input parameters for the workflow, specifying file paths for the RNA-Seq reads and the transcriptome.

---

```
println """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
```

- This block prints the title and input parameters to the console, providing immediate feedback about the configuration.

---

```
transcriptome_file = file(params.transcriptome)
```

- This line creates a file object for the transcriptome, making it easier to reference later in the script.

---

```
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
```

- This is the indexing process, which creates a binary index for the transcriptome using Salmon.
-  It uses a single directive (`conda 'bioconda::salmon'`): it specifies that the process should run in the `bioconda::salmon` environment (already present).
- Input and Output: It takes the transcriptome as input and produces an output named `index`.

---

```
Channel.fromFilePairs(params.reads)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }
```

- This block creates a channel from paired FASTQ files based on the specified pattern and checks for the existence of files.

---

```
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
```

- This process quantifies the expression levels of transcripts using the Salmon quantification tool (`salmon quant`).
- As input takes the index created in the index process and paired reads; as output it produces quantification results named after the sample IDs.

---

```
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
```
- This process runs FastQC on the input reads to assess their qualit and  generates a log directory for the FastQC results.

---

```
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
```

- This process runs MultiQC to aggregate the results from FastQC and Salmon quantification and produces an the MultiQC HTML report summarizing the results.

---

```
workflow {
    index_ch = index(transcriptome_file)
    fastqc_ch = fastqc(read_pairs_ch)
    quant_ch = quantification(index_ch, read_pairs_ch)
    multiqc(fastqc_ch, quant_ch)
}
```

- The workflow block defines the order of execution for the processes, linking the outputs of one process to the inputs of another, Establishing the flow of data throughout the pipeline.

---

We can now visualize the output, in `results`.

---

## Closing Remarks and Resources

<br>
<br>
<p align="center">
    <img src="https://upload.wikimedia.org/wikipedia/commons/2/23/Red_jungle_fowl.png" width="450">
</p>
<br>

Image source: [wikipedia](https://en.wikipedia.org/wiki/Red_junglefowl).


This is **BY NO MEANS** something simple to grasp. If anything, only a fool would try to explain Nextflow in 1 hour (I am, indeed, that fool).

Nextflow is complex, and it has a steep learning curve. There are things we did not cover:
- What's the `work` folder created?
- How can I check logs?
- Error handling??!?
- What's a `resume`?

Learning Nextflow is, however, extremely rewarding:
- Extremely well integrated with containers.
- Support for SLURM and PBS HPC systems.
- "Easily" sharable workflows.
- Honestly, pretty good on your CV.

Here are some great Nextflow resources for you to learn more:

- [General Documentation](https://www.nextflow.io/docs/latest/overview.html)
- [Nextflow training](https://training.nextflow.io/hands_on/)
- [Find well supported and extremely powerful pipelines at nf-co.re](https://nf-co.re/)
- [rna-seq specific pipeline using the aligner of your choice](https://nf-co.re/rnaseq/3.15.1/)
- [The full extended tutorial this was inspired by (note: this might only run on previous iterations of Nextflow, using dsl1)](https://github.com/seqeralabs/nextflow-tutorial)
- [CyVerse Webinar on NextFlow](https://www.youtube.com/watch?v=umd0ikOCxgM) (shoutout to Conner Copeland!)