# RNA-Seq and Differential Expression Analysis using DESeq2

<br>
<br>
<p align="center">
    <img src="https://raw.githubusercontent.com/CosiMichele/workshop-materials/refs/heads/main/2024/2408-Biosci_series/assets/deseq2_volcano.png" width="750">
</p>
<br>

---
>[!important]
> :clock1: **Schedule**
> - 2:00pm-2:15pm: Welcome and introdution to topic (What is DEA?)
> - 2:15pm-2:25pm: Experiment Design and Data Preparation (What's upstream?)
> - 2:25pm-2:00pm: Differential Expression Analysis with DESeq2 and visualizations

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Understanding a general DEA Workflow
> - Revisiting previously discussed tools and techniques
> - Understanding the basics of DESeq2 

<br>

---
---

## Topic overview

QTL (Quantitative Trait Locus) mapping is the process of identifying genomic regions associated with variation in quantitative traits. Unlike single-gene traits, quantitative traits (e.g., height, yield, disease resistance) are controlled by multiple genes and environmental factors. QTL mapping uses genetic markers to locate these regions.

Quantitative trait mapping is used in biomedical, agricultural, and evolutionary studies to find causal genes for quantitative traits, to aid crop and breed selection in agriculture, and to shed light on natural selection. Examples of quantitative traits include cholesterol level, plant yield, or egg size, all of which are continuous variables. 

The goal of quantitative trait locus (QTL) analysis is to identify genomic regions linked to a phenotype, to map these regions precisely, and to define the effects, number, and interactions of QTL.

QTL analysis in experimental crosses requires two or more strains that differ genetically with regard to a phenotype of interest. Genetic markers, such as SNPs or microsatellites, distinguish between parental strains in the experimental cross. 

Markers that are genetically linked to a phenotype will segregate more often with phenotype values (high or low values, for example), while unlinked markers will not be significantly associated with the phenotype. The markers themselves might be associated with the phenotype but are not causal. Rather, markers may be associated with the phenotype through linkage to nearby QTL. They serve as signposts indicating the neighborhood of a QTL that influences a phenotype. Covariates such as sex or diet can also influence the phenotype.