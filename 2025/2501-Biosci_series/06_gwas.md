# Introduction to GWAS

<br>
<br>
<p align="center">
    <img src="https://en.wikipedia.org/wiki/Genome-wide_association_study#/media/File:Manhattan_plot_from_a_GWAS_of_kidney_stone_disease.png" width="550">
</p>
<br>

Figure from [Genomics, Domestication, and Evolution of Forest Trees, Sederoff *et al.*, 2010](https://www.researchgate.net/publication/43073276_Genomics_Domestication_and_Evolution_of_Forest_Trees)

---
>[!important]
> :clock1: **Schedule**
> - 2:00pm-2:10pm: Welcome and introdution to topic (What is GWAS?)
> - 2:10pm-2:25pm: General pipeline for GWAS
> - 2:25pm-end: Example Analysis

>[!important]
> :heavy_exclamation_mark: **Requirements**
> - Basic command line knowledge
>- Access to a [Terminal](https://en.wikipedia.org/wiki/Unix_shell)
>    - Unix and Mac users already have access to the Terminal
>    - Windows users can use either [PowerShell](https://en.wikipedia.org/wiki/PowerShell) or the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install)
> - A registered CyVerse account (Register for a CyVerse account)

>[!important]
> :white_check_mark: **Expected Outcomes**
> - Understanding the general GWAS pipeline

<br>

---
---

>[!important]
> :heavy_exclamation_mark: **Materials for today's workshop are made available by the Marees *et al.*, [A tutorial on conducting genome-wide association studies: Quality control and statistical analysis](https://pubmed.ncbi.nlm.nih.gov/29484742/), *Int J Methods Psychiatr Res*, 2018. Our thanks and appreciation go to the rightful creators of the materials.**
> 
> The original material covers multiple hours of content; This workshop introduces the learner to what GWAS is conducted and the approaches required for experiments surrounding the topic.
>
> Additionally, this workshop adds an executable Jupyter Notebook learners can follow in order to carry out the entire pipeline.
>
> For a more in depth dive, visit the original material at **https://github.com/MareesAT/GWA_tutorial**.

>[!note]
> This workshop requires [plink](https://www.cog-genomics.org/plink/), an open-source tool designed for large scale genetic analysis (e.g., GWAS), and R.

---

## What's a GWAS?

<br>
<br>
<p align="center">
    <img src="https://avikarn.com/image/gwas/gwasflow.PNG" width="550">
</p>
<br>

From: *Genome-wide Association Study (GWAS) in TASSEL (GUI)*, from https://avikarn.com/2019-07-22-GWAS/ (honestly, a pretty good GUI friendly tutorial!)

Genome-Wide Association Studies (GWAS) are a method used in genetics to identify genetic variants associated with specific traits or diseases. The aim of GWAS is to identify **single nucleotide polymorphisms (SNPs)** of which the allele frequencies vary systematically as a function of phenotypic trait values (e.g., between cases with schizophrenia and healthy controls, or between individuals with high vs. low scores on neuroticism). Identification of trait-associated SNPs may subsequently reveal new insights into the biological mechanisms underlying these phenotypes.

The study requires a large number of data:

**DNA Collection** – DNA samples from (at least) two groups of people (or other organisms, to stay on topic we are going to be using humans as our species of interest):

- One group with the trait/disease (cases)
- One group without the trait/disease (controls)

**Genotyping** – each collected DNA needs to be genotyped in order to find single nucleotide polymorphisms (SNPs) across the population.

Each SNP is tested to see if it appears more frequently in cases than in controls. If a SNP is significantly more common in cases, it may be associated with the trait/disease.

>[!IMPORTANT]
> Identified SNPs don’t necessarily cause the trait but may be *linked* to genes influencing it. Further research is needed to determine their biological impact (*knock out! ... or knock in!*)

There are obvious limitations to GWAS:

- GWAS finds associations, not direct causation
- Some traits are influenced by many genes and environmental factors
- Rare genetic variants may be missed
- Population bias
- Missing Heritability

These are important points to keep in mind when designing a GWA experiment.

### How does GWAS compare to QTL analysis?

Quantitative Trait Locus (QTL) analysis, which we have seen last week, is another approach used to study the genetic basis of traits, particularly in model organisms (e.g., plants, mice). Unlike GWAS, which scans genomes in large populations, QTL analysis is often performed in controlled crosses of genetically distinct parents.

| Feature |	GWAS | QTL Analysis |
| :---: | :---: | :---: |
| Study Type | Observational (in natural populations) | Experimental (crossing of selected individuals) |
| Resolution | High (identifies SNPs) | Lower (identifies larger chromosomal regions) |
| Trait Type | Common, complex traits | Quantitative traits in controlled settings |
| Data Source | Natural human populations | Laboratory or agricultural species | 
| Population | Large, unrelated individuals | Genetically controlled crosses | 
| Limitations | Detects associations but not causation, small effect sizes, missing heritability | Lower resolution, requires controlled breeding | 

GWAS and QTL analysis can go hand in hand: **GWAS is best for studying complex traits in large populations, whilst QTL is best for controlled studies and finding major-effect genes.**

1. **QTL mapping can validate GWAS findings**:
    - GWAS identifies SNPs associated with a trait in natural populations, but these associations are often weak and require further validation.
    - QTL analysis, which uses controlled crosses, can help confirm whether a genomic region linked to a GWAS-identified SNP truly influences the trait.
    - Example: In Arabidopsis thaliana, a GWAS might find an SNP linked to drought tolerance. A QTL study in controlled crosses can determine if this SNP is consistently associated with drought resistance across different environments.
2. **QTL mapping helps identify causal variants**
    - GWAS often finds SNPs in non-coding regions, making it unclear whether they directly affect the trait.
    - QTL analysis provides a broader genomic region linked to the trait, which can then be fine-mapped to find the causal variant.
    - Example: In mice, a GWAS might link an SNP to obesity, but QTL analysis in controlled breeding experiments can narrow down the gene responsible.
3. **GWAS has "higher resolution", QTL has "higher power"**
    - GWAS scans millions of SNPs across the genome, offering high resolution but sometimes weak signals.
    - QTL studies focus on specific genomic regions and usually identify stronger genetic effects, though at lower resolution.
    - Combining them allows researchers to identify both small-effect and large-effect variants.
4. **Using QTL to improve GWAS in complex traits**
    - Many traits (e.g., height, yield, disease resistance) are polygenic, meaning they result from multiple genes interacting.
    - QTL studies can identify major-effect loci, while GWAS can help detect smaller-effect loci distributed across the genome.
    - Example: In livestock genetics, QTL mapping is used to find major genes controlling growth, while GWAS helps detect polygenic effects from many small SNPs.

---

## The GWAS pipeline

This is not *set in stone*, there isn't a single catch-them-all pipeline. This workshop follows the general outline of the Marees tutorial and from a bird eye view, these are the "steps":

1. Quality Control (QC): essential to produce reliable GWAS outcome (*like any other research!*).
2. Association Testing: inferring if there is a connection between SNPs to phenotypic traits.
3. Vizualization of results (through a Manhattan plot).

Within QC, the original author performs the following:

- Filtering SNPs & Samples: Removing low-quality data based on missingness, **minor allele frequency (MAF)**, and **Hardy-Weinberg Equilibrium (HWE)**.
- Sex Check: Ensuring consistency between reported and genetic sex.
- Population Stratification: Addressing ethnic differences using **Multidimensional Scaling (MDS)** and **Principal Component Analysis (PCA)**.
- Heterozygosity Check: Identifying sample contamination.
- Relatedness Check: Removing cryptic relatedness using **identity-by-descent (IBD)** measures.

Most of these steps are covered in "Step 1" and "Step 2" in the tutorial.