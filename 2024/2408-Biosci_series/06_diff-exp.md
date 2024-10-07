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
> - 3:00pm-3:10pm: 
> - 3:10pm-3:45pm:
> - 3:45pm-4:00pm:

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

This workshop focuses on **RNA-seq analysis** and **differential expression analysis (DEA)**, core techniques for examining gene expression. 

We’ll start by revisiting what a typical RNA-seq analysis workflow looks like, taking raw sequencing data and processing it to uncover genes that show significant changes in expression between conditions. 

In the bulk of the session, you’ll see how to handle RNA-seq data using R and popular packages like **DESeq2** for analyzing differential gene expression, as well as **ggplot2** and **EnhancedVolcano** for visualizing results in meaningful ways. We’ll also cover some of the practical steps, including **filtering low-count genes**, applying statistical tests to find **significant differences**, and generating visual summaries with **PCA plots**, **heatmaps**, and **volcano plots**.

> [!IMPORTANT]
> This workshop uses is an adaptation of the "RNA-seq and Differential Expression" Workshop by Texas A&M HPC group ([ACES](https://hprc.tamu.edu/aces/)) for the University of Arizona, first held in April 2024. The workshop was held over a period of 3 hours; This Workshop willl focus on the DESeq2 section in order to respect time limits. 
>
> The model organism used is the house mouse ([Mus musculus, GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/)); RNA-seq data was obtained by Texas A&M, trimmed and aligned by us.
> 
> For this workshop, we have made available just the files needed for the DESeq2 analysis.
---
---

## What is Differential Expression Analysis?

<p align="center">
    <img src="https://www.sc-best-practices.org/_images/differential_gene_expression.jpg" width="750">
</p>

Image source: [Single Cell Best Practices](https://www.sc-best-practices.org/conditions/differential_gene_expression.html).

**Differential gene expression analysis** (DEA) attempts to infer genes that are statistically significantly over- or underexpressed between any compared groups (commonly between healthy and condition per cell type). The outcome of such an analysis could be genesets which effect and potentially explain any observed phenotypes. These genesets can then be examined more closely with respect to, for example, affected pathways or induced cell-cell communication changes.

In simple terms: *DEA is a method to identify which genes are more or less active under different conditions.*

**Why is it important?**

This analysis is crucial because it helps scientists understand how genes function in various situations. For instance, if a gene is very active in cancer cells but not in normal cells, it might be important for cancer. This understanding can lead to new treatments or a better grasp of diseases.

**How does it work?**

The process involves collecting data on gene expression levels from different samples. Then, statistical tests determine if the differences in expression are significant, meaning they are likely not due to random chance. This is similar to comparing test scores between two classes to see if one class performed better. 

Let's break DEA down into 4 steps:

1. **Measure Expression**: First, RNA-seq provides a snapshot of gene activity by counting how many RNA molecules (reads) correspond to each gene in the samples. This read count reflects the gene’s expression level in each sample.
2. **Identify Differences**: DEA compares these expression levels across groups. For each gene, it assesses whether the average expression in one condition (e.g., treated) is statistically higher or lower than in the other condition (e.g., control).
3. **Apply Statistical Testing**: To decide if the differences are meaningful, DEA uses statistical tests to see if the observed changes are likely due to the conditions rather than random variation. This analysis outputs measures like p-values and adjusted p-values (accounting for multiple testing) to indicate the significance of these changes.
4. **Interpret Fold-Changes**: If a gene shows a large fold-change with a low p-value, it’s considered "differentially expressed" because its expression is significantly different between the conditions. These genes may then be candidates for further research, as their expression differences might reveal important biological insights about the system under study.

>[!NOTE]
> **What's a Fold-Change?**
> 
> **Fold change** is a measure describing how much a quantity changes between an original and a subsequent measurement. It is defined as the ratio between the two quantities; for quantities A and B the fold change of B with respect to A is B/A.
>
> Fold-Changes are usually observed via Volcano plots, post Log2 transformation: instead of working directly with the raw fold-change, we take the logarithm (base 2) of this value. This has two main benefits:
> - Symmetry: Log2 fold-change makes it easier to interpret increases and decreases equally. For example:
>   - A log2 fold-change of +1 means the gene expression doubled ($2^1 = 2$).
>   - A log2 fold-change of -1 means the gene expression halved ($2^-1 = 0.5$).
> - Range Reduction: It reduces extreme values, making data easier to handle and visualize.
>
> Example:
> <p align="center">
>    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/7/74/Volcano_eg.jpg/1280px-Volcano_eg.jpg" width="550">
>   <figcaption> Volcano plot showing metabolomic data. The red arrows indicate points-of-interest that display both large magnitude fold-changes (x axis) and high statistical significance (-log10 of p value, y axis). The dashed red line shows where p = 0.05 with points above the line having p < 0.05 and points below the line having p > 0.05. This plot is colored such that those points having a fold-change less than 2 (log2 = 1) are shown in gray. </figcaption>
> </p>
> Image source: [Wikipedia](https://en.wikipedia.org/wiki/Fold_change). 


- [content 1](#section-1-subsection-1) 
- [content 2](#section-1-subsection-2)

### Section 1 Subsection 1

---

### Section 1 Subsection 2

---
---

```
###########################################################
# 0. LOAD NECESSARY LIBRARIES AND SET WORKING DIRECTORY  ##
###########################################################

library(ggplot2)                                          # For data visualization
library(pheatmap)                                         # For heatmap visualization
library(DESeq2)                                           # For differential gene expression analysis
library(EnhancedVolcano)                                  # For volcano plot visualization of DESeq2 results

# Set working directory
setwd("/home/jovyan/data-store/r_tutorial/counts")

###########################################################
# 1. LOADING FILES, CREATING DATASET, FILTER LOW COUNTS ###
###########################################################

# Check files in the directory
list.files()               # Lists all files in the working directory
system("cat sampleTable.csv") # Displays contents of sampleTable.csv

# Load sample table with sample metadata and conditions
sampleTable <- read.csv("sampleTable.csv", header=TRUE)   # Reads CSV file
sampleTable <- as.data.frame(sampleTable)                 # Ensures it is in data frame format
sampleTable$condition <- factor(sampleTable$condition)    # Sets 'condition' as a factor
sampleTable                                               # Displays the sample table

# Create DESeq2 dataset using counts and design information
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = ".",
                                  design = ~ condition)   # Design formula for the experiment


# Show summary of DESeq2 dataset
dds

# Filter out genes with low counts across all samples
keep <- rowSums(counts(dds)) >= 10                        # Only retain genes with counts ≥ 10
dds <- dds[keep,]                                         # Updates dds to only include high-count genes

###########################################################
# 2. DIFFERENTIAL EXPRESSION ANALYSIS #####################
###########################################################

# Perform differential expression analysis
dds <- DESeq(dds)                                         # Runs DESeq2 pipeline for DE analysis
res <- results(dds)                                       # Extracts results from the DESeq2 analysis
res                                                       # Shows results with fold changes, p-values, etc.

# Count significantly differentially expressed genes
sum(res$padj <= 0.05, na.rm = TRUE)                       # Counts genes with an adjusted p-value ≤ 0.05

# Filter for significantly differentially expressed genes and save to CSV
sigGenes <- res[which(res$padj < 0.05), ]                 # Filters significant genes
write.csv(sigGenes, "Differentrially_Expressed.csv", 
          row.names = TRUE)                               # Saves results

###########################################################
# 3. VISUALIZATION: #######################################
#   3.1 PCA PLOT    #######################################
###########################################################

# Regularized log transformation to stabilize variance
logTran <- rlog(dds)               

# Calculate row variances and select top variable genes
rv <- rowVars(assay(logTran))                             # Calculates variance for each gene across samples
select <- order(rv, decreasing = TRUE)[seq_len(min(100, length(rv)))] # Selects top 100 variable genes

# Perform Principal Component Analysis (PCA) on selected genes
PCA <- prcomp(t(assay(logTran)[select, ]), scale = FALSE) # Runs PCA
summary(PCA)                                              # Shows PCA summary

# Calculate and prepare variance explained for plotting
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2), 1)    # Calculates percent variance for each PC
ggPCA_out <- as.data.frame(PCA$x)                         # Converts PCA output to data frame
ggPCA_out <- cbind(ggPCA_out, sampleTable)                # Adds sample information to PCA data

# Plot PCA with ggplot2
ggplot(ggPCA_out, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size=4) +
  labs(x = paste0("PC1 - variance explained: ", percentVar[1], "%"),
       y = paste0("PC2 - variance explained: ", percentVar[2], "%")) +
  theme_bw()                                              # Adds basic theme

###########################################################
# 3. VISUALIZATION:  ######################################
#   3.2 VOLCANO PLOT ######################################
###########################################################

# Create volcano plot for differential expression results
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,                          # Significance threshold
                FCcutoff = 1.0,                          # Fold-change threshold
                pointSize = 3.0,                         # Point size
                labSize = 4.0,                           # Label size
                colAlpha = 1/2,                          # Point transparency
                drawConnectors = FALSE,                  # Disable connectors for labels
                legendPosition = "right")                # Position legend

##########################################################
# 3. VISUALIZATION:  #####################################
#   3.3 HEATMAP      #####################################
##########################################################


# Filter top significant genes with strong log2 fold-change for heatmap
resorted_deresults <- res[order(res$padj),]              # Sort results by adjusted p-value
sig <- resorted_deresults[!is.na(resorted_deresults$padj) &
                            resorted_deresults$padj < 0.05 &
                            abs(resorted_deresults$log2FoldChange) >= 6.5, ]

selected <- rownames(sig)                                # Get names of selected genes
selected                                                 # Display selected genes

# Plot heatmap of selected gene expression values
transformed_readcounts <- normTransform(dds)             # Apply variance-stabilizing transformation
pheatmap(assay(transformed_readcounts)[selected,],
         cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = TRUE,
         labels_col = colData(dds)$sampleName)           # Plots heatmap with clustering of samples and genes

```