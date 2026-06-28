# R-SEQ-Flow v3
## RNA-seq Pipeline — Automatic Quality Control, Pseudo-Alignment Quantification, Differential Expression Analysis and Results Visualization

![License](https://img.shields.io/badge/License-MIT-green.svg)
![Platform](https://img.shields.io/badge/Linux-Ubuntu-orange)
![Language](https://img.shields.io/badge/Bash%20%7C%20R-blue)
![Version](https://img.shields.io/badge/version-v3.0-success)
![Maintained](https://img.shields.io/badge/Maintained-Yes-brightgreen)

---
## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running R-SEQ-Flow](#running-r-seq-flow)
- [Smart Caching & Resumability](#smart-caching--resumability)
- [Case Study Results](#case-study-results)
- [Output Structure](#output-structure)
- [Configuration](#configuration)
- [Input Format](#input-format)
- [Output Interpretation](#output-interpretation)
- [Methods](#methods)
- [License](#license)
- [Author](#author)
- [Related Resources](#related-resources)




---

---

## Overview

**R-SEQ-Flow** is a complete, **production-ready RNA-seq processing pipeline** that automates the entire workflow from raw sequencing reads (NCBI SRA) to differential expression analysis and results visualization.

### What It Does

- **SRA Download & Conversion**: Automatic FASTQ extraction from NCBI SRA with 3-level caching
- **Quality Control**: Pre/post-trimming FastQC analysis with detailed HTML reports
- **Adapter Trimming**: fastp with quality filtering and read length validation
- **Pseudo-Alignment Quantification**: Kallisto-based transcript abundance estimation
- **Differential Expression Analysis**: DESeq2 with gene-level aggregation, LFC shrinkage, and results visualizations
- **Resumable Execution**: Automatic checkpoint system—crash and resume without redundant re-processing
- **Shared Reference Architecture**: Download transcriptome and build Kallisto index once; reuse across unlimited experiments

### Features

 **11-step Bash pipeline** + **24-step R analysis**  
 **Smart 3-level cache system** (FASTQ → SRA → download)  
 **Per-sample and per-experiment checkpoints** for resumability  
 **Multi-threaded processing** (FastQC, fastp, Kallisto parallelized)  
 **Automatic PE/SE detection** from FASTQ files  
 **Comprehensive error logging** with per-tool logs  
 **Disk-space pre-flight checks**  
 **Plots generation** (PCA, volcano, heatmaps, MA plot, Dispersion)  
 **Full gene expression matrices** with apeglm LFC shrinkage  
 **Gene lists** (up-regulated, down-regulated, significant)  

---


## Prerequisites


Before running **R-SEQ-Flow**, ensure the following software is installed on your system.

### Operating System

- Linux (Ubuntu 20.04 or later recommended)
- Windows users can run the pipeline through Windows Subsystem for Linux (WSL2)

### Programming Languages

- Bash (GNU Bash 5.0 or later)
- R (version 4.3 or later)

### Required Bioinformatics Software

- SRA Toolkit
- FastQC
- Fastp
- Kallisto

### Required R Packages

- DESeq2
- tximport
- ggplot2
- pheatmap
- RColorBrewer
- EnhancedVolcano
- readr
- dplyr

### Recommended System Requirements

- Multi-core CPU
- 8 GB RAM minimum (16 GB recommended)
- At least 20 GB of free disk space

### Installation 

```bash
# 1. Clone repository
git clone https://github.com/abdullahakram013/R-SEQ-Flow.git
cd R-SEQ-Flow

# 2. Make executable
chmod +x scripts/R-seq.sh
./scripts/R-seq.sh

# 3. Run pipeline
bash scripts/R-seq.sh
```
---

## Running R-SEQ-Flow


All pipeline scripts are located in the `scripts/` directory.


R-SEQ-Flow consists of two main executable files.


| File | Purpose |
|------|---------|
| `R-seq.sh` | Main Bash pipeline that performs data download, quality control, trimming, pseudo-alignment, quantification, checkpoint management, and metadata generation. |
| `deseq2_analysis.R` | Performs downstream differential expression analysis, statistical testing, result export, and visualization using DESeq2. |

### Run the Complete Pipeline

bash
```
bash scripts/R-seq.sh
```


## Smart Caching & Resumability

### 3-Level Cache Architecture

**Level 1**: FASTQs cached → Skip all downloads  
**Level 2**: SRA file cached → Skip SRA download, convert locally  
**Level 3**: Nothing cached → Download + convert

### Crash Recovery

```bash
# If pipeline crashes at any Step :
bash scripts/R-seq.sh
# Enter same experiment name + SRA IDs
# Pipeline skips all completed steps (1-9), resumes at Step 10
```

---

## Case Study: RNA-seq Differential Expression Analysis (ASTHMA Dataset)

### ASTHMA Study (6 samples)

```bash
$ bash scripts/R-seq.sh

Experiment name : ASTHMA_EXPERIMENT
CONTROL SRA ID(s)  : SRR1039508 SRR1039512 SRR1039516
TREATED SRA ID(s)  : SRR1039509 SRR1039513 SRR1039517

[Processing 11 steps ...]

Continue to DESeq2 analysis? (yes/no): yes
```

### 11-Step Bash Pipeline

1. **Experiment Setup** — Metadata input validation
2. **Directory Creation** — Smart folder structure with resume detection
3. **Tool Installation** — `apt-get` install (runs once via checkpoint)
4. **SRA Download** — 3-level cache: FASTQ → SRA → download
5. **FastQC (raw)** — Pre-trimming quality reports
6. **Trimming** — fastp adapter removal and QC filtering
7. **FastQC (trimmed)** — Post-trimming validation
8. **Reference Download** — GENCODE v29 transcriptome (shared, once-only)
9. **Kallisto Index** — Build index from transcriptome (shared, once-only)
10. **Quantification** — Per-sample pseudo-alignment (Kallisto)
11. **Metadata** — Create DESeq2 sample sheet

### 24-Step R Analysis (Optional)

- Steps 1-7: Package loading, metadata parsing, Kallisto file validation
- Steps 8-11: Gene-level count matrix creation via tximport
- Steps 12-14: DESeq2 object creation, pre-filtering, normalization
- Steps 15-19: Differential expression testing, LFC shrinkage, results export
- Steps 20-24: Generating Plots for visualization (9 plots + 7 tables)

---

## Output Structure

```
ASTHMA_EXPERIMENT/
├── checkpoints/                    # Resumability flags
├── trimmed_data/                   # Trimmed FASTQ files
├── qc_results/                     # FastQC HTML reports
├── kallisto_output/SRR*/           # Pseudo-alignment (6 samples)
├── metadata/samples.txt            # DESeq2 sample sheet
├── logs/                           # Per-tool detailed logs
└── deseq2_results/                 # Differential expression outputs
    ├── all_genes_DESeq2_results.csv
    ├── significant_genes.csv       (padj < 0.05)
    ├── upregulated_genes.csv       (LFC > 1)
    ├── downregulated_genes.csv     (LFC < -1)
    ├── PCA_plot.png
    ├── Volcano_plot.png
    ├── Top50_heatmap.png
    ├── sample_distance_heatmap.png
    ├── MA_plot.png
    ├── dispersion_plot.png
    ├── Gene_expression_scatter.png
    ├── Gene_expression_histogram.png
    ├── Gene_expression_boxplot.png
    ├── summary_statistics.csv
    ├── FIGURE_CAPTIONS
    └── session_info.txt
```

---



**Shared Directory** (reused across all experiments):
```
SHARED/
├── fastq/                   # Downloaded FASTQs (78 GB for 6 samples)
├── reference/gencode.v29    # Transcriptome (750 MB)
├── kallisto_index/          # Index (built once, 2 GB)
└── checkpoints/             # Global flags
```

---

## Case Study Results

**Study**: Airway smooth muscle cells ± dexamethasone  
**Samples**: 3 control + 3 treated (paired-end, 75bp)  
**Runtime**: ~90 minutes (first run), ~2 minutes (resume)

```
Total genes in reference:       58,721
Total genes quantified:         17,450
Significantly DE genes (padj<0.05): 1,189
  Up-regulated (LFC > 1):         259
  Down-regulated (LFC < -1):      170
```

### Top DE Genes
```
Gene          log2FC   -log10(padj)   Description
DUSP1          +3.45      89          Dual-specificity phosphatase
FKBP5          +3.12      82          FK506 binding protein 5
HSPA1A         +2.89      73          Heat shock protein 70
```

---

## Configuration


```bash
THREADS=4              # CPU threads for parallel tools
BOOTSTRAP=100          # Kallisto bootstrap resamples
SHARED_DIR="/path/to/shared"  # Adjust to your system
```



---

## Input Format

**SRA IDs**: Space-separated NCBI SRA run identifiers  
Example: `SRR1039508 SRR1039512 SRR1039516`

Find SRA IDs: https://www.ncbi.nlm.nih.gov/sra

**Metadata** (auto-generated):
```
sample          condition    path
SRR1039508      control      /experiment/kallisto_output/SRR1039508/abundance.tsv
SRR1039509      treated      /experiment/kallisto_output/SRR1039509/abundance.tsv
```

---

## Output Interpretation

### CSV Files

**all_genes_DESeq2_results.csv** — All 17,450 genes
```csv
gene_id,gene_name,baseMean,log2FoldChange,lfcSE,pvalue,padj
ENSG00000223972,DDX11L1,0.52,0.18,1.42,0.898,0.989
```

**significant_genes.csv** — 1,189 genes with padj < 0.05
```csv
gene_id,gene_name,baseMean,log2FoldChange,lfcSE,pvalue,padj
ENSG00000030582,DUSP1,2145.3,3.45,0.142,1.2e-130,2.1e-89
```

**upregulated_genes.csv** — 259 genes with LFC > 1  
**downregulated_genes.csv** — 170 genes with LFC < -1


---

## Methods

**Quantification**: Kallisto v0.48+ (pseudoalignment)  
**Reference**: GENCODE v29 human transcriptome  
**Normalization**: DESeq2 size-factor estimation  
**Testing**: Wald test with Benjamini-Hochberg FDR correction  
**Shrinkage**: apeglm (Approximate Posterior Estimation)  
**QC**: FastQC + fastp  







## License

**MIT License** — Free for commercial and non-commercial use with attribution.  
See [LICENSE](LICENSE) file for full details.

---

## Author

**Mr.Abdullah Akram**  
Department of Bioinformatics and Biotechnology  | Government College University Faisalabad | Faisalabad, 38000, Pakistan 
- GitHub: [@abdullahakram013](https://github.com/abdullahakram013)
- Email: abdullahakram7652a@gmail.com


---

## Related Resources

- [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
- [GENCODE](https://www.gencodegenes.org/)
- [Kallisto](https://pachterlab.github.io/kallisto/)
- [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)

---



