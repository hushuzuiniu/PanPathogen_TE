# Pan-Pathogen Activation of Transposable Elements Shapes Host Immune Responses

This repository contains the computational analysis pipeline and code for the research project studying how various pathogens activate transposable elements (TEs) in host cells and their impact on immune responses.

## Authors

**Xiaoman Wei¹,³, Shu Hu²,³, Tianqi Xu², Chengjun Lu², Jie Cui¹,⁴***

¹ Shanghai Sci-Tech Innovation Center for Infection & Immunity, National Medical Center for Infectious Diseases, Huashan Hospital, Institute of Infection and Health, Fudan University, Shanghai, China  
² Fudan University, Shanghai, China  
³ These authors contributed equally  
⁴ Lead contact  

*Correspondence: jiecui@fudan.edu.cn

## Project Overview

Transposable elements (TEs) constitute a significant portion of mammalian genomes and play crucial roles in host-pathogen interactions. This project investigates how different pathogens (bacteria, viruses) systematically activate specific TE families and how these activations contribute to host immune responses.

## Repository Structure

### Data Processing Modules
- **`00.data/`** - Raw data storage and organization
- **`RNA-seq/`** - RNA sequencing data processing pipeline
- **`ATAC-seq/`** - ATAC sequencing data processing pipeline

### Analysis Modules
- **`01.Genomic_features_Gencode/`** - Genome annotation and feature classification (exon/intron/intergenic)
- **`02.DESeq2_analysis_Gene/`** - Differential gene expression analysis
- **`02.DESeq2_analysis_TE_loci/`** - Individual TE loci differential expression analysis
- **`02.DESeq2_analysis_TE_subfamily_feature/`** - TE subfamily analysis by genomic features
- **`03.Enrichment_analysis/`** - Functional enrichment analysis
- **`04.TE_evolution/`** - Evolutionary conservation analysis across 95 mammals
- **`05.TF_prediction/`** - Transcription factor binding motif prediction
- **`06.TE_GeneHancer/`** - TE interaction with regulatory elements
- **`07.TE_peptide/`** - TE-derived peptide and protein analysis

## Key Features

### Multi-Level TE Analysis
- **TE Loci Level**: Individual transposable element insertion sites
- **TE Subfamily Level**: Grouped analysis by TE family classifications
- **Genomic Context**: Analysis stratified by genomic location (exonic, intronic, intergenic)

### Cross-Pathogen Comparison
- Bacterial pathogens: Mycobacterium tuberculosis, E. coli, etc.
- Viral pathogens: Influenza A virus, SARS-CoV-2, RSV, etc.
- Systematic comparison of TE activation patterns

### Evolutionary Perspective
- Conservation analysis across 95 mammalian species
- LiftOver-based orthologous TE identification
- Evolutionary constraint assessment

## Requirements

### R Dependencies
```r
# Core packages
library(DESeq2)
library(openxlsx)
library(tidyverse)
library(ComplexHeatmap)

# Additional analysis packages
library(enrichR)
library(clusterProfiler)
library(ChIPseeker)
```

### Bioinformatics Tools
- **STAR** - RNA-seq alignment
- **TEcount** - Transposable element quantification
- **featureCounts** - Gene/TE loci quantification
- **fastp** - Quality control
- **bedtools** - Genomic interval operations
- **BWA** - DNA alignment (ATAC-seq)
- **MEME Suite** - Motif discovery
- **MACS2** - Peak calling (ATAC-seq)

### System Requirements
- Linux/Unix environment
- SLURM job scheduler (for HPC scripts)
- Sufficient storage for large genomic datasets
- Multi-core processing capability

## Quick Start

### 1. RNA-seq Processing Pipeline
```bash
# Quality control
bash RNA-seq/02.qc_mapping_quantify/step1.QC_by_fastp.sh

# STAR alignment
bash RNA-seq/02.qc_mapping_quantify/step2.GeneTE_Mapping_by_STAR_parallel.sh

# TE quantification
bash RNA-seq/02.qc_mapping_quantify/step3.Gene_and_TE_Subfamily_Quantification_by_TEcount_MultiMode_parallel.sh
bash RNA-seq/02.qc_mapping_quantify/step4.TE_Loci_Quantification_by_featureCounts_MultiMode.sh
```

### 2. Differential Expression Analysis
```r
# Load core functions
source("02.DESeq2_analysis_Gene/00.Functions.r")

# Run analysis
TE_Gene_analysis(
  metadata_file = "raw_data/sample_info.xlsx",
  raw_count_file = "raw_data/readscounts_matrix.txt",
  outDir = "DE_results/",
  dataset = "Virus_IAV_11",
  cond1 = "Pre_0m",
  cond2 = "Post_360m"
)
```

### 3. TF Motif Analysis
```bash
# Extract TE sequences
bash 05.TF_prediction/01.get_TE_loci_bed_and_sequences.sh

# Find enriched motifs
bash 05.TF_prediction/02.MEME_find_motif.sh
```

## Data Organization

### Input Data Structure
```
00.data/
├── sample_metadata.xlsx          # Sample information and experimental design
├── raw_counts/                   # Raw count matrices
│   ├── gene_counts.txt
│   ├── TE_loci_counts.txt
│   └── TE_subfamily_counts.txt
└── annotations/                  # Genomic annotations
    ├── hg38_TE_annotations.gtf
    └── gencode_annotations.gtf
```

### Sample Naming Convention
- Dataset files: `{Pathogen}_{ExperimentNumber}.txt` (e.g., `Bacteria_1.txt`, `Virus_RSV_1.txt`)
- Sample names: Descriptive identifiers with timepoint information
- Results: Include analysis date and researcher initials (e.g., `hushu_250322`)

## Analysis Workflow

1. **Data Preprocessing**: Quality control and alignment of RNA-seq/ATAC-seq data
2. **Quantification**: Count generation for genes, TE loci, and TE subfamilies
3. **Differential Analysis**: DESeq2-based identification of activated elements
4. **Functional Analysis**: Enrichment analysis and pathway mapping
5. **Evolutionary Analysis**: Conservation assessment across mammals
6. **Regulatory Analysis**: TF binding prediction and regulatory element overlap
7. **Integration**: Multi-omics integration and visualization

## Key Findings

- **Pan-pathogen TE activation**: Different pathogens activate overlapping sets of TE families
- **Immune response correlation**: TE activation correlates with host immune gene expression
- **Evolutionary conservation**: Recurrently activated TEs show specific evolutionary patterns
- **Regulatory potential**: TE-derived sequences contribute transcription factor binding sites

## Citation

If you use this code or data, please cite:

```bibtex
@article{wei2024panpathogen,
  title={Pan-Pathogen Activation of Transposable Elements Shapes Host Immune Responses},
  author={Wei, Xiaoman and Hu, Shu and Xu, Tianqi and Lu, Chengjun and Cui, Jie},
  journal={[Journal Name]},
  year={2024},
  publisher={[Publisher]}
}
```

## Support

For questions about the analysis pipeline or code:
- Create an issue in this repository
- Contact: jiecui@fudan.edu.cn

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Shanghai Sci-Tech Innovation Center for Infection & Immunity
- Fudan University
- All contributors to the data analysis and methodology development
## Updates

- Initial release of comprehensive TE analysis pipeline
