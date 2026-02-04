# ATAC-seq Analysis Pipeline

A streamlined workflow for analyzing ATAC-seq accessibility from Seurat single-cell objects using deepTools.

## Overview

This pipeline extracts ATAC-seq data from Seurat objects, converts them to genome browser formats, and generates publication-ready heatmaps showing chromatin accessibility patterns.

## Workflow Steps

1. **Extract Data from Seurat** → `01_seurat_extraction.sh`
   - Extracts normalized ATAC-seq data from Seurat object
   - Creates bedGraph files for each treatment group
   - Supports both legacy and modern Seurat/Signac versions

2. **Convert to BigWig** → `02_bedgraph_to_bigwig.sh`
   - Converts bedGraph to bigWig format
   - Validates chromosome boundaries
   - Creates log2FC tracks

3. **Compute Matrices** → `03_computematrices.sh`
   - Generates TSS-centered accessibility matrices
   - Optional gene body analysis
   - Custom gene list analysis

4. **Generate Visualizations** → `04_generate_plots.sh`
   - Creates heatmaps and profile plots
   - Multiple visualization styles
   - Publication-ready figures

## Requirements

### R Packages
- Seurat
- Signac
- Matrix
- data.table

### Command-line Tools
- [deepTools](https://deeptools.readthedocs.io/) (computeMatrix, plotHeatmap, plotProfile)
- [UCSC Tools](http://hgdownload.soe.ucsc.edu/admin/exe/) (bedGraphToBigWig)

### Reference Files
- Genome annotation (BED format): TSS coordinates for genes of interest
- Chromosome sizes file (e.g., mm10.chrom.sizes from UCSC)

## Input Files

- **Seurat object**: RDS file containing single-cell ATAC-seq data with normalized counts
- **Gene lists**: CSV files with gene names for focused analysis (optional)
- **Custom regions**: BED files for specific genomic regions (optional)

## Quick Start

```bash
# 1. Set your project directory
export PROJECT_DIR=/path/to/your/project
cd $PROJECT_DIR

# 2. Create required directories
mkdir -p bigwig matrices heatmaps profiles reference

# 3. Edit configuration variables in scripts
# Update file paths in each script to match your data

# 4. Run the workflow
bash 01_seurat_extraction.sh
bash 02_bedgraph_to_bigwig.sh
bash 03_computematrices.sh
bash 04_generate_plots.sh
```

## Customization

### Treatment Groups
Edit the `create_bedgraph()` function calls to match your metadata column names and treatment groups.

### Window Sizes
Default TSS window: ±3kb (adjustable with `-b` and `-a` parameters in computeMatrix)

### Gene Lists
Place your gene list CSV files in the project directory and update paths in scripts.

## Output Files

```
project_dir/
├── bigwig/                    # Browser-viewable tracks
│   ├── *_S3normalized.bw
│   └── *_log2FC.bw
├── matrices/                  # Computed data matrices
│   ├── TSS_accessibility_matrix.gz
│   └── *_matrix.gz
├── heatmaps/                  # Visualization outputs
│   └── *.pdf
└── profiles/                  # Line plot profiles
    └── *.pdf
```

## Notes

- Scripts are designed for mouse genome (mm10), but can be adapted for other species
- Assumes S3/TF-IDF normalized ATAC-seq data in Seurat object
- BigWig files can be loaded into IGV or UCSC Genome Browser

## Citation

If you use this pipeline, please cite:
- [deepTools](https://doi.org/10.1093/nar/gkw257): Ramírez et al., Nucleic Acids Research (2016)
- [Seurat](https://doi.org/10.1016/j.cell.2021.04.048): Hao et al., Cell (2021)
- [Signac](https://doi.org/10.1038/s41592-021-01282-5): Stuart et al., Nature Methods (2021)
