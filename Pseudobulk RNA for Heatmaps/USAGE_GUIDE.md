# ATAC-seq Analysis Pipeline - Usage Guide

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Quick Start](#quick-start)
3. [Detailed Workflow](#detailed-workflow)
4. [Script Descriptions](#script-descriptions)
5. [Customization Guide](#customization-guide)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Usage](#advanced-usage)

---

## Prerequisites

### Required Software

#### R Packages
```R
install.packages("data.table")
install.packages("BiocManager")
BiocManager::install("Seurat")
BiocManager::install("Signac")
BiocManager::install("Matrix")
```

#### Python Tools (deepTools)
```bash
pip install deeptools
# or
conda install -c bioconda deeptools
```

#### UCSC Tools
```bash
# Download bedGraphToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod +x bedGraphToBigWig
# Move to a directory in your PATH
```

### Required Input Files

1. **Seurat Object** (`.rds` file)
   - Must contain single-cell ATAC-seq data
   - Should have normalized accessibility values (S3/TF-IDF)
   - Metadata column with treatment/condition labels

2. **Gene Annotation** (BED format)
   ```
   chr1    3204562    3204563    Xkr4    0    -
   chr1    4399322    4399323    Lypla1  0    +
   ```
   Format: `chromosome    start    end    gene_name    score    strand`

3. **Chromosome Sizes File**
   - Downloaded automatically by script
   - Or manually from: `http://hgdownload.cse.ucsc.edu/goldenPath/[genome]/bigZips/`

---

## Quick Start

```bash
# 1. Clone or download the scripts
git clone [your-repo-url]
cd atac-seq-pipeline

# 2. Make scripts executable
chmod +x *.sh

# 3. Prepare your data
mkdir -p bigwig matrices heatmaps profiles reference

# 4. Edit configuration in each script
# Update file paths and parameters in the CONFIGURATION section

# 5. Run the pipeline
./01_seurat_extraction.sh
./02_bedgraph_to_bigwig.sh
./03_computematrices.sh
./04_generate_plots.sh
```

---

## Detailed Workflow

### Step 1: Extract Data from Seurat Object

**Script:** `01_seurat_extraction.sh`

**What it does:**
- Loads Seurat object containing scATAC-seq data
- Extracts normalized accessibility values
- Groups cells by treatment condition
- Calculates mean accessibility per peak
- Creates bedGraph files for each condition
- Optionally creates log2FC bedGraph

**Key parameters to edit:**
```bash
SEURAT_OBJECT="your_seurat_object.rds"
METADATA_COLUMN="Treatment"
CONTROL_GROUP="Control"
TREATMENT_GROUP="Treatment"
```

**Output files:**
- `ATAC_Control_normalized.bedGraph`
- `ATAC_Treatment_normalized.bedGraph`
- `ATAC_log2FC_Treatment_vs_Control.bedGraph`

### Step 2: Convert to BigWig Format

**Script:** `02_bedgraph_to_bigwig.sh`

**What it does:**
- Sorts bedGraph files by genomic coordinates
- Validates and fixes chromosome boundaries
- Converts to bigWig format (compact, indexed)
- Creates browser-viewable tracks

**Key parameters:**
```bash
GENOME="mm10"  # or hg38, etc.
```

**Output files:**
- `bigwig/ATAC_Control_normalized.bw`
- `bigwig/ATAC_Treatment_normalized.bw`
- `bigwig/ATAC_log2FC_Treatment_vs_Control.bw`

### Step 3: Compute Accessibility Matrices

**Script:** `03_computematrices.sh`

**What it does:**
- Centers analysis on gene TSS regions
- Computes accessibility in windows around TSS
- Creates data matrices for visualization
- Optionally analyzes gene bodies

**Key parameters:**
```bash
TSS_WINDOW_UPSTREAM=3000    # bp upstream of TSS
TSS_WINDOW_DOWNSTREAM=3000  # bp downstream
BIN_SIZE=10                 # resolution
```

**Output files:**
- `matrices/TSS_accessibility_matrix.gz`
- `matrices/genebody_accessibility_matrix.gz` (optional)
- `matrices/differential_TSS_matrix.gz` (optional)

### Step 4: Generate Visualizations

**Script:** `04_generate_plots.sh`

**What it does:**
- Creates heatmaps showing accessibility patterns
- Generates profile plots (average signal)
- Multiple visualization styles
- Publication-ready figures

**Key parameters:**
```bash
COLOR_MAP="viridis"  # or plasma, magma, RdBu, etc.
DPI=300
PLOT_FORMAT="pdf"    # or png, svg
```

**Output files:**
- `heatmaps/*.pdf` - Heatmap visualizations
- `profiles/*.pdf` - Profile plots

---

## Script Descriptions

### Main Pipeline Scripts

| Script | Purpose | Runtime | Memory |
|--------|---------|---------|---------|
| `01_seurat_extraction.sh` | Extract data from Seurat | ~10-30 min | 32-64 GB |
| `02_bedgraph_to_bigwig.sh` | Format conversion | ~5-15 min | 16-32 GB |
| `03_computematrices.sh` | Compute data matrices | ~30-60 min | 16-32 GB |
| `04_generate_plots.sh` | Create visualizations | ~10-20 min | 8-16 GB |

### Optional Analysis Scripts

| Script | Purpose | Use Case |
|--------|---------|----------|
| `03b_custom_genelist_analysis.sh` | Specific genes | Focus on genes of interest |
| `03c_custom_regions_analysis.sh` | Custom genomic regions | Analyze binding sites, peaks |
| `03d_dual_genelist_comparison.sh` | Compare two gene lists | Up vs down-regulated DEGs |

---

## Customization Guide

### Adjusting Analysis Windows

**For promoter-focused analysis:**
```bash
TSS_WINDOW_UPSTREAM=1000
TSS_WINDOW_DOWNSTREAM=500
```

**For distal regulatory element analysis:**
```bash
TSS_WINDOW_UPSTREAM=10000
TSS_WINDOW_DOWNSTREAM=10000
```

### Changing Resolution

**High resolution (slower, more detail):**
```bash
BIN_SIZE=1  # 1 bp bins
```

**Standard resolution:**
```bash
BIN_SIZE=10  # 10 bp bins
```

**Low resolution (faster, less detail):**
```bash
BIN_SIZE=50  # 50 bp bins
```

### Color Schemes for Heatmaps

```bash
# Sequential (for accessibility)
COLOR_MAP="viridis"    # purple to yellow
COLOR_MAP="plasma"     # purple to pink
COLOR_MAP="Blues"      # white to blue

# Diverging (for differential analysis)
COLOR_MAP="RdBu_r"     # red to blue (reversed)
COLOR_MAP="RdYlBu"     # red-yellow-blue
```

### Adapting for Different Organisms

**Human (hg38):**
```bash
GENOME="hg38"
standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
```

**Mouse (mm10):**
```bash
GENOME="mm10"
standard_chrs <- paste0("chr", c(1:19, "X", "Y", "M"))
```

**Other organisms:**
Update chromosome names and counts accordingly.

---

## Troubleshooting

### Common Issues

**1. "Seurat object not found"**
```bash
# Check file path and existence
ls -lh your_seurat_object.rds
```

**2. "No genes found in annotation"**
- Verify gene names match between CSV and BED files
- Check for case sensitivity issues
- Ensure BED file is properly formatted

**3. "computeMatrix failed"**
- Check that bigWig files are valid: `file bigwig/*.bw`
- Ensure BED file has no errors
- Try reducing number of threads if memory issues

**4. "Chromosome boundaries exceeded"**
- This is handled automatically by scripts
- If persists, check your genome build matches data

**5. Memory errors**
```bash
# Reduce memory usage by:
# - Processing fewer genes at once
# - Increasing bin size (lower resolution)
# - Using fewer threads
```

### Validation Checks

**Check bedGraph format:**
```bash
head -20 ATAC_Control_normalized.bedGraph
# Should show: chr start end score
```

**Check bigWig files:**
```bash
# Count entries
bigWigInfo bigwig/ATAC_Control_normalized.bw
```

**Check matrix files:**
```bash
# Extract first few rows
zcat matrices/TSS_accessibility_matrix.gz | head
```

---

## Advanced Usage

### Analyzing Multiple Conditions

Modify `01_seurat_extraction.sh` to loop over multiple conditions:

```bash
CONDITIONS=("Control" "Treatment1" "Treatment2" "Treatment3")

for condition in "${CONDITIONS[@]}"; do
    create_bedgraph "$condition" "$OUTPUT_PREFIX"
done
```

### Batch Processing Multiple Samples

Create a wrapper script:

```bash
#!/bin/bash
SAMPLES=("sample1" "sample2" "sample3")

for sample in "${SAMPLES[@]}"; do
    export SEURAT_OBJECT="${sample}_seurat.rds"
    export OUTPUT_PREFIX="${sample}_ATAC"
    
    ./01_seurat_extraction.sh
    ./02_bedgraph_to_bigwig.sh
    ./03_computematrices.sh
    ./04_generate_plots.sh
done
```

### Creating Custom Gene Lists

**From R/Seurat:**
```R
# Get differentially expressed genes
markers <- FindMarkers(seurat_obj, ident.1 = "condition1", ident.2 = "condition2")
deg_genes <- rownames(markers[markers$p_val_adj < 0.05, ])
write.csv(deg_genes, "deg_genes.csv", row.names = FALSE)
```

**From Python/Scanpy:**
```python
# Get marker genes
sc.tl.rank_genes_groups(adata, groupby='condition')
deg_df = sc.get.rank_genes_groups_df(adata, group='condition1')
deg_df[['names']].to_csv('deg_genes.csv', index=False)
```

### Integrating with Existing Workflows

**Load in R for further analysis:**
```R
library(data.table)
matrix_data <- fread("matrices/TSS_accessibility_matrix.gz")
```

**Use bigWig files in IGV:**
1. Open IGV
2. File → Load from File
3. Select your `.bw` files
4. Browse to genes of interest

---

## Best Practices

1. **Always validate intermediate outputs**
   - Check file sizes are reasonable
   - Inspect first few lines of each output
   - Verify gene counts match expectations

2. **Keep original data separate**
   - Don't overwrite source Seurat object
   - Store outputs in dedicated directories

3. **Document your analysis**
   - Save script configurations
   - Note any modifications made
   - Record software versions

4. **Use version control**
   ```bash
   git init
   git add *.sh README.md
   git commit -m "Initial pipeline setup"
   ```

5. **Test with subset first**
   - Try with small gene list initially
   - Verify pipeline works end-to-end
   - Then scale to full analysis

---

## Citation

If you use this pipeline in your research, please cite:

- **deepTools**: Ramírez et al., Nucleic Acids Research (2016)
- **Seurat**: Hao et al., Cell (2021)
- **Signac**: Stuart et al., Nature Methods (2021)

---

## Support

For questions or issues:
1. Check the [Troubleshooting](#troubleshooting) section
2. Review deepTools documentation: https://deeptools.readthedocs.io/
3. Open an issue on GitHub (if applicable)

---

## License

These scripts are provided as-is for academic and research use.
