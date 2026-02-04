# Example Data Formats

This file provides examples of the input file formats required for the pipeline.

## Gene Annotation File (BED format)

**File:** `reference/genes.bed`

**Format:** Tab-separated values with 6 columns
```
chromosome    start    end    gene_name    score    strand
```

**Example:**
```
chr1    3204562    3204563    Xkr4     0    -
chr1    4399322    4399323    Lypla1   0    +
chr1    4481008    4481009    Tcea1    0    +
chr1    4763278    4763279    Rgs20    0    -
chr1    4807892    4807893    Klhl17   0    +
chr2    74693581   74693582   Actb     0    -
chr6    125166095  125166096  Gapdh    0    +
chr3    108107280  108107281  Sox2     0    +
chr17   35409847   35409848   Oct4     0    +
chr6    122611448  122611449  Nanog    0    +
```

**Notes:**
- TSS position is used (single base pair for reference-point mode)
- Score column can be 0 or any value
- Strand: `+` for forward, `-` for reverse

### Creating Gene Annotation from GTF/GFF

**Using bedtools:**
```bash
# Extract TSS from GTF
awk '$3=="gene" {
    if($7=="+") 
        print $1"\t"$4-1"\t"$4"\t"$10"\t0\t"$7;
    else 
        print $1"\t"$5-1"\t"$5"\t"$10"\t0\t"$7;
}' genes.gtf | sed 's/"//g' | sed 's/;//g' > genes_tss.bed
```

**From R using GenomicFeatures:**
```R
library(GenomicFeatures)
library(rtracklayer)

# Load GTF
txdb <- makeTxDbFromGFF("genes.gtf")

# Extract genes
genes <- genes(txdb)

# Get TSS
tss <- promoters(genes, upstream=0, downstream=1)

# Export to BED
export.bed(tss, "genes_tss.bed")
```

---

## Gene List CSV Format

**File:** `genes_of_interest.csv`

**Format:** CSV with header, gene names in first column

**Example:**
```csv
gene,log2FoldChange,pvalue
Sox2,2.45,0.001
Nanog,1.89,0.005
Oct4,2.12,0.002
Pax6,-1.67,0.008
Foxg1,1.34,0.015
```

**Notes:**
- Only first column (gene names) is used
- Additional columns are ignored but can be kept for reference
- Gene names must match those in annotation file
- Case-sensitive matching (can be relaxed in script)

---

## Custom Regions BED Format

**File:** `custom_regions.bed`

**Format:** Tab-separated, minimum 3 columns (chr, start, end)

**Basic format (3 columns):**
```
chr1    1000000    1001000
chr1    2000000    2001000
chr2    5000000    5001000
```

**Extended format (6 columns):**
```
chr1    1000000    1001000    peak1    100    +
chr1    2000000    2001000    peak2    200    +
chr2    5000000    5001000    peak3    150    -
```

**Use cases:**
- ChIP-seq peaks (e.g., ERα binding sites)
- ATAC-seq peaks
- Enhancer regions
- Known regulatory elements
- CNVs or structural variants

**Example - ERα binding sites:**
```
chr1    45796588    45796988    ERa_peak_1    856    .
chr1    46123456    46123856    ERa_peak_2    742    .
chr2    12345678    12346078    ERa_peak_3    923    .
```

---

## Chromosome Sizes File

**File:** `reference/mm10.chrom.sizes` or `reference/hg38.chrom.sizes`

**Format:** Tab-separated
```
chromosome    length
```

**Example (mm10):**
```
chr1    195471971
chr2    182113224
chr3    160039680
chr4    156508116
chr5    151834684
chr6    149736546
chr7    145441459
chr8    129401213
chr9    124595110
chr10   130694993
chr11   122082543
chr12   120129022
chr13   120421639
chr14   124902244
chr15   104043685
chr16   98207768
chr17   94987271
chr18   90702639
chr19   61431566
chrX    171031299
chrY    91744698
chrM    16299
```

**Download:**
```bash
# Mouse (mm10)
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# Human (hg38)
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

---

## Seurat Object Requirements

The Seurat object must contain:

1. **ATAC-seq assay** (ChromatinAssay from Signac)
   - Assay name: "ATAC" or "peaks"
   - Normalized data (S3/TF-IDF)

2. **Metadata with treatment groups**
   ```R
   # Check metadata columns
   colnames(seurat_obj@meta.data)
   
   # Should include treatment/condition column:
   # "Treatment", "Condition", "Group", etc.
   ```

3. **Peak naming convention**
   - Format 1: `chr1:1000-2000`
   - Format 2: `chr1_1000_2000`
   - Must be parseable to extract chr, start, end

**Example structure:**
```R
> seurat_obj
An object of class Seurat 
30000 features across 5000 samples within 1 assay 
Active assay: ATAC (30000 features, 0 variable features)

> head(seurat_obj@meta.data)
             nCount_ATAC nFeature_ATAC Treatment
CELL_1           5234        2341       Control
CELL_2           4567        2103       Control
CELL_3           6789        2876       Treatment
CELL_4           5432        2234       Control
CELL_5           7890        3012       Treatment

> rownames(seurat_obj)[1:5]
[1] "chr1:1000-2000"  "chr1:3000-4000"  "chr1:5000-6000"
[4] "chr2:1000-2000"  "chr2:3000-4000"
```

---

## Output File Formats

### bedGraph Format

**File:** `ATAC_Control_normalized.bedGraph`

```
track type=bedGraph name="Control" description="Normalized ATAC-seq"
chr1    1000    2000    0.5432
chr1    3000    4000    1.2345
chr1    5000    6000    0.8765
```

### BigWig Format

Binary format - use IGV or command-line tools to inspect:

```bash
# Get info about bigWig file
bigWigInfo ATAC_Control_normalized.bw

# Extract specific region
bigWigToBedGraph ATAC_Control_normalized.bw -chrom=chr1 -start=1000 -end=2000 stdout
```

### Matrix Format (compressed)

DeepTools matrix format (gzipped):

```bash
# View first few lines
zcat TSS_accessibility_matrix.gz | head -20

# Full decompression (large files!)
gunzip -c TSS_accessibility_matrix.gz > TSS_accessibility_matrix.txt
```

---

## Validation Commands

### Check BED file format:
```bash
# Should have 6 columns, tab-separated
awk -F'\t' '{print NF; exit}' genes.bed
# Output: 6

# Check for proper strand column
awk -F'\t' '{print $6}' genes.bed | sort | uniq
# Output: + and/or -
```

### Check CSV format:
```bash
# Count genes
tail -n +2 genes_of_interest.csv | wc -l

# View first gene
tail -n +2 genes_of_interest.csv | head -1 | cut -d',' -f1
```

### Check chromosome names match:
```bash
# From BED file
cut -f1 genes.bed | sort | uniq

# From bigWig file  
bigWigInfo ATAC_Control_normalized.bw | grep "chr"
```

---

## Common Format Issues

### Issue: Gene names don't match
```bash
# Compare gene names between files
cut -f4 genes.bed | sort > bed_genes.txt
tail -n +2 genes.csv | cut -d',' -f1 | sort > csv_genes.txt
comm -3 bed_genes.txt csv_genes.txt  # Shows differences
```

### Issue: Chromosome naming inconsistency
```bash
# BED has "chr1" but data has "1"
# Fix with sed:
sed 's/^/chr/' genes.bed > genes_fixed.bed

# Or opposite:
sed 's/^chr//' genes.bed > genes_fixed.bed
```

### Issue: Windows line endings (^M characters)
```bash
# Convert to Unix format
dos2unix genes.csv
# or
sed -i 's/\r$//' genes.csv
```
