#!/bin/bash
# 03b_custom_genelist_analysis.sh
# Compute matrices for specific genes of interest from CSV files

set -euo pipefail

echo "=== CUSTOM GENE LIST ANALYSIS ==="
echo "Start time: $(date)"

# ============================================================================
# CONFIGURATION - Edit these variables for your analysis
# ============================================================================

# CSV file with genes of interest (first column should be gene names)
GENE_CSV="genes_of_interest.csv"

# Gene annotation file (BED format)
GENES_BED="reference/genes.bed"

# BigWig files
CONTROL_BW="bigwig/ATAC_Control_normalized.bw"
TREATMENT_BW="bigwig/ATAC_Treatment_normalized.bw"

# Analysis parameters
TSS_WINDOW_UPSTREAM=5000
TSS_WINDOW_DOWNSTREAM=5000
BIN_SIZE=5  # High resolution for focused analysis

# Sample labels
CONTROL_LABEL="Control"
TREATMENT_LABEL="Treatment"

# Number of threads
THREADS=8

# ============================================================================
# EXTRACT GENES FROM CSV
# ============================================================================

echo ""
echo "Extracting genes from CSV: $GENE_CSV"

if [[ ! -f "$GENE_CSV" ]]; then
    echo "ERROR: CSV file not found: $GENE_CSV"
    exit 1
fi

# Extract gene names (skip header, get 1st column, remove quotes)
tail -n +2 "$GENE_CSV" | cut -d',' -f1 | sed 's/"//g' > extracted_genes.txt

total_genes=$(wc -l < extracted_genes.txt)
echo "  Extracted $total_genes gene names from CSV"
echo "  First 10 genes:"
head -10 extracted_genes.txt | sed 's/^/    /'

# ============================================================================
# CREATE SPECIFIC GENE BED FILE
# ============================================================================

echo ""
echo "Finding genes in annotation..."

if [[ ! -f "$GENES_BED" ]]; then
    echo "ERROR: Gene annotation not found: $GENES_BED"
    exit 1
fi

# Create BED file with specific genes
> specific_genes.bed
found_count=0

while IFS= read -r gene; do
    # Try exact match first
    if grep -w "^[^\t]*\t[^\t]*\t[^\t]*\t${gene}\t" "$GENES_BED" >> specific_genes.bed 2>/dev/null; then
        echo "  ✓ Found: $gene"
        ((found_count++))
    # Try case-insensitive match
    elif grep -i -w "^[^\t]*\t[^\t]*\t[^\t]*\t${gene}\t" "$GENES_BED" >> specific_genes.bed 2>/dev/null; then
        echo "  ✓ Found (case-insensitive): $gene"
        ((found_count++))
    else
        echo "  ✗ Not found: $gene"
    fi
done < extracted_genes.txt

echo ""
echo "Results: Found $found_count of $total_genes genes in annotation"
echo "Created specific_genes.bed with $(wc -l < specific_genes.bed) entries"

if [[ $found_count -eq 0 ]]; then
    echo "ERROR: No genes found in annotation"
    exit 1
fi

# ============================================================================
# COMPUTE MATRIX FOR SPECIFIC GENES
# ============================================================================

echo ""
echo "Computing high-resolution TSS matrix for specific genes..."
echo "  Window: ±${TSS_WINDOW_UPSTREAM}bp"
echo "  Bin size: ${BIN_SIZE}bp (high resolution)"
echo "  Genes: $found_count"

mkdir -p matrices

computeMatrix reference-point \
    --referencePoint TSS \
    -b $TSS_WINDOW_UPSTREAM -a $TSS_WINDOW_DOWNSTREAM \
    -S "$CONTROL_BW" "$TREATMENT_BW" \
    -R specific_genes.bed \
    --skipZeros \
    --missingDataAsZero \
    --numberOfProcessors $THREADS \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --binSize $BIN_SIZE \
    -o matrices/specific_genes_TSS_matrix.gz

if [[ $? -eq 0 ]]; then
    echo "✓ Matrix created: matrices/specific_genes_TSS_matrix.gz"
    ls -lh matrices/specific_genes_TSS_matrix.gz
else
    echo "ERROR: Matrix computation failed"
    exit 1
fi

# ============================================================================
# OPTIONAL: CREATE HEATMAP
# ============================================================================

echo ""
echo "Generating heatmap for specific genes..."

mkdir -p heatmaps

plotHeatmap \
    -m matrices/specific_genes_TSS_matrix.gz \
    -o heatmaps/specific_genes_TSS_heatmap.pdf \
    --plotTitle "ATAC-seq at Genes of Interest" \
    --sortRegions descend \
    --sortUsing mean \
    --dpi 300 \
    --plotFileFormat pdf \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --yAxisLabel "Genes of Interest" \
    --xAxisLabel "Distance from TSS (bp)" \
    --colorMap viridis

echo "✓ Heatmap created: heatmaps/specific_genes_TSS_heatmap.pdf"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== ANALYSIS COMPLETE ==="
echo "Input genes: $total_genes"
echo "Found in annotation: $found_count"
echo "Output files:"
echo "  - specific_genes.bed"
echo "  - matrices/specific_genes_TSS_matrix.gz"
echo "  - heatmaps/specific_genes_TSS_heatmap.pdf"
echo ""
echo "Completed: $(date)"

# Clean up temporary file
rm extracted_genes.txt
