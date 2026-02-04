#!/bin/bash
# 03d_dual_genelist_comparison.sh
# Compare ATAC-seq accessibility at two different gene lists
# (e.g., upregulated vs downregulated DEGs)

set -euo pipefail

echo "=== DUAL GENE LIST COMPARISON ==="
echo "Start time: $(date)"

# ============================================================================
# CONFIGURATION - Edit these variables for your analysis
# ============================================================================

# CSV files with two gene lists
# First column should contain gene names
GENELIST1_CSV="upregulated_genes.csv"
GENELIST2_CSV="downregulated_genes.csv"

# Labels for the two gene groups
GROUP1_LABEL="Upregulated DEGs"
GROUP2_LABEL="Downregulated DEGs"

# Gene annotation file (BED format)
GENES_BED="reference/genes.bed"

# BigWig files
CONTROL_BW="bigwig/ATAC_Control_normalized.bw"
TREATMENT_BW="bigwig/ATAC_Treatment_normalized.bw"

# Analysis parameters
TSS_WINDOW_UPSTREAM=5000
TSS_WINDOW_DOWNSTREAM=5000
BIN_SIZE=5  # High resolution

# Sample labels
CONTROL_LABEL="Control"
TREATMENT_LABEL="Treatment"

# Number of threads
THREADS=8

# ============================================================================
# FUNCTION: Extract genes from CSV and create BED file
# ============================================================================

extract_genes_to_bed() {
    local csv_file=$1
    local output_bed=$2
    local group_label=$3
    
    echo ""
    echo "Processing: $group_label"
    
    if [[ ! -f "$csv_file" ]]; then
        echo "ERROR: CSV file not found: $csv_file"
        return 1
    fi
    
    # Extract gene names (skip header, get 1st column, remove quotes)
    local temp_genes="temp_${group_label// /_}_genes.txt"
    tail -n +2 "$csv_file" | cut -d',' -f1 | sed 's/"//g' > "$temp_genes"
    
    local total_genes=$(wc -l < "$temp_genes")
    echo "  Extracted $total_genes genes from CSV"
    
    # Find genes in annotation
    > "$output_bed"
    local found_count=0
    
    while IFS= read -r gene; do
        # Try exact match
        if grep -w "^[^\t]*\t[^\t]*\t[^\t]*\t${gene}\t" "$GENES_BED" >> "$output_bed" 2>/dev/null; then
            ((found_count++))
        # Try case-insensitive match
        elif grep -i -w "^[^\t]*\t[^\t]*\t[^\t]*\t${gene}\t" "$GENES_BED" >> "$output_bed" 2>/dev/null; then
            ((found_count++))
        fi
    done < "$temp_genes"
    
    echo "  Found $found_count of $total_genes genes in annotation"
    
    # Clean up
    rm "$temp_genes"
    
    if [[ $found_count -eq 0 ]]; then
        echo "ERROR: No genes found for $group_label"
        return 1
    fi
    
    return 0
}

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

echo "Configuration:"
echo "  Group 1: $GROUP1_LABEL ($GENELIST1_CSV)"
echo "  Group 2: $GROUP2_LABEL ($GENELIST2_CSV)"
echo "  Window: ±${TSS_WINDOW_UPSTREAM}bp around TSS"
echo "  Bin size: ${BIN_SIZE}bp"

if [[ ! -f "$GENES_BED" ]]; then
    echo "ERROR: Gene annotation not found: $GENES_BED"
    exit 1
fi

# ============================================================================
# EXTRACT GENES FROM BOTH CSV FILES
# ============================================================================

mkdir -p reference

GROUP1_BED="reference/group1_genes.bed"
GROUP2_BED="reference/group2_genes.bed"

if ! extract_genes_to_bed "$GENELIST1_CSV" "$GROUP1_BED" "$GROUP1_LABEL"; then
    exit 1
fi

if ! extract_genes_to_bed "$GENELIST2_CSV" "$GROUP2_BED" "$GROUP2_LABEL"; then
    exit 1
fi

group1_count=$(wc -l < "$GROUP1_BED")
group2_count=$(wc -l < "$GROUP2_BED")

echo ""
echo "Gene lists prepared:"
echo "  $GROUP1_LABEL: $group1_count genes"
echo "  $GROUP2_LABEL: $group2_count genes"

# ============================================================================
# COMPUTE DUAL MATRIX
# ============================================================================

echo ""
echo "Computing dual gene list matrix..."
echo "This will create a heatmap with two distinct gene groups"

mkdir -p matrices

computeMatrix reference-point \
    --referencePoint TSS \
    -b $TSS_WINDOW_UPSTREAM -a $TSS_WINDOW_DOWNSTREAM \
    -S "$CONTROL_BW" "$TREATMENT_BW" \
    -R "$GROUP1_BED" "$GROUP2_BED" \
    --skipZeros \
    --missingDataAsZero \
    --numberOfProcessors $THREADS \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --regionsLabel "$GROUP1_LABEL" "$GROUP2_LABEL" \
    --binSize $BIN_SIZE \
    -o matrices/dual_genelist_matrix.gz \
    --outFileSortedRegions matrices/dual_genelist_sorted.bed

if [[ $? -ne 0 ]]; then
    echo "First attempt failed, trying without region labels..."
    computeMatrix reference-point \
        --referencePoint TSS \
        -b $TSS_WINDOW_UPSTREAM -a $TSS_WINDOW_DOWNSTREAM \
        -S "$CONTROL_BW" "$TREATMENT_BW" \
        -R "$GROUP1_BED" "$GROUP2_BED" \
        --skipZeros \
        --missingDataAsZero \
        --numberOfProcessors $THREADS \
        --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
        --binSize $BIN_SIZE \
        -o matrices/dual_genelist_matrix.gz \
        --outFileSortedRegions matrices/dual_genelist_sorted.bed
fi

if [[ $? -eq 0 ]]; then
    echo "✓ Dual matrix created: matrices/dual_genelist_matrix.gz"
else
    echo "ERROR: Matrix computation failed"
    exit 1
fi

# ============================================================================
# GENERATE COMPARATIVE HEATMAP
# ============================================================================

echo ""
echo "Generating comparative heatmap..."

mkdir -p heatmaps

plotHeatmap \
    -m matrices/dual_genelist_matrix.gz \
    -o heatmaps/dual_genelist_comparison.pdf \
    --plotTitle "ATAC-seq Accessibility: Gene List Comparison" \
    --sortRegions descend \
    --sortUsing mean \
    --dpi 300 \
    --plotFileFormat pdf \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --regionsLabel "$GROUP1_LABEL" "$GROUP2_LABEL" \
    --xAxisLabel "Distance from TSS (bp)" \
    --colorMap viridis

echo "✓ Heatmap created: heatmaps/dual_genelist_comparison.pdf"

# Profile plot
mkdir -p profiles

plotProfile \
    -m matrices/dual_genelist_matrix.gz \
    -o profiles/dual_genelist_profile.pdf \
    --plotTitle "Accessibility Comparison: Two Gene Lists" \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --regionsLabel "$GROUP1_LABEL" "$GROUP2_LABEL" \
    --plotType lines \
    --perGroup \
    --dpi 300 \
    --plotFileFormat pdf \
    --yAxisLabel "Mean Accessibility"

echo "✓ Profile created: profiles/dual_genelist_profile.pdf"

# ============================================================================
# OPTIONAL: CREATE INDIVIDUAL MATRICES
# ============================================================================

echo ""
echo "Creating individual matrices for each group..."

# Group 1 matrix
computeMatrix reference-point \
    --referencePoint TSS \
    -b $TSS_WINDOW_UPSTREAM -a $TSS_WINDOW_DOWNSTREAM \
    -S "$CONTROL_BW" "$TREATMENT_BW" \
    -R "$GROUP1_BED" \
    --skipZeros \
    --missingDataAsZero \
    --numberOfProcessors $THREADS \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --binSize $BIN_SIZE \
    -o "matrices/group1_matrix.gz" && \
    echo "✓ Group 1 matrix created"

# Group 2 matrix
computeMatrix reference-point \
    --referencePoint TSS \
    -b $TSS_WINDOW_UPSTREAM -a $TSS_WINDOW_DOWNSTREAM \
    -S "$CONTROL_BW" "$TREATMENT_BW" \
    -R "$GROUP2_BED" \
    --skipZeros \
    --missingDataAsZero \
    --numberOfProcessors $THREADS \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --binSize $BIN_SIZE \
    -o "matrices/group2_matrix.gz" && \
    echo "✓ Group 2 matrix created"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== DUAL GENE LIST ANALYSIS COMPLETE ==="
echo ""
echo "Input:"
echo "  $GROUP1_LABEL: $group1_count genes"
echo "  $GROUP2_LABEL: $group2_count genes"
echo ""
echo "Output files:"
echo "  - matrices/dual_genelist_matrix.gz (combined)"
echo "  - matrices/group1_matrix.gz (group 1 only)"
echo "  - matrices/group2_matrix.gz (group 2 only)"
echo "  - heatmaps/dual_genelist_comparison.pdf"
echo "  - profiles/dual_genelist_profile.pdf"
echo ""
echo "The dual heatmap shows:"
echo "  - Top section: $GROUP1_LABEL"
echo "  - Bottom section: $GROUP2_LABEL"
echo "  - Side-by-side comparison of $CONTROL_LABEL vs $TREATMENT_LABEL"
echo ""
echo "Completed: $(date)"
