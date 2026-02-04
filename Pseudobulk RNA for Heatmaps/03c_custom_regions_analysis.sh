#!/bin/bash
# 03c_custom_regions_analysis.sh
# Compute matrices for custom genomic regions (e.g., binding sites, peaks)

set -euo pipefail

echo "=== CUSTOM GENOMIC REGIONS ANALYSIS ==="
echo "Start time: $(date)"

# ============================================================================
# CONFIGURATION - Edit these variables for your analysis
# ============================================================================

# BED file with custom regions of interest
# Format: chr start end [name] [score] [strand]
REGIONS_BED="custom_regions.bed"

# BigWig files
CONTROL_BW="bigwig/ATAC_Control_normalized.bw"
TREATMENT_BW="bigwig/ATAC_Treatment_normalized.bw"

# Analysis parameters
WINDOW_UPSTREAM=3000      # Base pairs upstream of region center
WINDOW_DOWNSTREAM=3000    # Base pairs downstream of region center
BIN_SIZE=10               # Base pairs per bin

# Reference point (options: center, TSS, TES)
REFERENCE_POINT="center"

# Sample labels
CONTROL_LABEL="Control"
TREATMENT_LABEL="Treatment"

# Region label (for plots)
REGION_LABEL="Custom Regions"

# Number of threads
THREADS=8

# Output names
OUTPUT_PREFIX="custom_regions"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

echo ""
echo "Configuration:"
echo "  Regions file: $REGIONS_BED"
echo "  Control BigWig: $CONTROL_BW"
echo "  Treatment BigWig: $TREATMENT_BW"
echo "  Window: ±${WINDOW_UPSTREAM}bp around $REFERENCE_POINT"
echo "  Bin size: ${BIN_SIZE}bp"

# Check files exist
for file in "$REGIONS_BED" "$CONTROL_BW" "$TREATMENT_BW"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Report region count
region_count=$(wc -l < "$REGIONS_BED")
echo ""
echo "Found $region_count regions in BED file"
echo "First 5 regions:"
head -5 "$REGIONS_BED" | sed 's/^/  /'

# ============================================================================
# COMPUTE MATRIX
# ============================================================================

echo ""
echo "Computing accessibility matrix for custom regions..."

mkdir -p matrices

computeMatrix reference-point \
    --referencePoint $REFERENCE_POINT \
    -b $WINDOW_UPSTREAM -a $WINDOW_DOWNSTREAM \
    -S "$CONTROL_BW" "$TREATMENT_BW" \
    -R "$REGIONS_BED" \
    --skipZeros \
    --missingDataAsZero \
    --numberOfProcessors $THREADS \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --binSize $BIN_SIZE \
    -o "matrices/${OUTPUT_PREFIX}_matrix.gz"

if [[ $? -eq 0 ]]; then
    echo "✓ Matrix created: matrices/${OUTPUT_PREFIX}_matrix.gz"
    ls -lh "matrices/${OUTPUT_PREFIX}_matrix.gz"
else
    echo "ERROR: Matrix computation failed"
    exit 1
fi

# ============================================================================
# GENERATE VISUALIZATIONS
# ============================================================================

echo ""
echo "Generating visualizations..."

mkdir -p heatmaps profiles

# Heatmap
plotHeatmap \
    -m "matrices/${OUTPUT_PREFIX}_matrix.gz" \
    -o "heatmaps/${OUTPUT_PREFIX}_heatmap.pdf" \
    --plotTitle "ATAC-seq Accessibility at ${REGION_LABEL}" \
    --sortRegions descend \
    --sortUsing mean \
    --dpi 300 \
    --plotFileFormat pdf \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --regionsLabel "$REGION_LABEL" \
    --xAxisLabel "Distance from ${REFERENCE_POINT} (bp)" \
    --yAxisLabel "$REGION_LABEL" \
    --colorMap viridis

echo "✓ Heatmap created: heatmaps/${OUTPUT_PREFIX}_heatmap.pdf"

# Profile plot
plotProfile \
    -m "matrices/${OUTPUT_PREFIX}_matrix.gz" \
    -o "profiles/${OUTPUT_PREFIX}_profile.pdf" \
    --plotTitle "Average Accessibility at ${REGION_LABEL}" \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --plotType lines \
    --perGroup \
    --dpi 300 \
    --plotFileFormat pdf \
    --yAxisLabel "Mean Accessibility" \
    --xAxisLabel "Distance from ${REFERENCE_POINT} (bp)"

echo "✓ Profile created: profiles/${OUTPUT_PREFIX}_profile.pdf"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== ANALYSIS COMPLETE ==="
echo "Regions analyzed: $region_count"
echo "Output files:"
echo "  - matrices/${OUTPUT_PREFIX}_matrix.gz"
echo "  - heatmaps/${OUTPUT_PREFIX}_heatmap.pdf"
echo "  - profiles/${OUTPUT_PREFIX}_profile.pdf"
echo ""
echo "Completed: $(date)"

# ============================================================================
# EXAMPLE: Analyzing ERα binding sites
# ============================================================================
# To analyze ERα binding sites or other ChIP-seq peaks:
# 1. Prepare a BED file with your peaks:
#    chr1  1000  2000  peak1
#    chr1  3000  4000  peak2
# 2. Update REGIONS_BED variable to point to your BED file
# 3. Set meaningful REGION_LABEL (e.g., "ERα Binding Sites")
# 4. Run this script
