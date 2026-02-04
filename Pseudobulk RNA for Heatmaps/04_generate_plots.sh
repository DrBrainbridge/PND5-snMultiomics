#!/bin/bash
# 04_generate_plots.sh
# Generate heatmaps and profile plots from computed matrices

set -euo pipefail

echo "=== GENERATING VISUALIZATIONS ==="
echo "Start time: $(date)"

# ============================================================================
# CONFIGURATION - Edit these variables for your analysis
# ============================================================================

# Sample labels
CONTROL_LABEL="Control"
TREATMENT_LABEL="Treatment"

# Plot aesthetics
COLOR_MAP="viridis"  # Options: viridis, plasma, inferno, magma, RdYlBu, etc.
DPI=300              # Resolution for output files
PLOT_FORMAT="pdf"    # Options: pdf, png, svg

# Output directories
HEATMAP_DIR="heatmaps"
PROFILE_DIR="profiles"

# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Create output directories
mkdir -p "$HEATMAP_DIR" "$PROFILE_DIR"

# Check for deepTools
if ! command -v plotHeatmap &> /dev/null; then
    echo "ERROR: plotHeatmap not found"
    echo "Please install deepTools: pip install deeptools"
    exit 1
fi

# ============================================================================
# 1. TSS ACCESSIBILITY PLOTS
# ============================================================================

if [[ -f "matrices/TSS_accessibility_matrix.gz" ]]; then
    echo ""
    echo "Generating TSS accessibility visualizations..."
    
    # Comparison heatmap
    plotHeatmap \
        -m matrices/TSS_accessibility_matrix.gz \
        -o "${HEATMAP_DIR}/TSS_accessibility_comparison.${PLOT_FORMAT}" \
        --plotTitle "ATAC-seq Accessibility at TSS" \
        --zMin 0 0 --zMax auto auto \
        --sortRegions descend \
        --sortUsing mean \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT \
        --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
        --colorMap $COLOR_MAP \
        --xAxisLabel "Distance from TSS (bp)" \
        --yAxisLabel "Genes"
    
    echo "  ✓ Created: TSS_accessibility_comparison.${PLOT_FORMAT}"
    
    # Clustered heatmap
    plotHeatmap \
        -m matrices/TSS_accessibility_matrix.gz \
        -o "${HEATMAP_DIR}/TSS_accessibility_clustered.${PLOT_FORMAT}" \
        --plotTitle "Clustered ATAC-seq Accessibility" \
        --sortRegions descend \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT \
        --colorMap $COLOR_MAP \
        --xAxisLabel "Distance from TSS (bp)"
    
    echo "  ✓ Created: TSS_accessibility_clustered.${PLOT_FORMAT}"
    
    # Profile plot
    plotProfile \
        -m matrices/TSS_accessibility_matrix.gz \
        -o "${PROFILE_DIR}/TSS_accessibility_profile.${PLOT_FORMAT}" \
        --plotTitle "ATAC-seq Accessibility Profile at TSS" \
        --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
        --regionsLabel "All Genes" \
        --plotType lines \
        --perGroup \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT \
        --yAxisLabel "Mean Accessibility"
    
    echo "  ✓ Created: TSS_accessibility_profile.${PLOT_FORMAT}"
else
    echo "✗ TSS matrix not found, skipping TSS plots"
fi

# ============================================================================
# 2. GENE BODY PLOTS (if available)
# ============================================================================

if [[ -f "matrices/genebody_accessibility_matrix.gz" ]]; then
    echo ""
    echo "Generating gene body visualizations..."
    
    plotHeatmap \
        -m matrices/genebody_accessibility_matrix.gz \
        -o "${HEATMAP_DIR}/genebody_accessibility.${PLOT_FORMAT}" \
        --plotTitle "ATAC-seq Accessibility across Gene Bodies" \
        --sortRegions descend \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT \
        --colorMap $COLOR_MAP \
        --xAxisLabel "Gene Body Position"
    
    echo "  ✓ Created: genebody_accessibility.${PLOT_FORMAT}"
    
    plotProfile \
        -m matrices/genebody_accessibility_matrix.gz \
        -o "${PROFILE_DIR}/genebody_accessibility_profile.${PLOT_FORMAT}" \
        --plotTitle "ATAC-seq Profile across Gene Bodies" \
        --plotType lines \
        --perGroup \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT
    
    echo "  ✓ Created: genebody_accessibility_profile.${PLOT_FORMAT}"
fi

# ============================================================================
# 3. DIFFERENTIAL ACCESSIBILITY PLOTS (if available)
# ============================================================================

if [[ -f "matrices/differential_TSS_matrix.gz" ]]; then
    echo ""
    echo "Generating differential accessibility plots..."
    
    plotHeatmap \
        -m matrices/differential_TSS_matrix.gz \
        -o "${HEATMAP_DIR}/differential_accessibility.${PLOT_FORMAT}" \
        --plotTitle "Differential Accessibility (Treatment vs Control)" \
        --zMin -2 --zMax 2 \
        --sortRegions descend \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT \
        --colorMap RdBu_r \
        --xAxisLabel "Distance from TSS (bp)" \
        --yAxisLabel "Genes"
    
    echo "  ✓ Created: differential_accessibility.${PLOT_FORMAT}"
    
    plotProfile \
        -m matrices/differential_TSS_matrix.gz \
        -o "${PROFILE_DIR}/differential_accessibility_profile.${PLOT_FORMAT}" \
        --plotTitle "Differential Accessibility Profile" \
        --plotType lines \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT \
        --yAxisLabel "Log2 Fold Change"
    
    echo "  ✓ Created: differential_accessibility_profile.${PLOT_FORMAT}"
fi

# ============================================================================
# 4. CUSTOM GENE LIST PLOTS (if available)
# ============================================================================

if [[ -f "matrices/specific_genes_TSS_matrix.gz" ]]; then
    echo ""
    echo "Generating specific gene visualizations..."
    
    plotHeatmap \
        -m matrices/specific_genes_TSS_matrix.gz \
        -o "${HEATMAP_DIR}/specific_genes_TSS.${PLOT_FORMAT}" \
        --plotTitle "ATAC-seq at Genes of Interest" \
        --sortRegions descend \
        --sortUsing mean \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT \
        --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
        --colorMap $COLOR_MAP \
        --yAxisLabel "Genes of Interest"
    
    echo "  ✓ Created: specific_genes_TSS.${PLOT_FORMAT}"
    
    plotProfile \
        -m matrices/specific_genes_TSS_matrix.gz \
        -o "${PROFILE_DIR}/specific_genes_profile.${PLOT_FORMAT}" \
        --plotTitle "Accessibility Profile at Genes of Interest" \
        --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
        --plotType lines \
        --perGroup \
        --dpi $DPI \
        --plotFileFormat $PLOT_FORMAT
    
    echo "  ✓ Created: specific_genes_profile.${PLOT_FORMAT}"
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== VISUALIZATION COMPLETE ==="

# Count generated files
heatmap_count=$(find "$HEATMAP_DIR" -type f 2>/dev/null | wc -l)
profile_count=$(find "$PROFILE_DIR" -type f 2>/dev/null | wc -l)

echo "Generated files:"
echo "  Heatmaps: $heatmap_count files in $HEATMAP_DIR/"
echo "  Profiles: $profile_count files in $PROFILE_DIR/"

# List all output files
if [[ $heatmap_count -gt 0 ]]; then
    echo ""
    echo "Heatmap files:"
    ls -lh "$HEATMAP_DIR"/* | awk '{print "  " $9 " (" $5 ")"}'
fi

if [[ $profile_count -gt 0 ]]; then
    echo ""
    echo "Profile files:"
    ls -lh "$PROFILE_DIR"/* | awk '{print "  " $9 " (" $5 ")"}'
fi

echo ""
echo "Visualization completed: $(date)"
echo ""
echo "All plots are publication-ready!"
echo "Files can be opened with standard PDF/image viewers."
