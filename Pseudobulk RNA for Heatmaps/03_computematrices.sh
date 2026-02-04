#!/bin/bash
# 03_computematrices.sh
# Compute accessibility matrices around gene TSS using deepTools

set -euo pipefail

echo "=== TSS MATRIX COMPUTATION ==="
echo "Start time: $(date)"

# ============================================================================
# CONFIGURATION - Edit these variables for your analysis
# ============================================================================

# BigWig files from previous step
CONTROL_BW="bigwig/ATAC_Control_normalized.bw"
TREATMENT_BW="bigwig/ATAC_Treatment_normalized.bw"
LOG2FC_BW="bigwig/ATAC_log2FC_Treatment_vs_Control.bw"

# Gene annotation file (BED format: chr start end name score strand)
GENES_BED="reference/genes.bed"

# TSS analysis parameters
TSS_WINDOW_UPSTREAM=3000    # Base pairs upstream of TSS
TSS_WINDOW_DOWNSTREAM=3000  # Base pairs downstream of TSS
BIN_SIZE=10                 # Base pairs per bin

# Gene body analysis parameters (optional)
GENEBODY_SIZE=5000          # Scale gene bodies to this size
GENEBODY_UPSTREAM=2000
GENEBODY_DOWNSTREAM=2000

# Computational resources
THREADS=8  # Number of CPU threads to use

# Sample labels for plots
CONTROL_LABEL="Control"
TREATMENT_LABEL="Treatment"

# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Create output directory
mkdir -p matrices

# Validate input files
echo "Checking input files..."
for file in "$CONTROL_BW" "$TREATMENT_BW"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file missing: $file"
        exit 1
    fi
    echo "✓ Found: $file ($(du -h "$file" | cut -f1))"
done

# Check for gene annotation
if [[ ! -f "$GENES_BED" ]]; then
    echo ""
    echo "WARNING: Gene annotation file not found: $GENES_BED"
    echo "Creating example file with common genes..."
    
    mkdir -p reference
    cat > "$GENES_BED" << 'EOF'
chr1	3204562	3204563	Xkr4	0	-
chr1	4399322	4399323	Lypla1	0	+
chr1	4481008	4481009	Tcea1	0	+
chr1	4763278	4763279	Rgs20	0	-
chr2	74693581	74693582	Actb	0	-
chr6	125166095	125166096	Gapdh	0	+
chr3	108107280	108107281	Sox2	0	+
EOF
    echo "✓ Created example annotation with $(wc -l < "$GENES_BED") genes"
    echo "  Please replace with your complete gene annotation"
fi

bed_lines=$(wc -l < "$GENES_BED")
echo "✓ Using $bed_lines genes from: $GENES_BED"

# Check for deepTools
if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: computeMatrix not found"
    echo "Please install deepTools: pip install deeptools"
    exit 1
fi

# ============================================================================
# 1. TSS-CENTERED ANALYSIS
# ============================================================================

echo ""
echo "Computing TSS-centered matrix..."
echo "  Window: ±${TSS_WINDOW_UPSTREAM}bp around TSS"
echo "  Bin size: ${BIN_SIZE}bp"
echo "  Genes: $bed_lines"

computeMatrix reference-point \
    --referencePoint TSS \
    -b $TSS_WINDOW_UPSTREAM -a $TSS_WINDOW_DOWNSTREAM \
    -S "$CONTROL_BW" "$TREATMENT_BW" \
    -R "$GENES_BED" \
    --skipZeros \
    --missingDataAsZero \
    --numberOfProcessors $THREADS \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --binSize $BIN_SIZE \
    -o matrices/TSS_accessibility_matrix.gz \
    --outFileSortedRegions matrices/TSS_sorted_regions.bed

if [[ $? -eq 0 ]]; then
    echo "✓ TSS matrix created: matrices/TSS_accessibility_matrix.gz"
else
    echo "ERROR: TSS matrix computation failed"
    exit 1
fi

# ============================================================================
# 2. GENE BODY ANALYSIS (optional)
# ============================================================================

echo ""
echo "Computing gene body matrix..."

if computeMatrix scale-regions \
    -S "$CONTROL_BW" "$TREATMENT_BW" \
    -R "$GENES_BED" \
    -m $GENEBODY_SIZE \
    -b $GENEBODY_UPSTREAM -a $GENEBODY_DOWNSTREAM \
    --skipZeros \
    --missingDataAsZero \
    --numberOfProcessors $THREADS \
    --samplesLabel "$CONTROL_LABEL" "$TREATMENT_LABEL" \
    --binSize $BIN_SIZE \
    -o matrices/genebody_accessibility_matrix.gz; then
    echo "✓ Gene body matrix created"
else
    echo "✗ Gene body analysis failed (optional)"
fi

# ============================================================================
# 3. DIFFERENTIAL ACCESSIBILITY (optional)
# ============================================================================

if [[ -f "$LOG2FC_BW" ]]; then
    echo ""
    echo "Computing differential accessibility matrix..."
    
    if computeMatrix reference-point \
        --referencePoint TSS \
        -b $TSS_WINDOW_UPSTREAM -a $TSS_WINDOW_DOWNSTREAM \
        -S "$LOG2FC_BW" \
        -R "$GENES_BED" \
        --missingDataAsZero \
        --numberOfProcessors $THREADS \
        --samplesLabel "Log2FC" \
        --binSize $BIN_SIZE \
        -o matrices/differential_TSS_matrix.gz; then
        echo "✓ Differential matrix created"
    else
        echo "✗ Differential analysis failed (optional)"
    fi
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== RESULTS SUMMARY ==="
for matrix in matrices/*.gz; do
    if [[ -f "$matrix" ]]; then
        echo "✓ $(basename "$matrix") ($(du -h "$matrix" | cut -f1))"
    fi
done

echo ""
echo "Sorted regions file:"
if [[ -f "matrices/TSS_sorted_regions.bed" ]]; then
    echo "✓ matrices/TSS_sorted_regions.bed ($(wc -l < matrices/TSS_sorted_regions.bed) regions)"
fi

echo ""
echo "Matrix computation completed: $(date)"
echo "Next step: Run 04_generate_plots.sh"
