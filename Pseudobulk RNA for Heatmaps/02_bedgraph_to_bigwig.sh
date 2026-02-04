#!/bin/bash
# 02_bedgraph_to_bigwig.sh
# Convert bedGraph files to bigWig format for genome browser visualization

set -euo pipefail

echo "Starting bedGraph to bigWig conversion"
echo "Start time: $(date)"

# ============================================================================
# CONFIGURATION - Edit these variables for your analysis
# ============================================================================

# Genome and reference files
GENOME="mm10"  # Change to your genome (hg38, mm10, etc.)
CHROM_SIZES="reference/${GENOME}.chrom.sizes"

# Input file pattern (matches output from step 01)
BEDGRAPH_PATTERN="ATAC_*_normalized.bedGraph"

# Output directory for bigWig files
BIGWIG_DIR="bigwig"

# Path to UCSC bedGraphToBigWig tool
# Download from: http://hgdownload.soe.ucsc.edu/admin/exe/
BEDGRAPHTOBIGWIG="bedGraphToBigWig"  # Assumes it's in PATH

# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Create output directory
mkdir -p "$BIGWIG_DIR"
mkdir -p reference

# Download chromosome sizes if not present
if [[ ! -f "$CHROM_SIZES" ]]; then
    echo "Downloading ${GENOME} chromosome sizes..."
    wget -q -O "$CHROM_SIZES" \
        "http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZigs/mm10.chrom.sizes"
    echo "✓ Downloaded chromosome sizes"
fi

# Check for bedGraphToBigWig tool
if ! command -v $BEDGRAPHTOBIGWIG &> /dev/null; then
    echo "ERROR: bedGraphToBigWig not found in PATH"
    echo "Please download from: http://hgdownload.soe.ucsc.edu/admin/exe/"
    exit 1
fi

# Function to convert bedGraph to bigWig
convert_to_bigwig() {
    local bedgraph_file=$1
    local basename=$(basename "$bedgraph_file" .bedGraph)
    local bigwig_file="${BIGWIG_DIR}/${basename}.bw"
    
    echo ""
    echo "Processing: $basename"
    
    # Skip track line, sort, and prepare file
    grep -v "^track" "$bedgraph_file" | \
        sort -k1,1 -k2,2n > "${bedgraph_file}.sorted"
    
    # Fix coordinates that exceed chromosome boundaries
    awk -v OFS='\t' '
        NR==FNR {
            chrom_size[$1]=$2
            next
        }
        {
            if ($1 in chrom_size) {
                # Fix end coordinate if it exceeds chromosome size
                if ($3 > chrom_size[$1]) {
                    $3 = chrom_size[$1]
                }
                # Only output valid intervals
                if ($2 < $3) {
                    print
                }
            }
        }
    ' "$CHROM_SIZES" "${bedgraph_file}.sorted" > "${bedgraph_file}.fixed"
    
    # Report filtering
    original_lines=$(wc -l < "${bedgraph_file}.sorted")
    fixed_lines=$(wc -l < "${bedgraph_file}.fixed")
    echo "  Lines: $original_lines → $fixed_lines (filtered: $((original_lines - fixed_lines)))"
    
    # Convert to bigWig
    $BEDGRAPHTOBIGWIG "${bedgraph_file}.fixed" "$CHROM_SIZES" "$bigwig_file"
    
    # Validate output
    if [[ ! -f "$bigwig_file" ]] || [[ ! -s "$bigwig_file" ]]; then
        echo "  ERROR: bigWig conversion failed"
        return 1
    fi
    
    echo "  ✓ Created: $bigwig_file ($(du -h "$bigwig_file" | cut -f1))"
    
    # Clean up temporary files
    rm "${bedgraph_file}.sorted" "${bedgraph_file}.fixed"
    
    return 0
}

# Process all matching bedGraph files
echo ""
echo "Finding bedGraph files matching: $BEDGRAPH_PATTERN"
file_count=0

for bedgraph in $BEDGRAPH_PATTERN; do
    if [[ -f "$bedgraph" ]]; then
        convert_to_bigwig "$bedgraph"
        ((file_count++))
    fi
done

# Also convert log2FC files if present
for bedgraph in ATAC_log2FC*.bedGraph; do
    if [[ -f "$bedgraph" ]]; then
        echo ""
        echo "Converting log2FC file..."
        
        # Sort and fix coordinates
        sort -k1,1 -k2,2n "$bedgraph" > "${bedgraph}.sorted"
        
        awk -v OFS='\t' '
            NR==FNR {
                chrom_size[$1]=$2
                next
            }
            {
                if ($1 in chrom_size) {
                    if ($3 > chrom_size[$1]) {
                        $3 = chrom_size[$1]
                    }
                    if ($2 < $3) {
                        print
                    }
                }
            }
        ' "$CHROM_SIZES" "${bedgraph}.sorted" > "${bedgraph}.fixed"
        
        # Convert to bigWig
        basename=$(basename "$bedgraph" .bedGraph)
        $BEDGRAPHTOBIGWIG "${bedgraph}.fixed" "$CHROM_SIZES" "${BIGWIG_DIR}/${basename}.bw"
        
        rm "${bedgraph}.sorted" "${bedgraph}.fixed"
        echo "  ✓ Created: ${BIGWIG_DIR}/${basename}.bw"
        ((file_count++))
    fi
done

# Summary
echo ""
echo "============================================"
echo "Conversion completed: $(date)"
echo "Files converted: $file_count"
echo ""
echo "Output files in: $BIGWIG_DIR/"
ls -lh $BIGWIG_DIR/*.bw 2>/dev/null || echo "No bigWig files found"
echo ""
echo "Next step: Run 03_computematrices.sh"
