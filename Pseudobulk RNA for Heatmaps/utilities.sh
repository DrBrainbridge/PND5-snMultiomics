#!/bin/bash
# utilities.sh
# Helper functions and utilities for ATAC-seq analysis pipeline

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

# Check if required software is installed
check_software() {
    local missing=0
    
    echo "Checking required software..."
    
    # Check R
    if command -v R &> /dev/null; then
        echo "✓ R found: $(R --version | head -1)"
    else
        echo "✗ R not found"
        ((missing++))
    fi
    
    # Check deepTools
    if command -v computeMatrix &> /dev/null; then
        echo "✓ deepTools found: $(computeMatrix --version)"
    else
        echo "✗ deepTools not found"
        ((missing++))
    fi
    
    if command -v plotHeatmap &> /dev/null; then
        echo "✓ plotHeatmap found"
    else
        echo "✗ plotHeatmap not found"
        ((missing++))
    fi
    
    # Check UCSC tools
    if command -v bedGraphToBigWig &> /dev/null; then
        echo "✓ bedGraphToBigWig found"
    else
        echo "✗ bedGraphToBigWig not found"
        ((missing++))
    fi
    
    if [[ $missing -gt 0 ]]; then
        echo ""
        echo "ERROR: $missing required tool(s) missing"
        echo "Please install missing software before proceeding"
        return 1
    else
        echo ""
        echo "✓ All required software found"
        return 0
    fi
}

# Validate BED file format
validate_bed() {
    local bed_file=$1
    
    echo "Validating BED file: $bed_file"
    
    if [[ ! -f "$bed_file" ]]; then
        echo "ERROR: File not found"
        return 1
    fi
    
    # Check number of columns
    local ncols=$(head -1 "$bed_file" | awk '{print NF}')
    if [[ $ncols -lt 3 ]]; then
        echo "ERROR: BED file must have at least 3 columns"
        return 1
    fi
    echo "✓ Columns: $ncols"
    
    # Check for valid chromosomes
    local chr_count=$(cut -f1 "$bed_file" | grep -c "^chr")
    local total_lines=$(wc -l < "$bed_file")
    echo "✓ Lines with 'chr' prefix: $chr_count / $total_lines"
    
    # Check strand column if present
    if [[ $ncols -ge 6 ]]; then
        local strand_check=$(cut -f6 "$bed_file" | grep -c "^[+-]")
        echo "✓ Valid strand annotations: $strand_check / $total_lines"
    fi
    
    echo "✓ BED file appears valid"
    return 0
}

# Validate bigWig file
validate_bigwig() {
    local bw_file=$1
    
    echo "Validating bigWig: $bw_file"
    
    if [[ ! -f "$bw_file" ]]; then
        echo "ERROR: File not found"
        return 1
    fi
    
    if command -v bigWigInfo &> /dev/null; then
        bigWigInfo "$bw_file"
        echo "✓ bigWig file is valid"
    else
        # Basic check - bigWig files start with specific magic number
        local magic=$(od -An -t x1 -N 4 "$bw_file" | tr -d ' ')
        if [[ "$magic" == "26fc8f88" ]]; then
            echo "✓ File has valid bigWig signature"
        else
            echo "ERROR: File does not appear to be a valid bigWig"
            return 1
        fi
    fi
    
    return 0
}

# ============================================================================
# FILE CONVERSION FUNCTIONS
# ============================================================================

# Extract TSS from GTF file
gtf_to_tss_bed() {
    local gtf_file=$1
    local output_bed=$2
    
    echo "Converting GTF to TSS BED: $gtf_file -> $output_bed"
    
    awk -F'\t' '
    $3=="gene" {
        split($9, a, ";")
        for(i in a) {
            if(a[i] ~ /gene_name/) {
                split(a[i], b, "\"")
                gene=b[2]
            }
        }
        if($7=="+") {
            print $1"\t"$4-1"\t"$4"\t"gene"\t0\t"$7
        } else {
            print $1"\t"$5-1"\t"$5"\t"gene"\t0\t"$7
        }
    }' "$gtf_file" > "$output_bed"
    
    echo "✓ Created $output_bed with $(wc -l < "$output_bed") genes"
}

# Convert between chromosome naming conventions
convert_chr_names() {
    local input_file=$1
    local output_file=$2
    local mode=$3  # "add" or "remove"
    
    if [[ "$mode" == "add" ]]; then
        echo "Adding 'chr' prefix..."
        sed 's/^\([0-9XYM]\)/chr\1/' "$input_file" > "$output_file"
    elif [[ "$mode" == "remove" ]]; then
        echo "Removing 'chr' prefix..."
        sed 's/^chr//' "$input_file" > "$output_file"
    else
        echo "ERROR: mode must be 'add' or 'remove'"
        return 1
    fi
    
    echo "✓ Converted: $output_file"
}

# ============================================================================
# DATA INSPECTION FUNCTIONS
# ============================================================================

# Get summary statistics from bigWig
bigwig_stats() {
    local bw_file=$1
    
    if ! command -v bigWigInfo &> /dev/null; then
        echo "ERROR: bigWigInfo not found"
        return 1
    fi
    
    echo "BigWig Statistics for: $bw_file"
    echo "----------------------------------------"
    bigWigInfo "$bw_file"
}

# Compare gene lists
compare_gene_lists() {
    local list1=$1
    local list2=$2
    
    echo "Comparing gene lists..."
    echo "List 1: $list1"
    echo "List 2: $list2"
    echo ""
    
    # Extract gene names
    tail -n +2 "$list1" | cut -d',' -f1 | sort > /tmp/genes1.txt
    tail -n +2 "$list2" | cut -d',' -f1 | sort > /tmp/genes2.txt
    
    local count1=$(wc -l < /tmp/genes1.txt)
    local count2=$(wc -l < /tmp/genes2.txt)
    local common=$(comm -12 /tmp/genes1.txt /tmp/genes2.txt | wc -l)
    local unique1=$(comm -23 /tmp/genes1.txt /tmp/genes2.txt | wc -l)
    local unique2=$(comm -13 /tmp/genes1.txt /tmp/genes2.txt | wc -l)
    
    echo "List 1: $count1 genes"
    echo "List 2: $count2 genes"
    echo "Common: $common genes"
    echo "Unique to list 1: $unique1 genes"
    echo "Unique to list 2: $unique2 genes"
    
    # Clean up
    rm /tmp/genes1.txt /tmp/genes2.txt
}

# Extract genes in specific region
extract_genes_in_region() {
    local bed_file=$1
    local chr=$2
    local start=$3
    local end=$4
    
    echo "Extracting genes from $bed_file in region ${chr}:${start}-${end}"
    
    awk -v chr="$chr" -v start="$start" -v end="$end" \
        '$1==chr && $2>=start && $3<=end {print $4}' "$bed_file"
}

# ============================================================================
# QUALITY CONTROL FUNCTIONS
# ============================================================================

# Check matrix file
check_matrix() {
    local matrix_file=$1
    
    echo "Checking matrix: $matrix_file"
    
    if [[ ! -f "$matrix_file" ]]; then
        echo "ERROR: File not found"
        return 1
    fi
    
    # Check if it's gzipped
    if file "$matrix_file" | grep -q "gzip"; then
        local lines=$(zcat "$matrix_file" | wc -l)
        echo "✓ Compressed matrix with $lines lines"
        
        # Show first few lines
        echo ""
        echo "First 5 lines:"
        zcat "$matrix_file" | head -5
    else
        local lines=$(wc -l < "$matrix_file")
        echo "✓ Uncompressed matrix with $lines lines"
        
        # Show first few lines
        echo ""
        echo "First 5 lines:"
        head -5 "$matrix_file"
    fi
}

# Generate QC report
generate_qc_report() {
    local output_file="QC_report.txt"
    
    echo "Generating QC report..."
    
    cat > "$output_file" << EOF
ATAC-seq Analysis Pipeline - Quality Control Report
Generated: $(date)
========================================================

INPUT FILES:
-----------
EOF
    
    # Check for bigWig files
    if ls bigwig/*.bw 1> /dev/null 2>&1; then
        echo "" >> "$output_file"
        echo "BigWig files:" >> "$output_file"
        for bw in bigwig/*.bw; do
            size=$(du -h "$bw" | cut -f1)
            echo "  - $(basename "$bw"): $size" >> "$output_file"
        done
    fi
    
    # Check for matrices
    if ls matrices/*.gz 1> /dev/null 2>&1; then
        echo "" >> "$output_file"
        echo "Matrix files:" >> "$output_file"
        for mat in matrices/*.gz; do
            size=$(du -h "$mat" | cut -f1)
            lines=$(zcat "$mat" | wc -l)
            echo "  - $(basename "$mat"): $size ($lines lines)" >> "$output_file"
        done
    fi
    
    # Check for plots
    if ls heatmaps/*.pdf 1> /dev/null 2>&1; then
        echo "" >> "$output_file"
        echo "Heatmap files:" >> "$output_file"
        for plot in heatmaps/*.pdf; do
            size=$(du -h "$plot" | cut -f1)
            echo "  - $(basename "$plot"): $size" >> "$output_file"
        done
    fi
    
    echo "" >> "$output_file"
    echo "========================================================" >> "$output_file"
    echo "QC report complete: $output_file" >> "$output_file"
    
    cat "$output_file"
    echo ""
    echo "✓ QC report saved to: $output_file"
}

# ============================================================================
# CLEANUP FUNCTIONS
# ============================================================================

# Remove intermediate files
cleanup_intermediates() {
    echo "Cleaning up intermediate files..."
    
    # Remove temporary bedGraph files
    if ls *.bedGraph 1> /dev/null 2>&1; then
        echo "Removing bedGraph files..."
        rm -v *.bedGraph
    fi
    
    # Remove sorted/fixed temporary files
    if ls *.sorted 1> /dev/null 2>&1; then
        echo "Removing sorted files..."
        rm -v *.sorted
    fi
    
    if ls *.fixed 1> /dev/null 2>&1; then
        echo "Removing fixed files..."
        rm -v *.fixed
    fi
    
    echo "✓ Cleanup complete"
}

# Archive results
archive_results() {
    local archive_name="atac_analysis_$(date +%Y%m%d_%H%M%S).tar.gz"
    
    echo "Creating archive: $archive_name"
    
    tar -czf "$archive_name" \
        bigwig/ \
        matrices/ \
        heatmaps/ \
        profiles/ \
        reference/ \
        *.txt \
        *.md \
        2>/dev/null
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Archive created: $archive_name ($(du -h "$archive_name" | cut -f1))"
    else
        echo "ERROR: Archive creation failed"
        return 1
    fi
}

# ============================================================================
# USAGE INFORMATION
# ============================================================================

show_usage() {
    cat << EOF
ATAC-seq Pipeline Utilities
============================

Available functions:

Validation:
  check_software              - Check if required tools are installed
  validate_bed FILE           - Validate BED file format
  validate_bigwig FILE        - Check bigWig file validity

Conversion:
  gtf_to_tss_bed IN OUT       - Convert GTF to TSS BED
  convert_chr_names IN OUT MODE - Add/remove 'chr' prefix (MODE: add|remove)

Inspection:
  bigwig_stats FILE           - Get bigWig statistics
  compare_gene_lists F1 F2    - Compare two gene lists
  extract_genes_in_region BED CHR START END - Extract genes in region
  check_matrix FILE           - Inspect matrix file

Quality Control:
  generate_qc_report          - Generate QC summary
  cleanup_intermediates       - Remove temporary files
  archive_results             - Create results archive

Usage:
  source utilities.sh
  check_software
  validate_bed reference/genes.bed
  generate_qc_report

EOF
}

# ============================================================================
# MAIN
# ============================================================================

# If script is run directly, show usage
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    show_usage
fi
