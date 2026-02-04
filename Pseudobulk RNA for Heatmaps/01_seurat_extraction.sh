#!/bin/bash
# 01_seurat_extraction.sh
# Extract ATAC-seq accessibility data from Seurat object and create bedGraph files

set -euo pipefail

echo "Starting Seurat data extraction"
echo "Start time: $(date)"

# ============================================================================
# CONFIGURATION - Edit these variables for your analysis
# ============================================================================

# Path to your Seurat object (RDS file)
SEURAT_OBJECT="your_seurat_object.rds"

# Output prefix for bedGraph files
OUTPUT_PREFIX="ATAC"

# Treatment/condition column name in Seurat metadata
METADATA_COLUMN="Treatment"

# Treatment groups to extract
CONTROL_GROUP="Control"
TREATMENT_GROUP="Treatment"

# Minimum accessibility threshold (removes low signal peaks)
MIN_ACCESSIBILITY=0.01

# ============================================================================
# R SCRIPT - Extract and process ATAC-seq data
# ============================================================================

cat > extract_seurat_data.R << 'R_SCRIPT_END'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(data.table)
})

# Load command line arguments or use defaults
args <- commandArgs(trailingOnly = FALSE)
seurat_file <- Sys.getenv("SEURAT_OBJECT", "your_seurat_object.rds")
output_prefix <- Sys.getenv("OUTPUT_PREFIX", "ATAC")
metadata_col <- Sys.getenv("METADATA_COLUMN", "Treatment")
control_group <- Sys.getenv("CONTROL_GROUP", "Control")
treatment_group <- Sys.getenv("TREATMENT_GROUP", "Treatment")
min_accessibility <- as.numeric(Sys.getenv("MIN_ACCESSIBILITY", "0.01"))

cat("Configuration:\n")
cat("  Seurat object:", seurat_file, "\n")
cat("  Metadata column:", metadata_col, "\n")
cat("  Control group:", control_group, "\n")
cat("  Treatment group:", treatment_group, "\n\n")

# Load Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(seurat_file)

# Set ATAC assay as default
DefaultAssay(seurat_obj) <- ifelse("peaks" %in% names(seurat_obj@assays), "peaks", "ATAC")

cat("Seurat object info:\n")
cat("  Default assay:", DefaultAssay(seurat_obj), "\n")
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Features:", nrow(seurat_obj), "\n\n")

# Get normalized data - handle both old and new Seurat versions
current_assay <- seurat_obj[[DefaultAssay(seurat_obj)]]

if (inherits(current_assay, "ChromatinAssay") && !("layers" %in% slotNames(current_assay))) {
  # Older Seurat/Signac version
  cat("Using legacy Seurat version data access\n")
  if ("data" %in% slotNames(current_assay) && length(slot(current_assay, "data")) > 0) {
    normalized_data <- slot(current_assay, "data")
  } else {
    cat("Applying TF-IDF normalization\n")
    seurat_obj <- RunTFIDF(seurat_obj)
    normalized_data <- GetAssayData(seurat_obj, slot = "data")
  }
} else {
  # Newer Seurat version with layers
  if ("layers" %in% slotNames(current_assay) && "data" %in% names(current_assay@layers)) {
    normalized_data <- LayerData(seurat_obj, layer = "data")
  } else {
    seurat_obj <- RunTFIDF(seurat_obj)
    normalized_data <- GetAssayData(seurat_obj, slot = "data")
  }
}

cat("Normalized data dimensions:", dim(normalized_data), "\n\n")

# Extract and parse peak coordinates
peak_names <- rownames(normalized_data)
peak_coords <- data.table(peak_name = peak_names)

# Parse coordinates (supports chr:start-end or chr_start_end formats)
if (any(grepl(":", peak_names)) && any(grepl("-", peak_names))) {
  # Format: chr1:1000-2000
  peak_coords[, c("chr", "coords") := tstrsplit(peak_name, ":", fixed=TRUE)]
  peak_coords[, c("start", "end") := tstrsplit(coords, "-", fixed=TRUE)]
  peak_coords$coords <- NULL
} else if (any(grepl("_", peak_names))) {
  # Format: chr1_1000_2000
  parts <- strsplit(peak_names, "_")
  peak_coords$chr <- sapply(parts, `[`, 1)
  peak_coords$start <- sapply(parts, `[`, 2)
  peak_coords$end <- sapply(parts, `[`, 3)
} else {
  stop("Could not parse peak coordinate format")
}

peak_coords$start <- as.numeric(peak_coords$start)
peak_coords$end <- as.numeric(peak_coords$end)
peak_coords <- peak_coords[!is.na(start) & !is.na(end)]

cat("Valid peaks:", nrow(peak_coords), "\n\n")

# Check for treatment groups
if (metadata_col %in% colnames(seurat_obj@meta.data)) {
  treatment_groups <- unique(seurat_obj@meta.data[[metadata_col]])
  cat("Treatment groups found:", paste(treatment_groups, collapse=", "), "\n")
} else {
  cat("Warning: Metadata column", metadata_col, "not found\n")
  cat("Available columns:", paste(colnames(seurat_obj@meta.data), collapse=", "), "\n")
}

# Function to create bedGraph files
create_bedgraph <- function(treatment_group, output_prefix) {
  cat("\nProcessing group:", treatment_group, "\n")
  
  # Filter cells by treatment
  if (metadata_col %in% colnames(seurat_obj@meta.data)) {
    selected_cells <- which(seurat_obj@meta.data[[metadata_col]] == treatment_group)
    if (length(selected_cells) == 0) {
      cat("Warning: No cells found for", treatment_group, "\n")
      return(NULL)
    }
    group_data <- normalized_data[peak_coords$peak_name, selected_cells, drop=FALSE]
    cat("  Selected cells:", length(selected_cells), "\n")
  } else {
    group_data <- normalized_data[peak_coords$peak_name, , drop=FALSE]
  }
  
  # Calculate mean accessibility per peak
  mean_accessibility <- Matrix::rowMeans(group_data)
  
  # Create bedGraph
  bedgraph_df <- data.table(
    chr = peak_coords$chr,
    start = peak_coords$start,
    end = peak_coords$end,
    score = mean_accessibility
  )
  
  # Filter and sort
  bedgraph_df <- bedgraph_df[is.finite(score) & !is.na(chr) & score > min_accessibility]
  bedgraph_df <- bedgraph_df[order(chr, start)]
  
  # Keep only standard chromosomes (adjust for your organism)
  standard_chrs <- paste0("chr", c(1:19, "X", "Y", "M"))  # Mouse
  bedgraph_df <- bedgraph_df[chr %in% standard_chrs]
  
  cat("  Final peaks:", nrow(bedgraph_df), "\n")
  cat("  Score range:", round(min(bedgraph_df$score), 4), "to", 
      round(max(bedgraph_df$score), 4), "\n")
  
  # Write bedGraph
  output_file <- paste0(output_prefix, "_", treatment_group, "_normalized.bedGraph")
  track_line <- paste0('track type=bedGraph name="', treatment_group, 
                       '" description="Normalized ATAC-seq"')
  cat(track_line, "\n", file = output_file)
  fwrite(bedgraph_df, output_file, sep = "\t", append = TRUE, col.names = FALSE)
  
  cat("  Written:", output_file, "\n")
  return(output_file)
}

# Generate bedGraph files
if (metadata_col %in% colnames(seurat_obj@meta.data)) {
  control_file <- create_bedgraph(control_group, output_prefix)
  treatment_file <- create_bedgraph(treatment_group, output_prefix)
  
  # Create log2FC file
  if (!is.null(control_file) && !is.null(treatment_file)) {
    cat("\nCreating log2FC file...\n")
    control_data <- Matrix::rowMeans(
      normalized_data[peak_coords$peak_name, 
                     seurat_obj@meta.data[[metadata_col]] == control_group, 
                     drop=FALSE])
    treatment_data <- Matrix::rowMeans(
      normalized_data[peak_coords$peak_name, 
                     seurat_obj@meta.data[[metadata_col]] == treatment_group, 
                     drop=FALSE])
    
    log2fc_data <- data.table(
      chr = peak_coords$chr,
      start = peak_coords$start,
      end = peak_coords$end,
      log2fc = log2((treatment_data + 0.01) / (control_data + 0.01))
    )
    
    log2fc_data <- log2fc_data[is.finite(log2fc)]
    output_log2fc <- paste0(output_prefix, "_log2FC_", treatment_group, "_vs_", control_group, ".bedGraph")
    fwrite(log2fc_data, output_log2fc, sep = "\t", col.names = FALSE)
    cat("  Written:", output_log2fc, "\n")
  }
}

cat("\nSeurat extraction completed successfully\n")
R_SCRIPT_END

# Export configuration variables for R script
export SEURAT_OBJECT
export OUTPUT_PREFIX
export METADATA_COLUMN
export CONTROL_GROUP
export TREATMENT_GROUP
export MIN_ACCESSIBILITY

# Execute R script
Rscript extract_seurat_data.R

# Validate outputs
echo ""
echo "Validating output files..."
for file in ${OUTPUT_PREFIX}_*.bedGraph; do
    if [[ -f "$file" ]]; then
        echo "✓ Created $file ($(wc -l < "$file") lines)"
    fi
done

echo ""
echo "Extraction completed: $(date)"
echo "Next step: Run 02_bedgraph_to_bigwig.sh"
