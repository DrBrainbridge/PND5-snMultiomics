---
title: "PND5 Control/DES Integrated RNA/ATAC"
output: html_notebook
editor_options: 
  chunk_output_type: inline
---

```{r load library}
library(Seurat)
library(SeuratDisk)
library(Signac, lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
library(irlba, lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(MAST)
library(celldex, lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
library(SingleR)
library(scran)
library(scMerge)
library(tidyverse)
library(limma)
library(topGO)
library(EnsDb.Mmusculus.v79)
library(BSgenome)
library("BSgenome.Mmusculus.UCSC.mm10", lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
library(EnhancedVolcano)
library(liana)
library(readxl)
library(VennDiagram)


library('MEDIPS', lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
library(rtracklayer)
library(hictoolsr)
library(dbscan)
```

Set up Seurats:
```{r Set up Seurats}
## Input files and create Seurat----
## load in RNAseq data----
PND5C.rna <- Read10X_h5("filtered_feature_bc_matrix.h5")
PND5HD.rna <- Read10X_h5("filtered_feature_bc_matrix.h5")
## load in ATAC fragments----
PND5C.ATAC <- "atac_fragments.tsv.gz"
PND5HD.ATAC <- "atac_fragments.tsv.gz"

##put control together----
options(Seurat.object.assay.version = "v5")

metadata <- read.csv(
  file = "per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1)

### create Seurat object for RNA adata----
PND5C <- CreateSeuratObject(
  counts = PND5C.rna$'Gene Expression',
  assay = "RNA",
  meta.data = metadata
)

###get gene chromatin annotations for mm10----
#get annotations from package
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

#create ATAC assay and add to Seurat object
PND5C[["ATAC"]] <- CreateChromatinAssay(
  counts = PND5C.rna$Peaks,
  sep = c(":", "-"),
  fragments = PND5C.ATAC,
  annotation = annotation
)

#Check the Seurat
PND5C
#124018 genes in 20000 cells

### Peak Calling----
DefaultAssay(PND5C) <- "ATAC"

#Call peaks using MACS2
peaksC <- CallPeaks(PND5C)

#remove peaks on nonstandard chromosomes 
peaksC <- keepStandardChromosomes(peaksC, pruning.mode = "coarse")

#remove peaks in genomic blacklist regions
peaksC <- subsetByOverlaps(x = peaksC, ranges = blacklist_mm10, invert = TRUE)

#quantify counts in each peak
macs2_counts_C <- FeatureMatrix(
  fragments = Fragments(PND5C), 
  features = peaksC,
  cells = colnames(PND5C)
)

#create a new assay using the MACS2 peak set & add to Seurat object
PND5C[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_C,
  fragments = PND5C.ATAC,
  annotation = annotation
)

### Quality Control----
DefaultAssay(PND5C) <- "ATAC"
PND5C <- NucleosomeSignal(PND5C)
PND5C <- TSSEnrichment(PND5C)
#expect enrichment of peaks at TSS

PND5C[['pct_reads_in_peaks']] <- PND5C$atac_peak_region_fragments / PND5C$atac_fragments * 100 
PND5C$pct_reads_in_peaks[is.na(PND5C$pct_reads_in_peaks)] <- 0
#usually want at least 15% fragments in

PND5C[['blacklist_fraction']] <- FractionCountsInRegion(
  object = PND5C,
  assay = 'peaks',
  regions = blacklist_mm10)
PND5C[['blacklistratio']] <- PND5C$blacklist_fraction / PND5C$atac_peak_region_fragments * 100
PND5C$blacklistratio[is.na(PND5C$blacklistratio)] <- 0
#blacklist regions defined by ENCODE with high signal in NGS
#cells with high ratio should be thrown out - technical artifacts

DefaultAssay(PND5C) <- "RNA"

PND5C[['percent.mt']] <- PercentageFeatureSet(PND5C, pattern = "^mt-")
#should be none because nucleus sequencing

PND5C[['percent.rp']] <- PercentageFeatureSet(PND5C, pattern = "^Rp")

VlnPlot(
  object = PND5C,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

###save before filter ----
saveRDS(PND5C, file = "PND5C_beforeQCfilter.Rds", compress = TRUE)

##put HighDES together----
rm(list=ls())

metadata <- read.csv(
  file = "per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1)

## create Seurat object for RNA adata----
PND5HD <- CreateSeuratObject(
  counts = PND5HD.rna$'Gene Expression',
  assay = "RNA",
  meta.data = metadata
)

## get gene chromatin annotations for mm10----

#get annotations from package
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

#create ATAC assay and add to Seurat object
PND5HD[["ATAC"]] <- CreateChromatinAssay(
  counts = PND5HD.rna$Peaks,
  sep = c(":", "-"),
  fragments = PND5HD.ATAC,
  annotation = annotation
)

# Peak Calling----
DefaultAssay(PND5HD) <- "ATAC"

#Call peaks using MACS2
peaksHD <- CallPeaks(PND5HD)

#remove peaks on nonstandard chromosomes 
peaksHD <- keepStandardChromosomes(peaksHD, pruning.mode = "coarse")

#remove peaks in genomic blacklist regions
peaksC <- subsetByOverlaps(x = peaksHD, ranges = blacklist_mm10, invert = TRUE)

#quantify counts in each peak
macs2_counts_HD <- FeatureMatrix(
  fragments = Fragments(PND5HD), 
  features = peaksHD,
  cells = colnames(PND5HD)
)

#create a new assay using the MACS2 peak set & add to Seurat object
PND5HD[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_HD,
  fragments = PND5HD.ATAC,
  annotation = annotation
)

# Quality Control----
DefaultAssay(PND5HD) <- "ATAC"
PND5HD <- NucleosomeSignal(PND5HD)
PND5HD <- TSSEnrichment(PND5HD)
#expect enrichment of peaks at TSS

PND5HD[['pct_reads_in_peaks']] <- PND5HD$atac_peak_region_fragments / PND5HD$atac_fragments * 100 
PND5HD$pct_reads_in_peaks[is.na(PND5HD$pct_reads_in_peaks)] <- 0
#usually want at least 15% fragments in

PND5HD[['blacklist_fraction']] <- FractionCountsInRegion(
  object = PND5HD,
  assay = 'peaks',
  regions = blacklist_mm10)
PND5HD[['blacklistratio']] <- PND5HD$blacklist_fraction / PND5HD$atac_peak_region_fragments * 100
PND5HD$blacklistratio[is.na(PND5HD$blacklistratio)] <- 0
#blacklist regions defined by ENCODE with high signal in NGS
#cells with high ratio should be thrown out - technical artifacts

DefaultAssay(PND5HD) <- "RNA"

PND5HD[['percent.mt']] <- PercentageFeatureSet(PND5HD, pattern = "^mt-")
#should be none because nucleus sequencing

PND5HD[['percent.rp']] <- PercentageFeatureSet(PND5HD, pattern = "^Rp")

VlnPlot(
  object = PND5HD,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

###save before filter ----
saveRDS(PND5HD, file = "PND5HD_beforeQCfilter.Rds", compress = TRUE)

# filter out low quality cells----
PND5C <- readRDS("PND5C_beforeQCfilter.Rds")
PND5HD <- readRDS("PND5HD_beforeQCfilter.Rds")

PND5C <- subset(
  x = PND5C,
  subset = 
    nCount_RNA < 10000 &
    nCount_RNA > 500 &
    nCount_ATAC < 10000 &
    nCount_ATAC > 500 &
    TSS.enrichment < 20 &
    TSS.enrichment > 1 &  
    nucleosome_signal < 0.6 &
    nucleosome_signal > 0.25 &
    percent.mt < 10 &
    percent.rp < 15 &
    pct_reads_in_peaks > 15
)
PND5C

PND5HD <- subset(
  x = PND5HD,
  subset = 
    nCount_RNA < 10000 &
    nCount_RNA > 500 &
    nCount_ATAC < 10000 &
    nCount_ATAC > 500 &
    TSS.enrichment < 20 &
    TSS.enrichment > 1 &  
    nucleosome_signal < 0.6 &
    nucleosome_signal > 0.25 &
    percent.mt < 10 &
    percent.rp < 15 &
    pct_reads_in_peaks > 15
)
PND5HD

## Save Seurat of filtered non processed----
saveRDS(PND5C, file = "PND5C_minanalysis.Rds", compress = TRUE)
saveRDS(PND5HD, file = "PND5HD_minanalysis.Rds", compress = TRUE)
```