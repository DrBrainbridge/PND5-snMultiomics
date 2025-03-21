---
title: "PND5 Control/DES Integrated RNA/ATAC"
output: html_notebook
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

```{r Initial Analysis Pipeline}
PND5C <- readRDS(file = "PND5C_minanalysis.Rds")
PND5HD <- readRDS(file = "PND5HD_minanalysis.Rds")

# Gene Expression data processing----
DefaultAssay(PND5C) <- "RNA"
DefaultAssay(PND5HD) <- "RNA"

##SC Transform (v2)----
options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb

PND5C <- SCTransform(PND5C, vst.flavor = "v2", verbose = FALSE, vars.to.regress = "percent.rp")
PND5C <- RunPCA(PND5C)
PND5HD <- SCTransform(PND5HD, vst.flavor = "v2", verbose = FALSE, vars.to.regress = "percent.rp")
PND5HD <- RunPCA(PND5HD)

##make plots to visualize dimensionality----
VizDimLoadings(PND5C, dims = 1:2, reduction = 'pca')
DimHeatmap(PND5C, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(PND5C, ndims = 50)

VizDimLoadings(PND5HD, dims = 1:2, reduction = 'pca')
DimHeatmap(PND5HD, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(PND5HD, ndims = 50)

##run UMAP based on RNA----
PND5C <- RunUMAP(PND5C, dims = 1:30)
PND5C <- FindNeighbors(PND5C, dims = 1:30)
PND5C <- FindClusters(PND5C, resolution = 0.4)
PND5HD <- RunUMAP(PND5HD, dims = 1:30)
PND5HD <- FindNeighbors(PND5HD, dims = 1:30)
PND5HD <- FindClusters(PND5HD, resolution = 0.4)

###print UMAP with clusters----
DimPlot(PND5C, label = TRUE, repel=TRUE)

DimPlot(PND5HD, label = TRUE, repel=TRUE)

```
