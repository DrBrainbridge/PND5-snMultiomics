---
title: "PND5 Control/DES Integrated RNA/ATAC"
output: html_notebook
---

This is an example workflow going from single cell RNA and ATAC files through
complete analysis. Paths are simplified as to not reflect server structure. All
analysis presented here was run in R 4.3.1 on the NIH server. Specific feature
and other plots are not included, but can be inferred from examples left.

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
library('scplotter')
library(TFBSTools)
  #devtools::install_github("da-bar/JASPAR2022")
library(JASPAR2022)
  #BiocManager::install("motifmatchr")
library(motifmatchr)
library(ggseqlogo)
  #BiocManager::install("GreenleafLab/chromVAR", lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
library(chromVAR, lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')

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

That may be too few DES cells, but we will proceed in doing dimensionality reduction.
Don't need to do this, but wanted to see if dimensions behaved similarly between
control and DES.

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
Everything looks normal. Any strangeness in the integrated dataset comes from
biological differences. 

Make integrated Seurat with filtered non processed Seurats
```{r Integrate}
PND5C <- readRDS(file = "PND5C_minanalysis.Rds")
PND5HD <- readRDS(file = "PND5HD_minanalysis.Rds")

#add treatment label
Clabel <- rep("Control", ncol(x = PND5C))
PND5C$Treatment <- Clabel
HDlabel <- rep("HighDES", ncol(x = PND5HD))
PND5HD$Treatment <- HDlabel

head(PND5HD,1)
head(PND5C,1)

##SC Transform (v2)----
options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb

PND5C <- SCTransform(PND5C, vst.flavor = "v2", verbose = FALSE, vars.to.regress = "percent.rp")
PND5C <- RunPCA(PND5C)
PND5HD <- SCTransform(PND5HD, vst.flavor = "v2", verbose = FALSE, vars.to.regress = "percent.rp")
PND5HD <- RunPCA(PND5HD)

#back to integration pipeline
CODESobjlist <- list(Control = PND5C, HighDES = PND5HD)
Cfeatures <- SelectIntegrationFeatures(CODESobjlist, nfeatures = 3000)
CODESobjlist <- PrepSCTIntegration(CODESobjlist, anchor.features = Cfeatures)
CODES.anchors <- FindIntegrationAnchors(CODESobjlist, normalization.method = "SCT",
                                        anchor.features = Cfeatures)
PND5.CODES.sct <- IntegrateData(anchorset = CODES.anchors,
                                 normalization.method = "SCT")

table(PND5.CODES.sct$Treatment)


PND5.CODES.sct <- RunPCA(PND5.CODES.sct, verbose = FALSE)

saveRDS(PND5.CODES.sct, file = "PND5CODES_intbeforeUMAP.Rds", compress = TRUE)

##make plots to visualize dimensionality----
ElbowPlot(PND5.CODES.sct, ndims = 50)

##finish dimensionality reduction & visualize
PND5.CODES.sct <- RunUMAP(PND5.CODES.sct, reduction = "pca", dims = 1:30)

PND5.CODES.sct <- FindNeighbors(PND5.CODES.sct, reduction = "pca", dims = 1:30)

PND5.CODES.sct <- FindClusters(PND5.CODES.sct, resolution = 0.4)

saveRDS(PND5.CODES.sct, file = "PND5CODES_SCT2integrate.Rds", compress = TRUE)

DimPlot(PND5.CODES.sct, label = TRUE, repel=TRUE, group.by = "integrated_snn_res.0.4")
DimPlot(PND5.CODES.sct, label = FALSE, repel=TRUE, group.by = "integrated_snn_res.0.4") & NoAxes() & NoLegend() & ggtitle(NULL)
DimPlot(PND5.CODES.sct, label = TRUE, repel=TRUE, group.by = 'Treatment')
DimPlot(PND5.CODES.sct, label = FALSE, repel=TRUE, group.by = 'Treatment') & NoAxes() & NoLegend() & ggtitle(NULL)

```

Now that integrated Seurat Control/DES is built, do differential analysis. 
Results should match bulk sequencing.

```{r Differential analysis}
PND5.CODES.sct <- readRDS("PND5CODES_SCT2integrate.Rds")

##find markers-----
#by treatment
Idents(PND5.CODES.sct) <- "Treatment"

PND5.CODES.sct <- PrepSCTFindMarkers(PND5.CODES.sct)

PND5.CODES.sct.markers <- FindMarkers(PND5.CODES.sct, assay = "SCT",
                                      ident.1="HighDES", ident.2="Control")

#saved with right order so DES is +
write.csv(PND5.CODES.sct.markers, "CoDESMarkerstreatment.csv", row.names = TRUE)

##volcano plot----
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
PND5.CODES.sct.markers <- read.csv("CoDESMarkerstreatment.csv")
head(PND5.CODES.sct.markers)

PND5.CODES.sct.markers <- PND5.CODES.sct.markers[order(PND5.CODES.sct.markers$avg_log2FC),]
topCon.markers <- head(PND5.CODES.sct.markers, 50)
topDES.markers <- tail(PND5.CODES.sct.markers, 50)
  
EnhancedVolcano(PND5.CODES.sct.markers, lab = PND5.CODES.sct.markers$X,
                x = 'avg_log2FC', y = 'p_val_adj', title = 'Control vs DES',
                pCutoff = 0.01, pointSize = 1.0, labSize = 5.0,
                ylim = c(0, -log10(10e-300)), drawConnectors = FALSE,
                selectLab = c(topCon.markers$X, topDES.markers$X))
#not ideal visualization because differences are so large.

```

Now I will annotate the clusters. 

```{r Annotate Integrated Full Uterus}
PND5.CODES.sct <- readRDS("PND5CODES_SCT2integrate.Rds")
options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb

##annotate----
###find markers per cluster----
####crashes R if using RNA, works if using SCT!-----
#if crashing, reset environment and redo (quit q() and do not save env)
Idents(PND5.CODES.sct) <- 'integrated_snn_res.0.4'
PND5.CODES.sct.markers.all <- FindAllMarkers(PND5.CODES.sct, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.CODES.sct.markers.all, "CoDESMarkers_byCluster.csv", row.names = TRUE)

###print plots of cell markers ----
DefaultAssay(PND5.CODES.sct) <- "SCT"
Idents(PND5.CODES.sct) <- "integrated_snn_res.0.4"

####Feature plots----
#####----PND3 markers----
#PND3 markers from Jia & Zhao
#Vim = Mesenchyme or Stroma (by Spencer paper)
#Hoxa10 = Mesenchyme (uterine)
#Top2a = proliferative
#Epcam = epithelium
#Msln = mesothelium
#Pecam1 = endothelium
#Rgs5 = pericyte
#Lyz2 = myeloid

FeaturePlot(object = PND5.CODES.sct,
            features = c('Hoxa10', 'Vim', 'Top2a', 'Epcam', 'Msln', 'Pecam1', 'Rgs5', 'Lyz2'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)

FeaturePlot(object = PND5.CODES.sct,
            features = c('Hoxa10', 'Vim', 'Top2a', 'Epcam', 'Msln', 'Pecam1', 'Rgs5', 'Lyz2'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE, keep.scale = 'all')

####Violin and Dot plots----
#####----PND3 markers----
VlnPlot(object = PND5.CODES.sct, features = c('Hoxa10', 'Vim', 'Top2a', 'Epcam', 'Msln', 'Pecam1', 'Rgs5', 'Lyz2'))

DotPlot(object = PND5.CODES.sct, features = c('Hoxa10', 'Vim', 'Top2a', 'Epcam', 'Msln', 'Pecam1', 'Rgs5', 'Lyz2')) + coord_flip()

##annotate----
table(PND5.CODES.sct$integrated_snn_res.0.4)
PND5.CODES.sct <- RenameIdents(object = PND5.CODES.sct,
                                  '0' = 'Mesenchyme',#Hoxa10
                                  '1' = 'Epithelium',#Epcam
                                  '2' = 'Mesenchyme',#Hoxa10 #could also be stromal fibroblast Mfap5
                                  '3' = 'Epithelium',#Epcam 
                                  '4' = 'Endothelium',#Pecam1
                                  '5' = 'Stroma',#Dcn & Mfap5
                                  '6' = 'Mesenchyme',#Hoxa10
                                  '7' = 'Pericyte',#Rgs5
                                  '8' = 'Epithelium',#Epcam
                                  '9' = 'Mesothelium',#Msln
                                  '10' = 'Myeloid'#Lyz2
)

PND5.CODES.sct[["Manual_annotation"]] <- Idents(PND5.CODES.sct)
table(PND5.CODES.sct$Manual_annotation)

##visualize annotated clusters----
Idents(PND5.CODES.sct) <- 'Manual_annotation'
DimPlot(PND5.CODES.sct, label = FALSE, repel=TRUE)

DimPlot(PND5.CODES.sct, label = FALSE, repel=TRUE, split.by = 'Treatment')

DimPlot(PND5.CODES.sct, label = FALSE, repel=TRUE)

##save annotated----
saveRDS(PND5.CODES.sct, file = "PND5CoDES_annotated.Rds", compress = TRUE)

##DEX and heatmap by annotated cluster & treatment----
table(PND5.CODES.sct$Manual_annotation)
table(PND5.CODES.sct$Treatment)

PND5.CODES.sct$TissueTreatment <- paste(PND5.CODES.sct$Manual_annotation, PND5.CODES.sct$Treatment)
table(PND5.CODES.sct$TissueTreatment)

options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb

Idents(PND5.CODES.sct) <- 'TissueTreatment'
PND5.CODES.sct <- PrepSCTFindMarkers(PND5.CODES.sct)

PND5.CODES.sct.markers.all <- FindAllMarkers(PND5.CODES.sct, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.CODES.sct.markers.all, "CoDESMarkers_byTissTrx.csv", row.names = TRUE)

###by tiss/trx----
PND5.CODES.sct.markers <- dplyr::filter(PND5.CODES.sct.markers.all, avg_log2FC > 1.5) %>%
  dplyr::filter(,p_val_adj == 0) %>%
  arrange(,avg_log2FC) %>%
  arrange(,cluster) %>%
  group_by(cluster) %>%
  slice_head(n = 20)

heatmap_tissue <- DoHeatmap(PND5.CODES.sct, group.by='TissueTreatment',
                     features = PND5.CODES.sct.markers$gene, size = 3, angle = 45) 

heatmap_tissue  

##pie chart of cell types====
#https://pwwang.github.io/scplotter/index.html
#devtools::install_github("pwwang/scplotter", force = TRUE)

PND5.CODES.sct <- readRDS("PND5CoDES_annotated.Rds")
sspalette <- c('#F3756E', '#C29B2E', '#4FB448', '#24B891', '#22B4E9', '#9489C1', '#D56DAB')
table(PND5.CODES.sct$Manual_annotation)

CellStatPlot(PND5.CODES.sct, plot = 'ring', group_by = 'Treatment', ident = 'Manual_annotation', palcolor = sspalette)

```

Subset epithelium, mesenchyme, stroma to look at epithelial-mesenchymal transition:
```{r Subset EpiMeStro}
PND5.CODES.sct <- readRDS("PND5CoDES_annotated.Rds")

##reset identities to seurat clusters----
DefaultAssay(PND5.CODES.sct) <- 'SCT'
Idents(PND5.CODES.sct) <- ('Manual_annotation')
head(Idents(PND5.CODES.sct),1)

options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb

##subset clusters of interest----
PND5.CoDES.EpiMeStro<- subset(x = PND5.CODES.sct, idents = c('Epithelium', 'Mesenchyme', 'Stroma'))
table(PND5.CoDES.EpiMeStro$Treatment)
table(PND5.CoDES.EpiMeStro$Manual_annotation)

## standard visualization and clustering----
DefaultAssay(PND5.CoDES.EpiMeStro) <- "RNA"
PND5.CoDES.EpiMeStro <- SCTransform(PND5.CoDES.EpiMeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 
  
DepthCor(PND5.CoDES.EpiMeStro, reduction = 'pca', n = 50)
#looks so much better than before - not much correllation at all 
ElbowPlot(PND5.CoDES.EpiMeStro, ndims = 50)

PND5.CoDES.EpiMeStro <- FindNeighbors(PND5.CoDES.EpiMeStro, dims = c(1:30)) %>% 
  FindClusters(resolution = 0.8) %>%
  RunUMAP(dims = 1:30)
  #original was run on res = 0.4
  #move forward with rse = 0.6 because that makes bridge its own cluster
  #trying 0.8 to get multiple clusters in bridge - works well
  
DimPlot(PND5.CoDES.EpiMeStro, label = FALSE, repel=TRUE)

DimPlot(PND5.CoDES.EpiMeStro, label = FALSE, repel=TRUE, group.by="Treatment")

PND5.CoDES.EpiMeStro[["EpiMeStro_subcluster"]] <- Idents(PND5.CoDES.EpiMeStro)

##check if any cluster has low QC metrics----
Idents(PND5.CoDES.EpiMeStro) <- 'EpiMeStro_subcluster'
VlnPlot(
  object = PND5.CoDES.EpiMeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

saveRDS(PND5.CoDES.EpiMeStro, file = "PND5CoDESEpiMeStro.Rds", compress = TRUE)
```

Ok so that's interesting. There are now two separate connected populations likely going from mesenchyme to epithelium. I will manually annotate and explore:
```{r Annotate EpiMeStro}
PND5.CoDES.EpiMeStro <- readRDS("PND5CoDESEpiMeStro.Rds")

options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb

##annotate----
###find markers per cluster----
####crashes R if using RNA, works if using SCT!-----
#if crashing, reset environment and redo (quit q() and do not save env)
Idents(PND5.CoDES.EpiMeStro) <- 'EpiMeStro_subcluster'
PND5.CoDES.EpiMeStro <- PrepSCTFindMarkers(PND5.CoDES.EpiMeStro)

PND5.CoDES.EpiMeStro.markers.all <- FindAllMarkers(PND5.CoDES.EpiMeStro, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.CoDES.EpiMeStro.markers.all, "CoDES_EpiMeStro_Markers_byCluster_res.csv", row.names = TRUE)

###print plots of cell markers ----
DefaultAssay(PND5.CoDES.EpiMeStro) <- "SCT"
Idents(PND5.CoDES.EpiMeStro) <- 'EpiMeStro_subcluster'

######----Epithelial----
FeaturePlot(object = PND5.CoDES.EpiMeStro,
            features = c('Epcam', 'Cdh1', 'Krt18'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
VlnPlot(object = PND5.CoDES.EpiMeStro, features = c('Epcam', 'Cdh1', 'Krt18'))
DotPlot(object = PND5.CoDES.EpiMeStro, features = c('Epcam', 'Cdh1', 'Krt18')) + coord_flip()

#####----Stroma----
FeaturePlot(object = PND5.CoDES.EpiMeStro,
            features = c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
VlnPlot(object = PND5.CoDES.EpiMeStro, features = c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6'))
dev.off()
DotPlot(object = PND5.CoDES.EpiMeStro, features = c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6')) + coord_flip()

#####----PND3 markers----
#PND3 markers from Jia & Zhao
#Vim = Mesenchyme or Stroma (by Spencer paper)
#Hoxa10 = Mesenchyme
#Cnn1 = smooth muscle
#Tcf21 = smooth muscle
#Pkib = mesenchyme subpopulation
#Pcdh10 = mesenchyme subpopulation

FeaturePlot(object = PND5.CoDES.EpiMeStro,
            features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE, keep.scale = 'all')
DotPlot(object = PND5.CoDES.EpiMeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Msln', 'Pecam1', 'Rgs5', 'Lyz2')) + coord_flip()
VlnPlot(object = PND5.CoDES.EpiMeStro, features = c('Hoxa10', 'Vim', 'Cnn1', 'Tcf21', 'Pkib', 'Pcdh10'))

##match markers in DEGs
PND5.EpiMeStro.markers <- read.csv('CoDES_EpiMeStro_Markers_byCluster_res0.8.csv')
head(PND5.EpiMeStro.markers,1)

PND5.EpiMeStro.markers.tissue <- read.csv('CoDES_EpiMeStro_Markers_byTiss.csv')
head(PND5.EpiMeStro.markers.tissue,1)

PND5.EpiMeStro.markers.tissuetrx <- read.csv('CoDES_EpiMeStro_Markers_byTissTrx.csv')
head(PND5.EpiMeStro.markers.tissuetrx,1)

markers <- c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6', 'Dpt', 'Vcan', 'Col6a3',
             'Hoxa10', 'Vim', 'Cnn1', 'Tcf21', 'Pkib', 'Pcdh10',
             'Epcam', 'Cdh1', 'Krt18')
tissue <- c(replicate(8, 'stroma'), replicate(6, 'mesenchyme'), replicate(3, 'epithelium'))
markerdict <- data.frame(markers, tissue)

PND5.EpiMeStro.markers.match <- merge(x = PND5.EpiMeStro.markers, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.EpiMeStro.markers.match, "CoDES_EpiMeStro_Markers_byCluster_res0.8_matchmarkers.csv", row.names = FALSE)

PND5.EpiMeStro.markers.match.tiss <- merge(x = PND5.EpiMeStro.markers.tissue, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.EpiMeStro.markers.match.tiss, "CoDES_EpiMeStro_Markers_byCluster_matchmarkers_tissue.csv", row.names = FALSE)

PND5.EpiMeStro.markers.match.tisstrx <- merge(x = PND5.EpiMeStro.markers.tissuetrx, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.EpiMeStro.markers.match.tisstrx, "CoDES_EpiMeStro_Markers_byCluster_matchmarkers_tissuetrx.csv", row.names = FALSE)
```
It is not possible to annotate clusters at this resolution while integrated. 
Differences in signaling between exposed/non-exposed are greater than differences from tissue type. This is seen in top DEGs being focused on growth factors, Wnt signaling molecules, etc.

Some clusters are clear:
3 = stroma
4 = epithelium (don't use Smad, likely affected by changed growth signaling)
6 = epithelium 
7 = epithelium
8 = epithelium
13 = stroma
14 = epithelium
15 = epithelium

Some are mixtures of stroma and mesenchyme (0, 2, 5, 7, 9, 11, 12)
And some are indistinguishable with markers from all 3 tissues (1, 10)

Therefore, I will split the seruat by exposure and re-analyze 
```{r Split EpiMeStro Integrated to Control and DES}
PND5.CoDES.EpiMeStro <- readRDS("PND5CoDESEpiMeStro.Rds")

EMSpalette <- c('#F3756E', '#31B34B', '#6F95CD')

DimPlot(PND5.CoDES.EpiMeStro, label = FALSE, repel=TRUE, group.by = "Treatment") & NoAxes() & NoLegend() & ggtitle(NULL)

##reset identities to exposure----
DefaultAssay(PND5.CoDES.EpiMeStro) <- 'SCT'
Idents(PND5.CoDES.EpiMeStro) <- "Treatment"

##subset clusters of interest----
PND5.Co.EpiMeStro <- subset(x = PND5.CoDES.EpiMeStro, idents = 'Control')
PND5.DES.EpiMeStro <- subset(x = PND5.CoDES.EpiMeStro, idents = 'HighDES')
table(PND5.Co.EpiMeStro$Treatment)
  #to check cell number
table(PND5.DES.EpiMeStro$Treatment)
  #to check cell number
  
## standard visualization and clustering----
DefaultAssay(PND5.Co.EpiMeStro) <- "RNA"
PND5.Co.EpiMeStro <- SCTransform(PND5.Co.EpiMeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 
DefaultAssay(PND5.DES.EpiMeStro) <- "RNA"
PND5.DES.EpiMeStro <- SCTransform(PND5.DES.EpiMeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 
  
DepthCor(PND5.Co.EpiMeStro, reduction = 'pca', n = 50)
ElbowPlot(PND5.Co.EpiMeStro, ndims = 50)
DepthCor(PND5.DES.EpiMeStro, reduction = 'pca', n = 50)
ElbowPlot(PND5.DES.EpiMeStro, ndims = 50)

PND5.Co.EpiMeStro <- FindNeighbors(PND5.Co.EpiMeStro, dims = c(1:30)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:30)
PND5.DES.EpiMeStro <- FindNeighbors(PND5.DES.EpiMeStro, dims = c(1:30)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:30)

DimPlot(PND5.Co.EpiMeStro, label = FALSE, repel=TRUE)
DimPlot(PND5.DES.EpiMeStro, label = FALSE, repel=TRUE)

PND5.Co.EpiMeStro[["EpiMeStro_subcluster"]] <- Idents(PND5.Co.EpiMeStro)
PND5.DES.EpiMeStro[["EpiMeStro_subcluster"]] <- Idents(PND5.DES.EpiMeStro)

##check if any cluster has low QC metrics----
Idents(PND5.Co.EpiMeStro) <- 'EpiMeStro_subcluster'
VlnPlot(
  object = PND5.Co.EpiMeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

Idents(PND5.DES.EpiMeStro) <- 'EpiMeStro_subcluster'
VlnPlot(
  object = PND5.DES.EpiMeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

saveRDS(PND5.Co.EpiMeStro, file = "PND5CoEpiMeStro_beforeremovecluster6.Rds", compress = TRUE)
saveRDS(PND5.DES.EpiMeStro, file = "PND5DESEpiMeStro_beforeremovecluster4.Rds", compress = TRUE)

#cluster 6 in control has high RNA (and is impossible to identify - probably junk)
  #remove and recluster
PND5.Co.EpiMeStro <- subset(x = PND5.Co.EpiMeStro, idents = c('1', '2', '3', '4', '5', '7'))
table(PND5.Co.EpiMeStro$EpiMeStro_subcluster)
table(PND5.Co.EpiMeStro$Treatment)

## standard visualization and clustering----
DefaultAssay(PND5.Co.EpiMeStro) <- "RNA"
PND5.Co.EpiMeStro <- SCTransform(PND5.Co.EpiMeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 
  
DepthCor(PND5.Co.EpiMeStro, reduction = 'pca', n = 50)
ElbowPlot(PND5.Co.EpiMeStro, ndims = 50)

PND5.Co.EpiMeStro <- FindNeighbors(PND5.Co.EpiMeStro, dims = c(1:30)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:30)

DimPlot(PND5.Co.EpiMeStro, label = FALSE, repel=TRUE)
DimPlot(PND5.Co.EpiMeStro, label = FALSE, repel=TRUE) & NoAxes() & NoLegend() & ggtitle(NULL)

PND5.Co.EpiMeStro[["EpiMeStro_subcluster"]] <- Idents(PND5.Co.EpiMeStro)

##check if any cluster has low QC metrics----
Idents(PND5.Co.EpiMeStro) <- 'EpiMeStro_subcluster'
VlnPlot(
  object = PND5.Co.EpiMeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

saveRDS(PND5.Co.EpiMeStro, file = "PND5CoEpiMeStro.Rds", compress = TRUE)

#cluster 4 in DES has low RNA and high mitochrondrial DNA
  #remove and recluster
PND5.DES.EpiMeStro <- subset(x = PND5.DES.EpiMeStro, idents = c('1', '2', '3', '5', '6', '7', '8'))
table(PND5.DES.EpiMeStro$EpiMeStro_subcluster)
table(PND5.DES.EpiMeStro$Treatment)

## standard visualization and clustering----
DefaultAssay(PND5.DES.EpiMeStro) <- "RNA"
PND5.DES.EpiMeStro <- SCTransform(PND5.DES.EpiMeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 
  
DepthCor(PND5.DES.EpiMeStro, reduction = 'pca', n = 50)
ElbowPlot(PND5.DES.EpiMeStro, ndims = 50)

PND5.DES.EpiMeStro <- FindNeighbors(PND5.DES.EpiMeStro, dims = c(1:30)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:30)

DimPlot(PND5.DES.EpiMeStro, label = FALSE, repel=TRUE)
DimPlot(PND5.DES.EpiMeStro, label = FALSE, repel=TRUE) & NoAxes() & NoLegend() & ggtitle(NULL)

PND5.DES.EpiMeStro[["EpiMeStro_subcluster"]] <- Idents(PND5.DES.EpiMeStro)

##check if any cluster has low QC metrics----
Idents(PND5.DES.EpiMeStro) <- 'EpiMeStro_subcluster'
VlnPlot(
  object = PND5.DES.EpiMeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)
#cluster 2 also a bit low RNA counts but we will keep for now

saveRDS(PND5.DES.EpiMeStro, file = "PND5DESEpiMeStro.Rds", compress = TRUE)

##annotate----
###find markers per cluster----
####crashes R if using RNA, works if using SCT!-----
options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb
#if crashing, reset environment and redo (quit q() and do not save env)
head(PND5.Co.EpiMeStro,1)
Idents(PND5.Co.EpiMeStro) <- 'SCT_snn_res.0.4'
PND5.Co.EpiMeStro.markers.all <- FindAllMarkers(PND5.Co.EpiMeStro, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.Co.EpiMeStro.markers.all, "PND5CoEpiMeStro_Markers_byCluster.csv", row.names = TRUE)

Idents(PND5.DES.EpiMeStro) <- 'SCT_snn_res.0.4'
PND5.DES.EpiMeStro.markers.all <- FindAllMarkers(PND5.DES.EpiMeStro, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.DES.EpiMeStro.markers.all, "PND5DESEpiMeStro_Markers_byCluster.csv", row.names = TRUE)

PND5.Co.EpiMeStro.markers.all <- dplyr::filter(PND5.Co.EpiMeStro.markers.all,p_val_adj <= 0.05) %>%
    dplyr::filter(, avg_log2FC > 0.5)
PND5.DES.EpiMeStro.markers.all <- dplyr::filter(PND5.DES.EpiMeStro.markers.all,p_val_adj <= 0.05) %>%
    dplyr::filter(, avg_log2FC > 0.5)
  
stroma_markers <- c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6', 'Dpt', 'Vcan', 'Col6a3')
mesenchyme_markers <-  c('Hoxa10', 'Vim', 'Cnn1', 'Tcf21', 'Pkib', 'Pcdh10')
epithelium_markers <- c('Epcam', 'Cdh1', 'Krt18')
myocyte_markers <- c('Mef2a', 'Pdlim3', 'Acta2', 'Chrm3', 'Myh11')
markers <- c(stroma_markers, mesenchyme_markers,  epithelium_markers, myocyte_markers)
tissue <- c(replicate(length(stroma_markers), 'stroma'), replicate(length(mesenchyme_markers), 'mesenchyme'), replicate(length(epithelium_markers), 'epithelium'), replicate(length(myocyte_markers), 'myocyte'))
markerdict <- data.frame(markers, tissue)

PND5.Co.EpiMeStro.markers.match <- merge(x = PND5.Co.EpiMeStro.markers.all, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.Co.EpiMeStro.markers.match, "PND5CoEpiMeStro_Markers_byCluster_matchmarkers.csv", row.names = FALSE)
PND5.DES.EpiMeStro.markers.match <- merge(x = PND5.DES.EpiMeStro.markers.all, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.DES.EpiMeStro.markers.match, "PND5DESEpiMeStro_Markers_byCluster_matchmarkers.csv", row.names = FALSE)

###print plots of cell markers ----
####Control----
DefaultAssay(PND5.Co.EpiMeStro) <- "SCT"
Idents(PND5.Co.EpiMeStro) <- 'SCT_snn_res.0.4'

#####----PND3 markers----
#Vim = Mesenchyme or Stroma (by Spencer paper)
#Hoxa10 = Mesenchyme (uterine)
#Epcam = epithelium
#Krt18 = Epitehlium
#Dcn = stroma
#Vcan = stroma
#Myh11 = myocyte
FeaturePlot(object = PND5.Co.EpiMeStro,
            features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
VlnPlot(object = PND5.Co.EpiMeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'))
DotPlot(object = PND5.Co.EpiMeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11')) + coord_flip()

#make sure no other cell identity is more likely
FeaturePlot(object = PND5.Co.EpiMeStro,
            features = c('Hoxa10', 'Vim', 'Top2a', 'Epcam', 'Msln', 'Pecam1', 'Rgs5', 'Lyz2'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE, keep.scale = 'all')

####DES----
DefaultAssay(PND5.DES.EpiMeStro) <- "SCT"
Idents(PND5.DES.EpiMeStro) <- 'SCT_snn_res.0.4'

FeaturePlot(object = PND5.DES.EpiMeStro,
            features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
VlnPlot(object = PND5.DES.EpiMeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'))
DotPlot(object = PND5.DES.EpiMeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11')) + coord_flip()

#make sure no other cell identity is more likely
FeaturePlot(object = PND5.DES.EpiMeStro,
            features = c('Hoxa10', 'Vim', 'Top2a', 'Epcam', 'Msln', 'Pecam1', 'Rgs5', 'Lyz2'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE, keep.scale = 'all')

```

cluster 6 in control has stroma, myocyte, and epithelial markers, no statistically significant GO terms of its top 20 DEGs, and higher total_nCountRNA than the rest. Junk? or something else? 
otherwise, rest are clear by markers. (ended up removing and reclustering)

Pdlim3 may be confusing marker - also seen in cardiac fibroblasts so could be in stromal cells
AND 'PDLIM3 expression was demonstrated to increase in intestinal epithelial cells' (https://doi.org/10.1093/ecco-jcc/jjad212.1361)
but Cdh1 still seems uniquely epithelium and Dpt uniquely stromal (this is the junk I was seeing before in my integrated)

Vcan found in perivascular eMSC (https://doi.org/10.1095/biolreprod.111.095885); medium expression in human endometrial smooth muscle cells (https://www.proteinatlas.org/ENSG00000038427-VCAN/tissue/endometrium)

```{r Annotate Unintegrated EpiMeStro Control and DES}
PND5.Co.EpiMeStro <- readRDS("PND5CoEpiMeStro.Rds")
PND5.DES.EpiMeStro <- readRDS("PND5DESEpiMeStro.Rds")

PND5.Co.EpiMeStro <- RenameIdents(object = PND5.Co.EpiMeStro,
                              '0' = 'Epithelium',
                              '1' = 'Epithelium',
                              '2' = 'Mesenchyme/Myocyte', 
                              '3' = 'Mesenchyme/Stroma',
                              '4' = 'Mesenchyme/Stroma',
                              '5' = 'Epithelium')
PND5.Co.EpiMeStro[['EpiMeStro_Annotated']] <- Idents(PND5.Co.EpiMeStro)

EMSpalette <- c('#F3756E', '#6F95CD', '#31B34B')

DimPlot(PND5.Co.EpiMeStro, label = FALSE, repel=TRUE)
DimPlot(PND5.Co.EpiMeStro, label = FALSE, repel=TRUE, cols = EMSpalette) & NoAxes() & NoLegend() & ggtitle(NULL)

PND5.DES.EpiMeStro <- RenameIdents(object = PND5.DES.EpiMeStro,
                                   '0' = 'Epithelium',
                                   '1' = 'Mesenchyme/Stroma',
                                   '2' = 'Epithelium',
                                   '3' = 'Mesenchyme/Myocyte',
                                   '4' = 'Epithelium',
                                   '5' = 'Mesenchyme/Stroma',
                                   '6' = 'Epithelium',
                                   '7' = 'Mesenchyme/Myocyte')
PND5.DES.EpiMeStro[['EpiMeStro_Annotated']] <- Idents(PND5.DES.EpiMeStro)

DimPlot(PND5.DES.EpiMeStro, label = FALSE, repel=TRUE)
DimPlot(PND5.DES.EpiMeStro, label = FALSE, repel=TRUE) & NoAxes() & NoLegend() & ggtitle(NULL)

saveRDS(PND5.Co.EpiMeStro, file = "PND5CoEpiMeStro_Annotated.Rds", compress = TRUE)
saveRDS(PND5.DES.EpiMeStro, file = "PND5DESEpiMeStro_Annotated.Rds", compress = TRUE)

###find markers per tissue----
options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb

####control----
Idents(PND5.Co.EpiMeStro) <- 'EpiMeStro_Annotated'
PND5.Co.EpiMeStro.markers.tiss <- FindAllMarkers(PND5.Co.EpiMeStro, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.Co.EpiMeStro.markers.tiss, "PND5CoEpiMeStro_Markers_byTissue.csv", row.names = TRUE)

####DES----
Idents(PND5.DES.EpiMeStro) <- 'EpiMeStro_Annotated'

PND5.DES.EpiMeStro.markers.tiss <- FindAllMarkers(PND5.DES.EpiMeStro, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.DES.EpiMeStro.markers.tiss, "PND5DESEpiMeStro_Markers_byTissue.csv", row.names = TRUE)

###compare tissues between exposures-----
####first merge annotated Control and DES epimestro----
PND5.CoDES.EpiMeStro.merge <- merge(x = PND5.Co.EpiMeStro, y = PND5.DES.EpiMeStro)
saveRDS(PND5.CoDES.EpiMeStro.merge, file = "PND5CoDESEpiMeStro_merged.Rds", compress = TRUE)

####next split out tissue of interest and find markers----
#####Epithelium----
Idents(PND5.CoDES.EpiMeStro.merge) <- 'EpiMeStro_Annotated'

PND5.CoDES.Epi.merge <- subset(PND5.CoDES.EpiMeStro.merge, idents = 'Epithelium')
Idents(PND5.CoDES.Epi.merge) <- 'Treatment'
PND5.CoDES.Epi.merge <- PrepSCTFindMarkers(PND5.CoDES.Epi.merge)
PND5.CoDES.Epi.markers <- FindMarkers(PND5.CoDES.Epi.merge, ident.1 = 'HighDES',
                                      ident.2 = 'Control', test.use = 'MAST', 
                                      only.pos = FALSE, assay = 'SCT', 
                                      latent.vars='percent.mt')
write.csv(PND5.CoDES.Epi.markers, "PND5CoDESEpiMeStro_merge_Markers_Epithelium.csv", row.names = TRUE)

#attach average expression
Idents(PND5.CoDES.Epi.merge) <- 'Treatment'
PND5.CoDES.Epi.markerexp <- AverageExpression(PND5.CoDES.Epi.merge,
                                              features = PND5.CoDES.Epi.markers$X,
                                              return.seurat = FALSE,
                                              group.by = 'ident', assay = 'SCT',
                                              layer = 'scale.data')
head(PND5.CoDES.Epi.markerexp$SCT,1)
PND5.CoDES.Epi.markerexp <- as.data.frame(PND5.CoDES.Epi.markerexp$SCT)
PND5.CoDES.Epi.markerexp$X <- rownames(PND5.CoDES.Epi.markerexp)
head(PND5.CoDES.Epi.markerexp)

PND5.CoDES.Epi.markers.exp <- merge(PND5.CoDES.Epi.markers, PND5.CoDES.Epi.markerexp,
                                    by = 'X')
head(PND5.CoDES.Epi.markers.exp)
write.csv(PND5.CoDES.Epi.markers.exp, "PND5CoDESEpiMeStro_merge_MarkerswAvgExp_Epithelium.csv", row.names = TRUE)

#####DEX volcano----
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
EnhancedVolcano(PND5.CoDES.Epi.markers, lab = PND5.CoDES.Epi.markers$X,
                x = 'avg_log2FC', y = 'p_val_adj', title = 'Control vs DES',
                pCutoff = 0.01, FCcutoff = 1.5, ylim = c(0, -log10(10e-300)), 
                selectLab = PND5.CoDES.Epi.markers.updn$X) + coord_flip()
```


```{r EpiMeStro Unintegrated - Cell-Cell Communication with Liana}

#if(!require("circlize")){
#  install.packages("circlize", quiet = TRUE,
#                   repos = "http://cran.us.r-project.org")}

#CCC----
##CCI with liana----
#using tutorial from https://saezlab.github.io/liana/articles/liana_tutorial.html

#library(BiocManager)
#remotes::install_github('saezlab/liana')
#library(liana)

##set liana to use mouse genes----
#using tutorial at https://rdrr.io/github/saezlab/liana/f/vignettes/liana_ortho.Rmd

#convert LIANA's Consensus resource to murine symbols
op_resource <- select_resource('Consensus')[[1]]

#generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) #mouse
#can then run LIANA with orthologous resource using:
#resource = 'custom', # resource has to be set to 'custom' to work with external resources
#external_resource = ortholog_resource, # provide orthologous resource

###multiple CCC with liana----
PND5.Co.EpiMeStro <- readRDS("PND5CoEpiMeStro_Annotated.Rds")
PND5.DES.EpiMeStro <- readRDS("PND5DESEpiMeStro_Annotated.Rds")

Idents(PND5.Co.EpiMeStro) <- 'EpiMeStro_Annotated'
Idents(PND5.DES.EpiMeStro) <- 'EpiMeStro_Annotated'

#leave out myocyte
PND5.Co.EpiMeStro <- subset(PND5.Co.EpiMeStro, idents = c('Epithelium', 'Mesenchyme/Stroma'))
PND5.DES.EpiMeStro <- subset(PND5.DES.EpiMeStro, idents = c('Epithelium', 'Mesenchyme/Stroma'))

DefaultAssay(PND5.Co.EpiMeStro) <- 'SCT'
DefaultAssay(PND5.DES.EpiMeStro) <- 'SCT'

PND5.EpMeStro.C.CCC <- liana_wrap(PND5.Co.EpiMeStro,
                                  idents_col = 'EpiMeStro_Annotated',
                                  resource = 'custom',
                                  external_resource = ortholog_resource)
PND5.EpMeStro.HD.CCC <- liana_wrap(PND5.DES.EpiMeStro,
                                   idents_col = 'EpiMeStro_Annotated',
                                   resource = 'custom',
                                   external_resource = ortholog_resource)

PND5.EpMeStro.C.CCC
PND5.EpMeStro.HD.CCC

PND5.EpMeStro.C.CCC.aggregate <- liana_aggregate(PND5.EpMeStro.C.CCC)
PND5.EpMeStro.HD.CCC.aggregate <- liana_aggregate(PND5.EpMeStro.HD.CCC)

saveRDS(PND5.EpMeStro.C.CCC, file = "PND5CoDES_EpiMeStro_liana_Con.rds")
saveRDS(PND5.EpMeStro.C.CCC.aggregate, file = "PND5CoDES_EpiMeStro_liana_agg_Con.rds")
saveRDS(PND5.EpMeStro.HD.CCC, file = "PND5CoDES_EpiMeStro_liana_HD.rds")
saveRDS(PND5.EpMeStro.HD.CCC.aggregate, file = "PND5CoDES_EpiMeStro_liana_agg_HD.rds")

##visualize aggregate scores----
liana_dotplot(PND5.EpMeStro.C.CCC.aggregate, ntop = 20) +
  theme(
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    strip.text = element_text(size = 12) 
  )
liana_dotplot(PND5.EpMeStro.C.CCC.aggregate, ntop = 50) +
  theme(
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    strip.text = element_text(size = 12) 
  )

###save table----
write.csv(PND5.EpMeStro.C.CCC.aggregate, "PND5CoDES_EpiMeStro_liana_aggregate_Con.csv", row.names = TRUE)
write.csv(PND5.EpMeStro.HD.CCC.aggregate, "PND5CoDES_EpiMeStro_liana_aggregate_HD.csv", row.names = TRUE)

#annotate aggregate score table to look for patterns----
PND5.EpMeStro.C.CCC.aggregate <- read.csv("PND5CoDES_EpiMeStro_liana_aggregate_Con.csv")
PND5.EpMeStro.HD.CCC.aggregate <- read.csv("PND5CoDES_EpiMeStro_liana_aggregate_HD.csv")
##merge with manually created ligand:receptor annotations
LR_Dictionary <- read.csv("ligandreceptor_dictionary.csv")
  LR_Dictionary <- LR_Dictionary[,2:4]

head(PND5.EpMeStro.C.CCC.aggregate,1)
head(LR_Dictionary,1)

PND5.EpMeStro.C.CCC.aggregate.desc <- left_join(x = PND5.EpMeStro.C.CCC.aggregate, y=LR_Dictionary)
PND5.EpMeStro.HD.CCC.aggregate.desc <- left_join(x = PND5.EpMeStro.HD.CCC.aggregate, y=LR_Dictionary)

head(PND5.EpMeStro.C.CCC.aggregate.desc)
head(PND5.EpMeStro.HD.CCC.aggregate.desc)

###save table----
write.csv(PND5.EpMeStro.C.CCC.aggregate.desc, "PND5CoDES_EpiMeStro_liana_aggregate_Con_descriptive.csv", row.names = TRUE)
write.csv(PND5.EpMeStro.HD.CCC.aggregate.desc, "PND5CoDES_EpiMeStro_liana_aggregate_HD_descriptive.csv", row.names = TRUE)

###make chord diagram to show frequency of interactions----
PND5.EpMeStro.C.CCC.aggregate <- dplyr::filter(PND5.EpMeStro.C.CCC.aggregate, aggregate_rank <= 0.01)
PND5.EpMeStro.HD.CCC.aggregate <- dplyr::filter(PND5.EpMeStro.HD.CCC.aggregate, aggregate_rank <= 0.01)

chord_freq(PND5.EpMeStro.C.CCC.aggregate,
           source_groups = c("Epithelium", "Mesenchyme/Stroma"),
           target_groups = c("Epithelium", "Mesenchyme/Stroma"))
chord_freq(PND5.EpMeStro.HD.CCC.aggregate,
           source_groups = c("Epithelium", "Mesenchyme/Stroma"),
           target_groups = c("Epithelium", "Mesenchyme/Stroma"))
```

Now to do analysis based on ATAC

I'm not sure how to compare peaks between control and DES without them being in an integrated Seurat together. Could either merge and integrate or put markers on integrated cells (which the unintegrated are taken from) - doing fewer transformations seems better, so put annotation labels on integrated EpiMeStro
Might also be better place to make violins between control and des (like for wnts)

```{r unintegrated labels on integrated EpiMeStro}

PND5.Co.EpiMeStro <- readRDS("PND5CoEpiMeStro_Annotated.Rds")
PND5.DES.EpiMeStro <- readRDS("PND5DESEpiMeStro_Annotated.Rds")

PND5CoUnintLabel <- PND5.Co.EpiMeStro[['EpiMeStro_Annotated']]
PND5DESUnintLabel <- PND5.DES.EpiMeStro[['EpiMeStro_Annotated']]
head(PND5CoUnintLabel)
head(PND5DESUnintLabel)

PNDCODESUnintLabel <- bind_rows(PND5CoUnintLabel, PND5DESUnintLabel)
head(PNDCODESUnintLabel)

PND5.CoDES.EpiMeStro <- readRDS("PND5CoDESEpiMeStro.Rds")
PND5.CoDES.EpiMeStro[['EpiMeStro_Annotated']] <- PNDCODESUnintLabel
table(PND5.CoDES.EpiMeStro$EpiMeStro_Annotated)

Idents(PND5.CoDES.EpiMeStro) <- 'EpiMeStro_Annotated'

DimPlot(PND5.CoDES.EpiMeStro, label = TRUE, repel=TRUE, reduction = 'umap')

saveRDS(PND5.CoDES.EpiMeStro, file = "PND5CoDESEpiMeStro_annotatedfromUnintegrated_beforedrop.Rds", compress = TRUE)

PND5.CoDES.EpiMeStro <- subset(PND5.CoDES.EpiMeStro, idents = c('Epithelium', 'Mesenchyme/Stroma', 'Mesenchyme/Myocyte'))

DimPlot(PND5.CoDES.EpiMeStro, label = TRUE, repel=TRUE, reduction = 'umap')

saveRDS(PND5.CoDES.EpiMeStro, file = "PND5CoDESEpiMeStro_annotatedfromUnintegrated.Rds", compress = TRUE)

```

```{r reintegrated EpiMeStro - integrate ATAC}
#UMAP based on ATAC----
#from https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html#clustering-and-further-dimensionality-reduction

PND5.CoDES.EpiMeStro <- readRDS("PND5CoDESEpiMeStro_annotatedfromUnintegrated.Rds")

DefaultAssay(PND5.CoDES.EpiMeStro) <- "peaks"
PND5.CoDES.EpiMeStro <- FindTopFeatures(PND5.CoDES.EpiMeStro, min.cutoff = 'q0')
PND5.CoDES.EpiMeStro <- RunTFIDF(PND5.CoDES.EpiMeStro)
PND5.CoDES.EpiMeStro <- RunSVD(PND5.CoDES.EpiMeStro)

DepthCor(PND5.CoDES.EpiMeStro, n = 30)
  #first 3 correllated
ElbowPlot(PND5.CoDES.EpiMeStro, ndims = 50)

PND5.CoDES.EpiMeStro <- RunUMAP(PND5.CoDES.EpiMeStro, reduction = "lsi", dims = 4:30)
PND5.CoDES.EpiMeStro <- FindNeighbors(PND5.CoDES.EpiMeStro, reduction = 'lsi', dims = 4:30)
#kparam may need to adjust - just left to default
PND5.CoDES.EpiMeStro <- FindClusters(PND5.CoDES.EpiMeStro)

DimPlot(PND5.CoDES.EpiMeStro, label = TRUE, repel=TRUE, reduction = 'umap')
DimPlot(PND5.CoDES.EpiMeStro, label = TRUE, repel=TRUE, group.by = 'Treatment', reduction = 'umap')

saveRDS(PND5.CoDES.EpiMeStro, file = "PND5CoEpiMeStro_ATACUMAP.Rds", compress = TRUE)

##find overlapping peaks----
####per treatment----
Idents(PND5.CoDES.EpiMeStro) <- 'Treatment'
DefaultAssay(PND5.CoDES.EpiMeStro) <- 'peaks'

PND5.CO.APs <- AccessiblePeaks(
  PND5.CoDES.EpiMeStro,
  idents = 'Control'
)
PND5.HD.APs <- AccessiblePeaks(
  PND5.CoDES.EpiMeStro,
  idents = 'HighDES'
)

PND5.CODES.shared.APs <-  PND5.CO.APs%in%PND5.HD.APs

length(PND5.CODES.shared.APs[PND5.CODES.shared.APs == TRUE])
length(PND5.CO.APs)
length(PND5.HD.APs)

venn.diagram(
  x = list(PND5.CO.APs, PND5.HD.APs),
  category.names = c('Control', 'DES'),
  filename = "PND5HD_CODES_Epi_Venn_ATACPeaks.png",
  output = TRUE
)

##Motif analysis----
###from https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html#clustering-and-further-dimensionality-reduction
###control----
#get list of motif position matrices from JASPAR
PND5CoDES.pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(species = "Mus musculus", all_versions = FALSE))
#scan the DNA seq of each peak for the presence of each motif
motif.matrix.PND5CoDES <- CreateMotifMatrix(
  features = granges(PND5.CoDES.EpiMeStro$peaks),
  pwm = PND5CoDES.pfm,
  genome = 'mm10',
  use.counts = FALSE)
dim(motif.matrix.PND5CoDES)
as.matrix(motif.matrix.PND5CoDES[1:10, 1:10])
#create new Motif object to store the results
motif.PND5CoDES <- CreateMotifObject(
  data = motif.matrix.PND5CoDES,
  pwm = PND5CoDES.pfm)
#Add the Motif object to the assay
PND5.CoDES.EpiMeStro <- SetAssayData(
  object = PND5.CoDES.EpiMeStro,
  assay = 'peaks',
  layer = 'motifs',
  new.data = motif.PND5CoDES)
PND5.CoDES.EpiMeStro$peaks@motifs

#calculate stats for each peak:
  #GC content, length, dinulceotide frequencies, etc.
  #these are used in test of motif enrichment
PND5.CoDES.EpiMeStro$peaks <- RegionStats(object = PND5.CoDES.EpiMeStro$peaks, genome = BSgenome.Mmusculus.UCSC.mm10)

#find differential accessible peaks in each cluster
DA.PND5CoDES <- FindAllMarkers(
  object = PND5.CoDES.EpiMeStro,
  only.pos = TRUE,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments')
write.csv(DA.PND5CoDES, 'AllDApeaksPND5CoDES.csv', row.names = TRUE)

head(DA.PND5CoDES)

#get top DAs with lowest p-values
top.DA.peak.PND5CoDES <- rownames(DA.PND5CoDES[DA.PND5CoDES$p_val < 0.01, ])
head(top.DA.peak.PND5CoDES)
write.csv(top.DA.peak.PND5CoDES, 'topDApeaksPND5CoDES.csv', row.names = TRUE)

#find motifs enriched in top peaks
enriched.motifs.PND5CoDES <- FindMotifs(object = PND5.CoDES.EpiMeStro, features = top.DA.peak.PND5CoDES)
head(enriched.motifs.PND5CoDES)
write.csv(enriched.motifs.PND5CoDES, 'EnrichedMotifsPND5CoDES.csv', row.names = TRUE)
#visualize these motifs
MotifPlot(object = PND5.CoDES.EpiMeStro, motifs = head(rownames(enriched.motifs.PND5CoDES)))

#compute per-cell motif activity score
PND5.CoDES.EpiMeStro <- RunChromVAR(
  object = PND5.CoDES.EpiMeStro,
  genome = BSgenome.Mmusculus.UCSC.mm10)
PND5.CoDES.EpiMeStro$chromvar
GetAssayData(PND5.CoDES.EpiMeStro$chromvar)[1:10, 1:3]

#check specific motifs
DefaultAssay(PND5.CoDES.EpiMeStro) <- "chromvar"
FeaturePlot(
  object = PND5.CoDES.EpiMeStro,
  features = 'MA1684.1',
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 1)

####save RDS with motif scores----
saveRDS(PND5.CoDES.EpiMeStro, file = "PND5CoDESEpiMeStro_ATAC_chromvar.Rds", compress = TRUE)
###find differentially accessible peaks----
####per treatment----
Idents(PND5.CoDES.EpiMeStro) <- 'Treatment'
DefaultAssay(PND5.CoDES.EpiMeStro) <- 'peaks'
#lets try dropping the min% from 0.2 to 0.01
#test that used to take 1 min now takes 1hr
PND5CoDES.DApeaks.HDtoC <- FindMarkers(
  object = PND5.CoDES.EpiMeStro,
  ident.1 = 'HighDES',
  ident.2 = 'Control',
  only.pos = FALSE,
  min.pct = 0.01,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments',
  assay = 'peaks',
  verbose = TRUE
)

write.csv(PND5CoDES.DApeaks.HDtoC, "PND5CoDESEpiMeStro_DApeaks_HDtoC.csv", row.names = TRUE)

#find peaks with sig p-values
PND5CoDES.top.DApeak.HDtoC <- PND5CoDES.DApeaks.HDtoC[PND5CoDES.DApeaks.HDtoC$p_val < 0.05, ]
PND5CoDES.top.DApeak.HDtoC <- PND5CoDES.top.DApeak.HDtoC$X

head(PND5CoDES.DApeaks.HDtoC)

#find motifs enriched in top sig DAs
PND5CoDES.enrichedmotifs.HDtoC <- FindMotifs(
  object <- PND5.CoDES.EpiMeStro,
  features = PND5CoDES.top.DApeak.HDtoC,
  assay = 'peaks'
)
head(PND5CoDES.enrichedmotifs.HDtoC)

write.csv(PND5CoDES.enrichedmotifs.HDtoC, "PND5CoDESEpiMeStro_enrichedmotifs_HDtoC.csv", row.names = TRUE)

MotifPlot(
  object = PND5.CoDES.EpiMeStro,
  motifs = head(rownames(PND5CoDES.enrichedmotifs.HDtoC))
)

###annotate regions----
#find genes near or at genomic regions
PND5CoDES.DApeaks_open_HDtoC <- PND5CoDES.DApeaks.HDtoC %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC >= 1.5 | avg_log2FC <= -1.5) %>%
  arrange(avg_log2FC)
PND5CoDES.closestgene_HDtoC <- ClosestFeature(PND5.CoDES.EpiMeStro$peaks, regions = PND5CoDES.DApeaks_open_HDtoC$X)

write.csv(PND5CoDES.closestgene_HDtoC, "PND5CoDESEpiMeStro_DAclosestgene_HDtoC.csv", row.names = TRUE)

#check if working using a known diff accessible gene
PND5CoDES.closestgene_HDtoC[PND5CoDES.closestgene_HDtoC$gene_name == 'Padi1', ]
  
DefaultAssay(PND5.CoDES.EpiMeStro) <- 'peaks'
Idents(PND5.CoDES.EpiMeStro) <- 'EpiMeStro_Annotated'
  table(PND5.CoDES.EpiMeStro$EpiMeStro_Annotated)
  
#####combine tables----
PND5CoDES.DApeaks.HDtoC <- read.csv("PND5CoDESEpiMeStro_DApeaks_HDtoC.csv")

PND5CoDES.closestgene_peakinfo_HDtoC <- merge(PND5CoDES.DApeaks.HDtoC, PND5CoDES.closestgene_HDtoC, by.x = 'X', by.y = 'query_region')
head(PND5CoDES.closestgene_peakinfo_HDtoC)
  write.csv(PND5CoDES.closestgene_peakinfo_HDtoC, "PND5CoDESEpiMeStro_DAclosestgene_withpeakinfo_HDtoC.csv", row.names = TRUE)
  
###differential activity----
#calculate a per-cell motif activity score
PND5.CoDES.EpiMeStro <- RunChromVAR(
  object = PND5.CoDES.EpiMeStro,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

PND5.CoDES.EpiMeStro$chromvar
GetAssayData(PND5.CoDES.EpiMeStro$chromvar)[1:10,1:3]
DefaultAssay(PND5.CoDES.EpiMeStro) <- 'chromvar'
  FeaturePlot(
    object = PND5.CoDES.EpiMeStro,
    features = 'MA0489.2'
  )

  FeaturePlot(
    object = PND5.CoDES.EpiMeStro,
    features = 'MA0112.1'
  )

saveRDS(PND5.CoDES.EpiMeStro, file = "PND5CoDESEpiMeStro_ATAC_ChromVar.rds")

DefaultAssay(PND5.CoDES.EpiMeStro) <- 'chromvar'
Idents(PND5.CoDES.EpiMeStro) <- 'Treatment'

PND5CoDES.ATAC.differential.activity <- FindMarkers(
  object = PND5.CoDES.EpiMeStro,
  ident.1 = 'HighDES',
  ident.2 = 'Control',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = 'avg_diff',
  verbose = TRUE
)
  
MotifPlot(
  object = PND5.CoDES.EpiMeStro,
  motifs = head(rownames(PND5CoDES.ATAC.differential.activity)),
  assay = 'peaks'
)

write.csv(PND5CoDES.ATAC.differential.activity, "PND5CoDESEpiMeStro_DiffAct_HDtoC.csv", row.names = TRUE)

```

Now going to do separate analyses on Epithelium and MesStro

First, epithelium
break out epithelium to look for progenitor markers, basal markers, etc

```{r Subset Epithelium from reintegrated EpiMeStro}
PND5.CoDES.EpiMeStro <- readRDS("PND5CoDESEpiMeStro_annotatedfromUnintegrated.Rds")
Idents(PND5.CoDES.EpiMeStro) <- 'EpiMeStro_Annotated'
PND5.CoDES.Epith <- subset(PND5.CoDES.EpiMeStro, idents = 'Epithelium')
table(PND5.CoDES.Epith$EpiMeStro_Annotated)
ttable(PND5.CoDES.Epith$Treatment)

## standard visualization and clustering----
DefaultAssay(PND5.CoDES.Epith) <- "RNA"
PND5.CoDES.Epith <- SCTransform(PND5.CoDES.Epith) %>% 
  RunPCA(features = VariableFeatures(object = .)) 

DepthCor(PND5.CoDES.Epith, reduction = 'pca', n = 50)
  #first 2 correlate
ElbowPlot(PND5.CoDES.Epith, ndims = 50)

PND5.CoDES.Epith <- FindNeighbors(PND5.CoDES.Epith, dims = c(3:30)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 3:30)

DimPlot(PND5.CoDES.Epith, label = FALSE, repel=TRUE)
DimPlot(PND5.CoDES.Epith, label = FALSE, repel=TRUE, group.by="Treatment")

PND5.CoDES.Epith[["Epithelium_subcluster"]] <- Idents(PND5.CoDES.Epith)

##check if any cluster has low QC metrics----
Idents(PND5.CoDES.Epith) <- 'Epithelium_subcluster'
VlnPlot(
  object = PND5.CoDES.Epith,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)
  #5, 8, 9 look like junk. may need to remove if not unique

saveRDS(PND5.CoDES.Epith, file = "PND5CoDESEpith.Rds", compress = TRUE)

####DEGs----
DefaultAssay(PND5.CoDES.Epith) <- "SCT"
Idents(PND5.CoDES.Epith) <- 'Treatment'

PND5.CoDES.Epith <- PrepSCTFindMarkers(PND5.CoDES.Epith)

PND5.CODES.Epith.markers.all <- FindMarkers(PND5.CoDES.Epith, test.use = 'MAST',
                                            only.pos = FALSE, ident.1 = 'HighDES',
                                            ident.2 = 'Control',
                           assay = 'SCT', latent.vars='percent.mt')

write.csv(PND5.CODES.Epith.markers.all, "PND5CoDESEpith.DEG.DEStoCon.csv")
  
#####----Select Marker Genes----
DefaultAssay(PND5.CoDES.Epith) <- "SCT"
Idents(PND5.CoDES.Epith) <- 'Epithelium_subcluster'

FeaturePlot(object = PND5.CoDES.Epith,
            features = c('Mcm2', 'Mcm5', 'Mki67', 'Pcna', 'Top2a'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
VlnPlot(object = PND5.CoDES.Epith, features = c('Mcm2', 'Mcm5', 'Mki67', 'Pcna', 'Top2a'))
DotPlot(object = PND5.CoDES.Epith, features = c('Mcm2', 'Mcm5', 'Mki67', 'Pcna', 'Top2a')) + coord_flip()

FeaturePlot(object = PND5.CoDES.Epith,
            features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
DotPlot(object = PND5.CoDES.Epith, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11')) + coord_flip()
```

Clusters are separating entirely by treatment. I will need to subset and analyze DES and control separately.

```{r Epithelium Unintegrated}
PND5.CoDES.Epith <- readRDS("Epithelium/Unintegrated/PND5CoDESEpith.Rds")
Idents(PND5.CoDES.Epith) <- 'Treatment'
PND5.Co.Epith <- subset(PND5.CoDES.Epith, idents = 'Control')
PND5.DES.Epith <- subset(PND5.CoDES.Epith, idents = 'HighDES')
table(PND5.Co.Epith$EpiMeStro_Annotated)
table(PND5.DES.Epith$EpiMeStro_Annotated)
table(PND5.Co.Epith$Treatment)
table(PND5.DES.Epith$Treatment)

## standard visualization and clustering----
DefaultAssay(PND5.Co.Epith) <- "RNA"
PND5.Co.Epith <- SCTransform(PND5.Co.Epith) %>% 
  RunPCA(features = VariableFeatures(object = .)) 
DefaultAssay(PND5.DES.Epith) <- "RNA"
PND5.DES.Epith <- SCTransform(PND5.DES.Epith) %>% 
  RunPCA(features = VariableFeatures(object = .)) 

DepthCor(PND5.Co.Epith, reduction = 'pca', n = 50)
ElbowPlot(PND5.Co.Epith, ndims = 50)
DepthCor(PND5.DES.Epith, reduction = 'pca', n = 50)
ElbowPlot(PND5.DES.Epith, ndims = 50)

PND5.Co.Epith <- FindNeighbors(PND5.Co.Epith, dims = c(6:25)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 6:25)
PND5.DES.Epith <- FindNeighbors(PND5.DES.Epith, dims = c(4:25)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 4:25)

DimPlot(PND5.Co.Epith, label = FALSE, repel=TRUE)
DimPlot(PND5.DES.Epith, label = FALSE, repel=TRUE)
DimPlot(PND5.Co.Epith, label = FALSE, repel=TRUE) & NoAxes() & NoLegend() & ggtitle(NULL)
DimPlot(PND5.DES.Epith, label = FALSE, repel=TRUE) & NoAxes() & NoLegend() & ggtitle(NULL)

PND5.Co.Epith[["Epithelium_subcluster"]] <- Idents(PND5.Co.Epith)
PND5.DES.Epith[["Epithelium_subcluster"]] <- Idents(PND5.DES.Epith)


##check if any cluster has low QC metrics----
Idents(PND5.Co.Epith) <- 'Epithelium_subcluster'
Idents(PND5.DES.Epith) <- 'Epithelium_subcluster'
VlnPlot(
  object = PND5.Co.Epith,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)
VlnPlot(
  object = PND5.DES.Epith,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

saveRDS(PND5.Co.Epith, file = "PND5CoEpith.Rds", compress = TRUE)
saveRDS(PND5.DES.Epith, file = "PND5DESEpith.Rds", compress = TRUE)

##annotate----
###find markers per cluster----
####crashes R if using RNA, works if using SCT!-----
options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb
#if crashing, reset environment and redo (quit q() and do not save env)
head(PND5.Co.Epith,1)
Idents(PND5.Co.Epith) <- 'Epithelium_subcluster'
Idents(PND5.DES.Epith) <- 'Epithelium_subcluster'

###Control----
PND5.Co.Epith.markers.all <- FindAllMarkers(PND5.Co.Epith, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.Co.Epith.markers.all, "PND5CoEpithMarkers_byCluster.csv", row.names = TRUE)
  PND5.Co.Epith.markers.all <- read.csv("PND5CoEpithMarkers_byCluster.csv")

PND5.Co.Epith.markers.all <- dplyr::filter(PND5.Co.Epith.markers.all,p_val_adj <= 0.05) %>%
    dplyr::filter(, avg_log2FC > 0.5)

stroma_markers <- c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6', 'Dpt', 'Vcan', 'Col6a3')
mesenchyme_markers <-  c('Hoxa10', 'Vim', 'Cnn1', 'Tcf21', 'Pkib', 'Pcdh10')
epithelium_markers <- c('Epcam', 'Cdh1', 'Krt18')
myocyte_markers <- c('Mef2a', 'Pdlim3', 'Acta2', 'Chrm3', 'Myh11')
luminal_markers <- c('Calb1', 'Cited4', 'Krt8', 'Krt18') 
basal_markers <- c('Trp63', 'Krt5', 'Krt14', 'Bcl11b', 'Sox15')
gland_markers <- c('Foxa2', 'Cxcl15', 'Pcna', 'Aldh1a1')
markers <- c(stroma_markers, mesenchyme_markers,  epithelium_markers, 
             myocyte_markers, luminal_markers, basal_markers, gland_markers)
tissue <- c(replicate(length(stroma_markers), 'stroma'), 
            replicate(length(mesenchyme_markers), 'mesenchyme'), 
            replicate(length(epithelium_markers), 'epithelium'), 
            replicate(length(myocyte_markers), 'myocyte'),
            replicate(length(luminal_markers), 'luminal epithelium'),
            replicate(length(basal_markers), 'basal epithelium'),
            replicate(length(gland_markers), 'glandular epithelium'))
markerdict <- data.frame(markers, tissue)

PND5.Co.Epith.markers.match <- merge(x = PND5.Co.Epith.markers.all, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.Co.Epith.markers.match, "PND5CoEpith_Markers_byCluster_matchmarkers.csv", row.names = FALSE)

###DES----
PND5.DES.Epith.markers.all <- FindAllMarkers(PND5.DES.Epith, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.DES.Epith.markers.all, "PND5DESEpithMarkers_byCluster.csv", row.names = TRUE)

PND5.DES.Epith.markers.all <- dplyr::filter(PND5.DES.Epith.markers.all,p_val_adj <= 0.05) %>%
    dplyr::filter(, avg_log2FC > 0.5)

stroma_markers <- c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6', 'Dpt', 'Vcan', 'Col6a3')
mesenchyme_markers <-  c('Hoxa10', 'Vim', 'Cnn1', 'Tcf21', 'Pkib', 'Pcdh10')
epithelium_markers <- c('Epcam', 'Cdh1', 'Krt18')
myocyte_markers <- c('Mef2a', 'Pdlim3', 'Acta2', 'Chrm3', 'Myh11')
luminal_markers <- c('Calb1', 'Cited4', 'Krt8', 'Krt18') 
basal_markers <- c('Trp63', 'Krt5', 'Krt14', 'Bcl11b', 'Sox15')
gland_markers <- c('Foxa2', 'Cxcl15', 'Pcna', 'Aldh1a1')
markers <- c(stroma_markers, mesenchyme_markers,  epithelium_markers, 
             myocyte_markers, luminal_markers, basal_markers, gland_markers)
tissue <- c(replicate(length(stroma_markers), 'stroma'), 
            replicate(length(mesenchyme_markers), 'mesenchyme'), 
            replicate(length(epithelium_markers), 'epithelium'), 
            replicate(length(myocyte_markers), 'myocyte'),
            replicate(length(luminal_markers), 'luminal epithelium'),
            replicate(length(basal_markers), 'basal epithelium'),
            replicate(length(gland_markers), 'glandular epithelium'))
markerdict <- data.frame(markers, tissue)

PND5.DES.Epith.markers.match <- merge(x = PND5.DES.Epith.markers.all, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.DES.Epith.markers.match, "PND5DESEpith_Markers_byCluster_matchmarkers.csv", row.names = FALSE)

#merge co and DES for shared visualization----
DefaultAssay(PND5.Co.Epith) <- "SCT"
Idents(PND5.Co.Epith) <- 'Epithelium_subcluster'

DefaultAssay(PND5.DES.Epith) <- "SCT"
Idents(PND5.DES.Epith) <- 'Epithelium_subcluster'

PND5.merge.Epith <- merge(PND5.Co.Epith, PND5.DES.Epith, merge.data = TRUE, merge.dr = TRUE)
head(PND5.merge.Epith,1)
PND5.merge.Epith$SCT

saveRDS(PND5.merge.Epith, 'PND5mergedEpith.Rds')

#use to make violin, feature, & dot plots comparing gene levels between treatments

```

```{r Epithelium ATAC - Integrated for peak comparison between Treatments}
#UMAP based on ATAC----
#from https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html#clustering-and-further-dimensionality-reduction
PND5.CoDES.Epith <- readRDS("PND5CoDESEpith.Rds")

DefaultAssay(PND5.CoDES.Epith) <- "peaks"
PND5.CoDES.Epith <- FindTopFeatures(PND5.CoDES.Epith, min.cutoff = 'q0')
PND5.CoDES.Epith <- RunTFIDF(PND5.CoDES.Epith)
PND5.CoDES.Epith <- RunSVD(PND5.CoDES.Epith)

DepthCor(PND5.CoDES.Epith, n = 30)
ElbowPlot(PND5.CoDES.Epith, ndims = 50)

PND5.CoDES.Epith <- RunUMAP(PND5.CoDES.Epith, reduction = "lsi", dims = 5:30)
PND5.CoDES.Epith <- FindNeighbors(PND5.CoDES.Epith, reduction = 'lsi', dims = 5:30)
#kparam may need to adjust - just left to default
PND5.CoDES.Epith <- FindClusters(PND5.CoDES.Epith)

Idents(PND5.CoDES.MeStro) <- 'peaks_snn_res.0.8'
DimPlot(PND5.CoDES.Epith, label = TRUE, repel=TRUE, reduction = 'umap')
DimPlot(PND5.CoDES.Epith, label = TRUE, repel=TRUE, reduction = 'umap', split.by = 'Treatment')

Idents(PND5.CoDES.Epith) <- 'Treatment'
DimPlot(PND5.CoDES.Epith, label = TRUE, repel=TRUE, reduction = 'umap')

saveRDS(PND5.CoDES.Epith, file = "PND5CoEpiMeStro_ATACUMAP.Rds", compress = TRUE)

##find overlapping peaks----
####per treatment----
Idents(PND5.CoDES.Epith) <- 'Treatment'
DefaultAssay(PND5.CoDES.Epith) <- 'peaks'

PND5.CO.APs.E <- AccessiblePeaks(
  PND5.CoDES.Epith,
  idents = 'Control'
)
PND5.HD.APs.E <- AccessiblePeaks(
  PND5.CoDES.Epith,
  idents = 'HighDES'
)

PND5.CODES.shared.APs.E <-  PND5.CO.APs.E%in%PND5.HD.APs.E

length(PND5.CODES.shared.APs.E[PND5.CODES.shared.APs.E == TRUE])
length(PND5.CO.APs.E)
length(PND5.HD.APs.E)

library(VennDiagram)

venn.diagram(
  x = list(PND5.CO.APs.E, PND5.HD.APs.E),
  category.names = c('Control', 'DES'),
  filename = "PND5HD_CODES_Epi_Venn_ATACPeaks.png",
  output = TRUE
)

  
##Motif analysis----
###from https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html#clustering-and-further-dimensionality-reduction
#get list of motif position matrices from JASPAR
PND5CoDES.pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(species = "Mus musculus", all_versions = FALSE))
#scan the DNA seq of each peak for the presence of each motif
motif.matrix.PND5CoDES <- CreateMotifMatrix(
  features = granges(PND5.CoDES.Epith$peaks),
  pwm = PND5CoDES.pfm,
  genome = 'mm10',
  use.counts = FALSE)
dim(motif.matrix.PND5CoDES)
as.matrix(motif.matrix.PND5CoDES[1:10, 1:10])
#create new Motif object to store the results
motif.PND5CoDES <- CreateMotifObject(
  data = motif.matrix.PND5CoDES,
  pwm = PND5CoDES.pfm)
#Add the Motif object to the assay
PND5.CoDES.Epith <- SetAssayData(
  object = PND5.CoDES.Epith,
  assay = 'peaks',
  layer = 'motifs',
  new.data = motif.PND5CoDES)
PND5.CoDES.Epith$peaks@motifs

#calculate stats for each peak:
  #GC content, length, dinulceotide frequencies, etc.
  #these are used in test of motif enrichment
PND5.CoDES.Epith$peaks <- RegionStats(object = PND5.CoDES.Epith$peaks, genome = BSgenome.Mmusculus.UCSC.mm10)

#find differential accessible peaks in each cluster
PND5.CoDES.Epith <- PrepSCTFindMarkers(PND5.CoDES.Epith)
DA.PND5CoDES <- FindAllMarkers(
  object = PND5.CoDES.Epith,
  only.pos = TRUE,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments')
write.csv(DA.PND5CoDES, 'AllDApeaksPND5CoDES.csv', row.names = TRUE)

head(DA.PND5CoDES)

#get top DAs with lowest p-values
top.DA.peak.PND5CoDES <- rownames(DA.PND5CoDES[DA.PND5CoDES$p_val < 0.01, ])
head(top.DA.peak.PND5CoDES)
write.csv(top.DA.peak.PND5CoDES, 'topDApeaksPND5CoDES.csv', row.names = TRUE)

#find motifs enriched in top peaks
enriched.motifs.PND5CoDES <- FindMotifs(object = PND5.CoDES.Epith, features = top.DA.peak.PND5CoDES)
head(enriched.motifs.PND5CoDES)
write.csv(enriched.motifs.PND5CoDES, 'EnrichedMotifsPND5CoDES.csv', row.names = TRUE)

#visualize these motifs
MotifPlot(object = PND5.CoDES.Epith, motifs = head(rownames(enriched.motifs.PND5CoDES)))

#compute per-cell motif activity score
  #BiocManager::install("GreenleafLab/chromVAR", lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
PND5.CoDES.Epith <- RunChromVAR(
  object = PND5.CoDES.Epith,
  genome = BSgenome.Mmusculus.UCSC.mm10)
PND5.CoDES.Epith$chromvar
GetAssayData(PND5.CoDES.Epith$chromvar)[1:10, 1:3]

#check specific motifs
DefaultAssay(PND5.CoDES.Epith) <- "chromvar"
Idents(PND5.CoDES.Epith) <- 'Treatment'

FeaturePlot(
  object = PND5.CoDES.Epith,
  features = 'MA0489.2',
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  split.by = 'Treatment',
  pt.size = 1)

VlnPlot(object = PND5.CoDES.Epith, features = 'MA0489.2', group.by = 'Treatment', add.noise = FALSE)

####save RDS with motif scores----
saveRDS(PND5.CoDES.Epith, file = "PND5CoDESEpiMeStro_ATAC_chromvar.Rds", compress = TRUE)

###link peaks to genes----
DefaultAssay(PND5.CoDES.Epith) <- 'ATAC'
PND5.CoDES.Epith <- RegionStats(PND5.CoDES.Epith, genome = BSgenome.Mmusculus.UCSC.mm10)
PND5.CoDES.Epith <- LinkPeaks(
  object = PND5.CoDES.Epith,
  peak.assay = 'ATAC', 
  expression.assay = 'SCT'
)

saveRDS(PND5.CoDES.Epith, file = "PND5CoDESEpiMeStro_ATAC_linked_peak.rds")

###find differentially accessible peaks----
####per treatment----
Idents(PND5.CoDES.Epith) <- 'Treatment'
DefaultAssay(PND5.CoDES.Epith) <- 'peaks'
#lets try dropping the min% from 0.2 to 0.01
#test that used to take 1 min now takes 1hr
PND5CoDESE.DApeaks.HDtoC <- FindMarkers(
  object = PND5.CoDES.Epith,
  ident.1 = 'HighDES',
  ident.2 = 'Control',
  only.pos = FALSE,
  min.pct = 0.01,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments',
  assay = 'peaks',
  verbose = TRUE
)
PND5CoDESE.DApeaks.CtoHD<- FindMarkers(
  object = PND5.CoDES.Epith,
  ident.1 = 'Control',
  ident.2 = 'HighDES',
  only.pos = FALSE,
  min.pct = 0.01,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments',
  assay = 'peaks',
  verbose = TRUE
) #same output! redundant and unnecessary

write.csv(PND5CoDESE.DApeaks.CtoHD, "PND5CoDESEpiMeStro_ATAC_DApeaks_CtoHD.csv", row.names = TRUE)

write.csv(PND5CoDESE.DApeaks.HDtoC, "PND5CoDESEpiMeStro_ATAC_DApeaks_HDtoC.csv", row.names = TRUE)
  PND5CoDESE.DApeaks.HDtoC <- read.csv('PND5CoDESEpiMeStro_ATAC_DApeaks_HDtoC.csv')

#find peaks with sig p-values
PND5CoDESE.top.DApeak.HDtoC <- PND5CoDESE.DApeaks.HDtoC[PND5CoDESE.DApeaks.HDtoC$p_val < 0.05, ]
PND5CoDESE.top.DApeak.HDtoC.gain <- dplyr::filter(PND5CoDESE.top.DApeak.HDtoC, avg_log2FC <= 0.5)
PND5CoDESE.top.DApeak.HDtoC.loss <- dplyr::filter(PND5CoDESE.top.DApeak.HDtoC, avg_log2FC >= -0.5)

head(PND5CoDESE.top.DApeak.HDtoC.gain)
head(PND5CoDESE.top.DApeak.HDtoC.loss)

#find motifs enriched in top sig DAs
##in peaks gained w DES
PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks <- FindMotifs(
  object <- PND5.CoDES.Epith,
  features = PND5CoDESE.top.DApeak.HDtoC.gain$X,
  assay = 'peaks'
)
head(PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks)

write.csv(PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks, "PND5CoDESEpiMeStro_ATAC_enrichedmotifs_HDtoC_gainedpeaks.csv", row.names = TRUE)

DefaultAssay(PND5.CoDES.Epith) <- 'peaks'
MotifPlot(
  object = PND5.CoDES.Epith,
  motifs = head(PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks$motif, 10)
)

##in peaks lost w DES
PND5CoDESE.enrichedmotifs.HDtoC.lostpeaks <- FindMotifs(
  object <- PND5.CoDES.Epith,
  features = PND5CoDESE.top.DApeak.HDtoC.loss$X,
  assay = 'peaks'
)
head(PND5CoDESE.enrichedmotifs.HDtoC.lostpeaks)

write.csv(PND5CoDESE.enrichedmotifs.HDtoC.lostpeaks, "PND5CoDESEpiMeStro_ATAC_enrichedmotifs_HDtoC_lostpeaks.csv", row.names = TRUE)

MotifPlot(
  object = PND5.CoDES.Epith,
  motifs = head(rownames(PND5CoDESE.enrichedmotifs.HDtoC))
)

DefaultAssay(PND5.CoDES.Epith) <- 'peaks'
MotifPlot(
  object = PND5.CoDES.Epith,
  motifs = head(PND5CoDESE.enrichedmotifs.HDtoC.lostpeaks$motif, 10)
)

###annotate regions----
#find genes near or at genomic regions
head(PND5CoDESE.DApeaks.HDtoC)

PND5CoDESE.DApeaks_open_HDtoC <- PND5CoDESE.DApeaks.HDtoC %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC >= 1.5 | avg_log2FC <= -1.5) %>%
  arrange(avg_log2FC)
PND5CoDESEclosestgene_HDtoC <- ClosestFeature(PND5.CoDES.Epith$peaks, regions = PND5CoDESE.DApeaks_open_HDtoC$X)

write.csv(PND5CoDESEclosestgene_HDtoC, "PND5CoDESEpithDAclosestgene_HDtoC.csv", row.names = TRUE)
  PND5CoDESEclosestgene_HDtoC <- read.csv("PND5CoDESEpithDAclosestgene_HDtoC.csv")
  
PND5CoDESEclosestgene_HDtoC[PND5CoDESEclosestgene_HDtoC$gene_name == 'Padi1', ]

#####combine tables----
PND5CoDESEclosestgene_peakinfo_HDtoC <- merge(PND5CoDESE.DApeaks.HDtoC, PND5CoDESEclosestgene_HDtoC, by.x = 'X', by.y = 'query_region')
head(PND5CoDESEclosestgene_peakinfo_HDtoC)
  write.csv(PND5CoDESEclosestgene_peakinfo_HDtoC, "PND5CoDESEpithDAclosestgene_withpeakinfo_HDtoC.csv", row.names = TRUE)

#####volcano plot - differential activity ----
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#
PND5CoDESEclosestgene_peakinfo_HDtoC.sig <- dplyr::filter(PND5CoDESEclosestgene_peakinfo_HDtoC, p_val_adj <= 0.01)

PND5CoDESEclosestgene_peakinfo_HDtoC.up <- dplyr::filter(PND5CoDESEclosestgene_peakinfo_HDtoC.sig, avg_log2FC >= 1.5)
head(PND5CoDESEclosestgene_peakinfo_HDtoC.up)
length(PND5CoDESEclosestgene_peakinfo_HDtoC.up$X)

PND5CoDESEclosestgene_peakinfo_HDtoC.dn <- dplyr::filter(PND5CoDESEclosestgene_peakinfo_HDtoC.sig, avg_log2FC <= -0.5)
head(PND5CoDESEclosestgene_peakinfo_HDtoC.dn)
length(PND5CoDESEclosestgene_peakinfo_HDtoC.dn$X)

EnhancedVolcano(PND5CoDESEclosestgene_peakinfo_HDtoC, lab = PND5CoDESEclosestgene_peakinfo_HDtoC$JASPARnames.M,
                x = 'avg_log2FC', y = 'p_val_adj', title = 'Control vs DES',
                pCutoff = 0.05, pointSize = 1.0, labSize = 5.0,
                ylim = c(0, -log10(10e-300)), drawConnectors = FALSE,
                selectLab = c(PND5CoDESEclosestgene_peakinfo_HDtoC.dn$JASPARnames.M, PND5CoDESEclosestgene_peakinfo_HDtoC.up$JASPARnames.M))


  
#####differential activity between treatments----
DefaultAssay(PND5.CoDES.Epith) <- 'chromvar'
Idents(PND5.CoDES.Epith) <- 'Treatment'

PND5CoDESEATAC.differential.activity <- FindMarkers(
  object = PND5.CoDES.Epith,
  ident.1 = 'HighDES',
  ident.2 = 'Control',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = 'avg_diff',
  verbose = TRUE
)
  
MotifPlot(
  object = PND5.CoDES.Epith,
  motifs = head(rownames(PND5CoDESEATAC.differential.activity)),
  assay = 'peaks'
)

write.csv(PND5CoDESEATAC.differential.activity, "PND5CoDESEpithDiffAct_HDtoC.csv", row.names = TRUE)

diffAct <- dplyr::filter(PND5CoDESEATAC.differential.activity, p_val_adj <= 0.05)
diffAct <- dplyr::filter(diffAct, avg_diff >= 1.5 | avg_diff <= -1.5)
diffAct <- diffAct[order(diffAct$avg_diff),]
head(diffAct,12)
  
#convert JASPAR IDs to names
PND5CoDESEJASPARidsrownames <- PND5CoDESEATAC.differential.activity$X
JASPARnames.M <- ConvertMotifID(PND5.CoDES.Epith$peaks@motifs, id = PND5CoDESEJASPARidsrownames)
  
M.IDName <- data.frame(PND5CoDESEJASPARidsrownames, JASPARnames.M)
  
PND5CoDESEATAC.differential.activity.name <- merge(PND5CoDESEATAC.differential.activity, M.IDName, by.x = 'X', by.y = 'PND5CoDESEJASPARidsrownames')
head(PND5CoDESEATAC.differential.activity.name,1)
  
write.csv(PND5CoDESEATAC.differential.activity.name, "PND5CoDESEpithDiffAct_HDtoC.csv", row.names = TRUE)

#####volcano plot - differential activity ----
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
PND5CoDESEATAC.differential.activity.name.sig <- dplyr::filter(PND5CoDESEATAC.differential.activity.name, p_val_adj <= 0.05)

PND5CoDESEATAC.differential.activity.name.up <- dplyr::filter(PND5CoDESEATAC.differential.activity.name.sig, avg_diff >= 1)
head(PND5CoDESEATAC.differential.activity.name.up)

PND5CoDESEATAC.differential.activity.name.dn <- dplyr::filter(PND5CoDESEATAC.differential.activity.name.sig, avg_diff <= -1)
head(PND5CoDESEATAC.differential.activity.name.dn)
  
EnhancedVolcano(PND5CoDESEATAC.differential.activity.name, lab = PND5CoDESEATAC.differential.activity.name$JASPARnames.M,
                x = 'avg_diff', y = 'p_val_adj', title = 'Control vs DES',
                pCutoff = 0.05, pointSize = 1.0, labSize = 5.0,
                ylim = c(0, -log10(10e-300)), drawConnectors = FALSE,
                selectLab = c(PND5CoDESEATAC.differential.activity.name.dn$JASPARnames.M, PND5CoDESEATAC.differential.activity.name.up$JASPARnames.M))

#get cell ids for Browser Track output -----
DefaultAssay(PND5.CoDES.Epith) <- 'ATAC'
Idents(PND5.CoDES.Epith) <- 'Treatment'

PND5.CoDES.Epith.C <- subset(PND5.CoDES.Epith, idents='Control')
write(colnames(PND5.CoDES.Epith.C), file="EpiConCellIDlist_trim.txt", ncolumns=1);

PND5.CoDES.Epith.HD <- subset(PND5.CoDES.Epith, idents='HighDES')
write(colnames(PND5.CoDES.Epith.HD), file="EpiHDCellIDlist_trim.txt", ncolumns=1);

#TSS enrichment----
gene.ranges <- genes(EnsDb.Mmusculus.v79)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
tss.ranges <- GRanges(
  seqnames = seqnames(gene.ranges),
  ranges = IRanges(start = start(gene.ranges), width = 2),
  strand = strand(gene.ranges)
)
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges = keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

PND5.CoDES.Epith <-  TSSEnrichment(
  PND5.CoDES.Epith,
  tss.positions = tss.ranges,
  fast = FALSE #if don't set fast = FALSE, will not make position enrichment matrix
)

Idents(PND5.CoDES.Epith) <- 'Treatment'
TSSPlot(PND5.CoDES.Epith)

Idents(PND5.CoDES.Epith) <- 'peaks_snn_res.0.8'
TSSPlot(PND5.CoDES.Epith)

saveRDS(PND5.CoDES.Epith, file = "PND5CoDESEpi_ATAC_TSS.rds")

PND5.CoDES.Epith$peaks
  
#find DAPs within 200kb of TSS of DEGs----
#use closest gene file and match to DEGs

PND5CoDESEclosestgene_peakinfo_HDtoC <- dplyr::filter(PND5CoDESEclosestgene_peakinfo_HDtoC, distance <= 200000) %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC <= 1.2 | avg_log2FC >= -1.2)

PND5.CODES.Epith.markers.HDtoC <- dplyr::filter(PND5.CODES.Epith.markers.HDtoC, p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC <= 1.2 | avg_log2FC >= -1.2)

PND5.CODES.Epith.DEG.DAP.match <- inner_join(x = PND5.CODES.Epith.markers.HDtoC, y = PND5CoDESEclosestgene_peakinfo_HDtoC, by = c('X' = 'gene_name'), suffix = c('_DEG', '_DAP'))
head(PND5.CODES.Epith.DEG.DAP.match,1)

write.csv(PND5.CODES.Epith.DEG.DAP.match, "PND5CoEpith_DEGDAPmatch.csv", row.names = TRUE)

###then make coverage plots-----
Idents(PND5.CoDES.Epith) <- 'Treatment'
DefaultAssay(PND5.CoDES.Epith) <- 'ATAC'

CoveragePlot(
  object = PND5.CoDES.Epith,
  region = 'chr4-140774089-140949211',
  extend.upstream = 200000,
  extend.downstream = 200000,
  features = c('Padi1', 'Padi2', 'Padi3', 'Padi4'),
  expression.assay = 'SCT', 
  ymax = 30
)
```

Now to subset Mesenchyme/Stroma. Repeat same basic visualization and clustering
and basic analysis as in epithelium. Probe for genes of interest in custom
feature and violin plots.

```{r Subset MeStro from reintegrated EpiMeStro}
PND5.CoDES.EpiMeStro <- readRDS("PND5CoDESEpiMeStro_annotatedfromUnintegrated.Rds")
Idents(PND5.CoDES.EpiMeStro) <- 'EpiMeStro_Annotated'
PND5.CoDES.MeStro <- subset(PND5.CoDES.EpiMeStro, idents = 'Mesenchyme/Stroma')
table(PND5.CoDES.MeStro$EpiMeStro_Annotated)
table(PND5.CoDES.MeStro$Treatment)

## standard visualization and clustering----
DefaultAssay(PND5.CoDES.MeStro) <- "RNA"
PND5.CoDES.MeStro <- SCTransform(PND5.CoDES.MeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 

DepthCor(PND5.CoDES.MeStro, reduction = 'pca', n = 50)
ElbowPlot(PND5.CoDES.MeStro, ndims = 50)

PND5.CoDES.MeStro <- FindNeighbors(PND5.CoDES.MeStro, dims = c(1:20)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:20)

DimPlot(PND5.CoDES.MeStro, label = FALSE, repel=TRUE)

PND5.CoDES.MeStro[["MeStro_subcluster"]] <- Idents(PND5.CoDES.MeStro)

##check if any cluster has low QC metrics----
Idents(PND5.CoDES.MeStro) <- 'MeStro_subcluster'
VlnPlot(
  object = PND5.CoDES.MeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)
#1 and 4 worse, but from DES

saveRDS(PND5.CoDES.MeStro, file = "PND5CoDESMeStro.Rds", compress = TRUE)

####DEGs----
DefaultAssay(PND5.CoDES.MeStro) <- "SCT"
Idents(PND5.CoDES.MeStro) <- 'Treatment'

PND5.CoDES.MeStro <- PrepSCTFindMarkers(PND5.CoDES.MeStro)

PND5.CODES.MeStro.markers.all <- FindMarkers(PND5.CoDES.MeStro, test.use = 'MAST',
                                            only.pos = FALSE, ident.1 = 'HighDES',
                                            ident.2 = 'Control',
                           assay = 'SCT', latent.vars='percent.mt')

write.csv(PND5.CODES.MeStro.markers.all, "PND5CoDESMeStro.DEG.DEStoCon.csv")

####DEG volcano----
head(PND5.CODES.MeStro.markers.all,1)

ggplot(PND5.CODES.MeStro.markers.all, aes(avg_log2FC, -log10(p_val))) + 
  geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = 
  ifelse(p_val_adj < 0.01, X, "")), colour = "red", size = 3, max.overlaps=5)

#####----Select Marker Genes----
DefaultAssay(PND5.CoDES.MeStro) <- "SCT"
Idents(PND5.CoDES.MeStro) <- 'MeStro_subcluster'

FeaturePlot(object = PND5.CoDES.MeStro,
            features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
DotPlot(object = PND5.CoDES.MeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11')) + coord_flip()
```

Clusters are separating entirely by treatment. I will need to subset and analyze DES and control separately.

```{r Subset MeStro by Treatment}
PND5.CoDES.MeStro <- readRDS("PND5CoDESMeStro.Rds")
Idents(PND5.CoDES.MeStro) <- 'Treatment'
PND5.Co.MeStro <- subset(PND5.CoDES.MeStro, idents = 'Control')
PND5.DES.MeStro <- subset(PND5.CoDES.MeStro, idents = 'HighDES')
table(PND5.Co.MeStro$EpiMeStro_Annotated)
table(PND5.DES.MeStro$EpiMeStro_Annotated)
table(PND5.Co.MeStro$Treatment)
table(PND5.DES.MeStro$Treatment)

## standard visualization and clustering----
DefaultAssay(PND5.Co.MeStro) <- "RNA"
PND5.Co.MeStro <- SCTransform(PND5.Co.MeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 
DefaultAssay(PND5.DES.MeStro) <- "RNA"
PND5.DES.MeStro <- SCTransform(PND5.DES.MeStro) %>% 
  RunPCA(features = VariableFeatures(object = .)) 

DepthCor(PND5.Co.MeStro, reduction = 'pca', n = 50)
ElbowPlot(PND5.Co.MeStro, ndims = 50)
DepthCor(PND5.DES.MeStro, reduction = 'pca', n = 50)
ElbowPlot(PND5.DES.MeStro, ndims = 50)

PND5.Co.MeStro <- FindNeighbors(PND5.Co.MeStro, dims = c(1:30)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:30)
PND5.DES.MeStro <- FindNeighbors(PND5.DES.MeStro, dims = c(1:30)) %>% 
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims = 1:30)

DimPlot(PND5.Co.MeStro, label = FALSE, repel=TRUE)
DimPlot(PND5.DES.MeStro, label = FALSE, repel=TRUE)


PND5.Co.MeStro[["MeStro_subcluster"]] <- Idents(PND5.Co.MeStro)
PND5.DES.MeStro[["MeStro_subcluster"]] <- Idents(PND5.DES.MeStro)

##check if any cluster has low QC metrics----
Idents(PND5.Co.MeStro) <- 'MeStro_subcluster'
Idents(PND5.DES.MeStro) <- 'MeStro_subcluster'
VlnPlot(
  object = PND5.Co.MeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)
VlnPlot(
  object = PND5.DES.MeStro,
  features = c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment',
               'nucleosome_signal', 'percent.mt', 'percent.rp',
               'blacklistratio', 'pct_reads_in_peaks'),
  ncol = 3, pt.size = 0)

saveRDS(PND5.Co.MeStro, file = "PND5CoMeStro.Rds", compress = TRUE)
saveRDS(PND5.DES.MeStro, file = "PND5DESMeStro.Rds", compress = TRUE)

##annotate----
###find markers per cluster----
####crashes R if using RNA, works if using SCT!-----
options(future.globals.maxSize= 1572864000) #allocate 1.5 Gb
#if crashing, reset environment and redo (quit q() and do not save env)
head(PND5.Co.MeStro,1)
Idents(PND5.Co.MeStro) <- 'MeStro_subcluster'
Idents(PND5.DES.MeStro) <- 'MeStro_subcluster'

###Control----
PND5.Co.MeStro.markers.all <- FindAllMarkers(PND5.Co.MeStro, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.Co.MeStro.markers.all, "PND5CoMeStroMarkers_byCluster.csv", row.names = TRUE)

PND5.Co.MeStro.markers.all <- dplyr::filter(PND5.Co.MeStro.markers.all,p_val_adj <= 0.05) %>%
    dplyr::filter(, avg_log2FC > 0.5)

stroma_markers <- c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6', 'Dpt', 'Vcan', 'Col6a3')
mesenchyme_markers <-  c('Hoxa10', 'Vim', 'Cnn1', 'Tcf21', 'Pkib', 'Pcdh10')
epithelium_markers <- c('Epcam', 'Cdh1', 'Krt18')
myocyte_markers <- c('Mef2a', 'Pdlim3', 'Acta2', 'Chrm3', 'Myh11')
markers <- c(stroma_markers, mesenchyme_markers,  epithelium_markers, 
             myocyte_markers)
tissue <- c(replicate(length(stroma_markers), 'stroma'), 
            replicate(length(mesenchyme_markers), 'mesenchyme'), 
            replicate(length(epithelium_markers), 'epithelium'), 
            replicate(length(myocyte_markers), 'myocyte'))
markerdict <- data.frame(markers, tissue)

PND5.Co.MeStro.markers.match <- merge(x = PND5.Co.MeStro.markers.all, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.Co.MeStro.markers.match, "PND5CoMeStro_Markers_byCluster_matchmarkers.csv", row.names = FALSE)

###DES----
PND5.DES.MeStro.markers.all <- FindAllMarkers(PND5.DES.MeStro, test.use = 'MAST', only.pos = TRUE,
                           assay = 'SCT', latent.vars='percent.mt')
write.csv(PND5.DES.MeStro.markers.all, "PND5DESMeStroMarkers_byCluster.csv", row.names = TRUE)

PND5.DES.MeStro.markers.all <- dplyr::filter(PND5.DES.MeStro.markers.all,p_val_adj <= 0.05) %>%
    dplyr::filter(, avg_log2FC > 0.5)

stroma_markers <- c('Pdgfra', 'Dcn', 'Col15a1', 'Foxl2', 'Smad6', 'Dpt', 'Vcan', 'Col6a3')
mesenchyme_markers <-  c('Hoxa10', 'Vim', 'Cnn1', 'Tcf21', 'Pkib', 'Pcdh10')
epithelium_markers <- c('Epcam', 'Cdh1', 'Krt18')
myocyte_markers <- c('Mef2a', 'Pdlim3', 'Acta2', 'Chrm3', 'Myh11')
markers <- c(stroma_markers, mesenchyme_markers,  epithelium_markers, 
             myocyte_markers)
tissue <- c(replicate(length(stroma_markers), 'stroma'), 
            replicate(length(mesenchyme_markers), 'mesenchyme'), 
            replicate(length(epithelium_markers), 'epithelium'), 
            replicate(length(myocyte_markers), 'myocyte'))
markerdict <- data.frame(markers, tissue)

PND5.DES.MeStro.markers.match <- merge(x = PND5.DES.MeStro.markers.all, y = markerdict, by.x = 'gene', by.y = 'markers')
write.csv(PND5.DES.MeStro.markers.match, "PND5DESMeStro_Markers_byCluster_matchmarkers.csv", row.names = FALSE)


##print plots of cell markers----
###Control----
DefaultAssay(PND5.Co.MeStro) <- "SCT"
Idents(PND5.Co.MeStro) <- 'MeStro_subcluster'

#example
FeaturePlot(object = PND5.Co.MeStro,
            features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
DotPlot(object = PND5.Co.MeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11')) + coord_flip()

###DES----
DefaultAssay(PND5.DES.MeStro) <- "SCT"
Idents(PND5.DES.MeStro) <- 'MeStro_subcluster'

###All Tissue Markers----
FeaturePlot(object = PND5.DES.MeStro,
            features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11'),
            order = TRUE, label = TRUE, repel = TRUE,
            reduction = 'umap', combine = TRUE)
DotPlot(object = PND5.DES.MeStro, features = c('Hoxa10', 'Vim', 'Epcam', 'Krt18', 'Pdgfra', 'Vcan', 'Dcn', 'Myh11')) + coord_flip()

#cell cycle----
###load in cell cycle markers included with Seurat----
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

###assign cell cycle scores----
PND5.Co.MeStro <-  CellCycleScoring(PND5.Co.MeStro, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
PND5.DES.MeStro <-  CellCycleScoring(PND5.DES.MeStro, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

###visualize----
# Visualize the distribution of cell cycle markers across
RidgePlot(PND5.Co.MeStro, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2, group.by = "Phase")
RidgePlot(PND5.DES.MeStro, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2, group.by = "Phase")

#show on current UMAP
DimPlot(PND5.Co.MeStro, label = FALSE, repel=TRUE, group.by = "Phase")
DimPlot(PND5.DES.MeStro, label = FALSE, repel=TRUE, group.by = "Phase")

#merge co and DES for shared visualization---
DefaultAssay(PND5.Co.MeStro) <- "SCT"
Idents(PND5.Co.MeStro) <- 'MeStro_subcluster'

DefaultAssay(PND5.DES.MeStro) <- "SCT"
Idents(PND5.DES.MeStro) <- 'MeStro_subcluster'

PND5.merge.MeStro <- merge(PND5.Co.MeStro, PND5.DES.MeStro, merge.data = TRUE, merge.dr = TRUE)
head(PND5.merge.MeStro,1)
PND5.merge.MeStro$SCT

saveRDS(PND5.merge.MeStro, 'PND5mergedMeStro.Rds')

###merge for visualization between treatments----
Idents(PND5.merge.MeStro) <- 'Treatment'

#example
FeaturePlot(object = PND5.merge.MeStro,
            features = 'Mki67',
            order = TRUE, label = FALSE, repel = TRUE, split.by = 'Treatment',
            reduction = 'umap', keep.scale = 'all', combine = TRUE) + theme(legend.position = 'right') & NoAxes() & ggtitle(NULL)

```

Now to look at ATAC

```{r MeStro ATAC - Integrated for peak comparison between Treatments}
#UMAP based on ATAC----
#from https://nbis-workshop-MeStrogenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html#clustering-and-further-dimensionality-reduction
PND5.CoDES.MeStro <- readRDS("PND5CoDESMeStro.Rds")

DefaultAssay(PND5.CoDES.MeStro) <- "peaks"
PND5.CoDES.MeStro <- FindTopFeatures(PND5.CoDES.MeStro, min.cutoff = 'q0')
PND5.CoDES.MeStro <- RunTFIDF(PND5.CoDES.MeStro)
PND5.CoDES.MeStro <- RunSVD(PND5.CoDES.MeStro)

DepthCor(PND5.CoDES.MeStro, n = 30)
ElbowPlot(PND5.CoDES.MeStro, ndims = 50)

PND5.CoDES.MeStro <- RunUMAP(PND5.CoDES.MeStro, reduction = "lsi", dims = 3:30)
PND5.CoDES.MeStro <- FindNeighbors(PND5.CoDES.MeStro, reduction = 'lsi', dims = 3:30)
#kparam may need to adjust - just left to default
PND5.CoDES.MeStro <- FindClusters(PND5.CoDES.MeStro)

Idents(PND5.CoDES.MeStro) <- 'peaks_snn_res.0.8'
DimPlot(PND5.CoDES.MeStro, label = TRUE, repel=TRUE, reduction = 'umap')

Idents(PND5.CoDES.MeStro) <- 'Treatment'
DimPlot(PND5.CoDES.MeStro, label = TRUE, repel=TRUE, reduction = 'umap')

saveRDS(PND5.CoDES.MeStro, file = "PND5CoMeStroMeStro_ATACUMAP.Rds", compress = TRUE)

##find overlapping peaks----
####per treatment----
Idents(PND5.CoDES.MeStro) <- 'Treatment'
DefaultAssay(PND5.CoDES.MeStro) <- 'peaks'

PND5.CO.APs.MS <- AccessiblePeaks(
  PND5.CoDES.MeStro,
  idents = 'Control'
)
PND5.HD.APs.MS <- AccessiblePeaks(
  PND5.CoDES.MeStro,
  idents = 'HighDES'
)

PND5.CODES.shared.APs.MS <-  PND5.CO.APs.MS%in%PND5.HD.APs.MS

length(PND5.CODES.shared.APs.MS[PND5.CODES.shared.APs.MS == TRUE])
length(PND5.CO.APs.MS)
length(PND5.HD.APs.MS)

library(VennDiagram)

venn.diagram(
  x = list(PND5.CO.APs.MS, PND5.HD.APs.MS),
  category.names = c('Control', 'DES'),
  filename = "PND5HD_CODES_MeStro_Venn_ATACPeaks.png",
  output = TRUE
)

##Motif analysis----
###from https://nbis-workshop-MeStrogenomics.readthedocs.io/en/latest/content/tutorials/scAtacSeq/lab-sc_atac_seq.html#clustering-and-further-dimensionality-reduction
#get list of motif position matrices from JASPAR
PND5CoDES.pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(species = "Mus musculus", all_versions = FALSE))
#scan the DNA seq of each peak for the presence of each motif
motif.matrix.PND5CoDES <- CreateMotifMatrix(
  features = granges(PND5.CoDES.MeStro$peaks),
  pwm = PND5CoDES.pfm,
  genome = 'mm10',
  use.counts = FALSE)
dim(motif.matrix.PND5CoDES)
as.matrix(motif.matrix.PND5CoDES[1:10, 1:10])
#create new Motif object to store the results
motif.PND5CoDES <- CreateMotifObject(
  data = motif.matrix.PND5CoDES,
  pwm = PND5CoDES.pfm)
#Add the Motif object to the assay
PND5.CoDES.MeStro <- SetAssayData(
  object = PND5.CoDES.MeStro,
  assay = 'peaks',
  layer = 'motifs',
  new.data = motif.PND5CoDES)
PND5.CoDES.MeStro$peaks@motifs

#calculate stats for each peak:
  #GC content, length, dinulceotide frequencies, etc.
  #these are used in test of motif enrichment
PND5.CoDES.MeStro$peaks <- RegionStats(object = PND5.CoDES.MeStro$peaks, genome = BSgenome.Mmusculus.UCSC.mm10)

#find differential accessible peaks in each cluster
PND5.CoDES.MeStro <- PrepSCTFindMarkers(PND5.CoDES.MeStro)
DA.PND5CoDES <- FindAllMarkers(
  object = PND5.CoDES.MeStro,
  only.pos = TRUE,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments')
write.csv(DA.PND5CoDES, 'AllDApeaksPND5CoDES.csv', row.names = TRUE)

head(DA.PND5CoDES)

#get top DAs with lowest p-values
top.DA.peak.PND5CoDES <- rownames(DA.PND5CoDES[DA.PND5CoDES$p_val < 0.01, ])
head(top.DA.peak.PND5CoDES)
write.csv(top.DA.peak.PND5CoDES, 'topDApeaksPND5CoDES.csv', row.names = TRUE)

#find motifs enriched in top peaks
enriched.motifs.PND5CoDES <- FindMotifs(object = PND5.CoDES.MeStro, features = top.DA.peak.PND5CoDES)
head(enriched.motifs.PND5CoDES)
write.csv(enriched.motifs.PND5CoDES, 'EnrichedMotifsPND5CoDES.csv', row.names = TRUE)
#visualize these motifs
MotifPlot(object = PND5.CoDES.MeStro, motifs = head(rownames(enriched.motifs.PND5CoDES)))

#compute per-cell motif activity score
PND5.CoDES.MeStro <- RunChromVAR(
  object = PND5.CoDES.MeStro,
  genome = BSgenome.Mmusculus.UCSC.mm10)
PND5.CoDES.MeStro$chromvar
GetAssayData(PND5.CoDES.MeStro$chromvar)[1:10, 1:3]

#check specific motifs - sanity check
DefaultAssay(PND5.CoDES.MeStro) <- "chromvar"
FeaturePlot(
  object = PND5.CoDES.MeStro,
  features = 'MA1684.1',
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 1)

####save RDS with motif scores----
saveRDS(PND5.CoDES.MeStro, file = "PND5CoDESMeStroMeStro_ATAC_chromvar.Rds", compress = TRUE)

###link peaks to genes----
DefaultAssay(PND5.CoDES.MeStro) <- 'ATAC'
PND5.CoDES.MeStro <- RegionStats(PND5.CoDES.MeStro, genome = BSgenome.Mmusculus.UCSC.mm10)
PND5.CoDES.MeStro <- LinkPeaks(
  object = PND5.CoDES.MeStro,
  peak.assay = 'ATAC', 
  expression.assay = 'SCT'
)

saveRDS(PND5.CoDES.MeStro, file = "PND5CoDESMeStroMeStro_ATAC_linked_peak.rds")

###find differentially accessible peaks----
####per treatment----
Idents(PND5.CoDES.MeStro) <- 'Treatment'
DefaultAssay(PND5.CoDES.MeStro) <- 'peaks'
#lets try dropping the min% from 0.2 to 0.01
#test that used to take 1 min now takes 1hr
PND5CoDESE.DApeaks.HDtoC <- FindMarkers(
  object = PND5.CoDES.MeStro,
  ident.1 = 'HighDES',
  ident.2 = 'Control',
  only.pos = FALSE,
  min.pct = 0.01,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments',
  assay = 'peaks',
  verbose = TRUE
)

write.csv(PND5CoDESE.DApeaks.HDtoC, "PND5CoDESMeStro_ATAC_DApeaks_HDtoC.csv", row.names = TRUE)

#find peaks with sig p-values
PND5CoDESE.top.DApeak.HDtoC <- PND5CoDESE.DApeaks.HDtoC[PND5CoDESE.DApeaks.HDtoC$p_val < 0.05, ]
PND5CoDESE.top.DApeak.HDtoC.gain <- dplyr::filter(PND5CoDESE.top.DApeak.HDtoC, avg_log2FC <= 0.5)
PND5CoDESE.top.DApeak.HDtoC.loss <- dplyr::filter(PND5CoDESE.top.DApeak.HDtoC, avg_log2FC >= -0.5)

head(PND5CoDESE.top.DApeak.HDtoC.gain)
head(PND5CoDESE.top.DApeak.HDtoC.loss)

#find motifs enriched in top sig DAs
PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks <- FindMotifs(
  object <- PND5.CoDES.MeStro,
  features = PND5CoDESE.top.DApeak.HDtoC.gain$X,
  assay = 'peaks'
)
head(PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks)

PND5CoDESE.enrichedmotifs.HDtoC.lostpeaks <- FindMotifs(
  object <- PND5.CoDES.MeStro,
  features = PND5CoDESE.top.DApeak.HDtoC.loss$X,
  assay = 'peaks'
)
head(PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks)

write.csv(PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks, "PND5CoDESMeStro_ATAC_enrichedmotifs_HDtoC_gainedpeaks.csv", row.names = TRUE)
write.csv(PND5CoDESE.enrichedmotifs.HDtoC.lostpeaks, "PND5CoDESMeStro_ATAC_enrichedmotifs_HDtoC_lostpeaks.csv", row.names = TRUE)

MotifPlot(
  object = PND5.CoDES.MeStro,
  motifs = head(rownames(PND5CoDESE.enrichedmotifs.HDtoC))
)

DefaultAssay(PND5.CoDES.MeStro) <- 'peaks'
MotifPlot(
  object = PND5.CoDES.MeStro,
  motifs = head(PND5CoDESE.enrichedmotifs.HDtoC.gainedpeaks$motif.name, 10)
)
MotifPlot(
  object = PND5.CoDES.MeStro,
  motifs = head(PND5CoDESE.enrichedmotifs.HDtoC.lostpeaks$motif.name, 10)
)

###annotate regions----
#find genes near or at genomic regions
PND5CoDESE.DApeaks_open_HDtoC <- PND5CoDESE.DApeaks.HDtoC %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC >= 1.5 | avg_log2FC <= -1.5) %>%
  arrange(avg_log2FC)
PND5CoDESEclosestgene_HDtoC <- ClosestFeature(PND5.CoDES.MeStro$peaks, regions = PND5CoDESE.DApeaks_open_HDtoC$X)

head(PND5CoDESEclosestgene_HDtoC)

write.csv(PND5CoDESEclosestgene_HDtoC, "PND5CoDESMeStroDAclosestgene_HDtoC.csv", row.names = TRUE)

#####combine tables----
PND5CoDESE.DApeaks.HDtoC <- read.csv('PND5CoDESMeStro_ATAC_DApeaks_HDtoC.csv')
PND5CoDESEclosestgene_HDtoC <- read.csv("PND5CoDESMeStroDAclosestgene_HDtoC.csv")
  
PND5CoDESEclosestgene_peakinfo_HDtoC <- merge(PND5CoDESE.DApeaks.HDtoC, PND5CoDESEclosestgene_HDtoC, by.x = 'X', by.y = 'query_region')
head(PND5CoDESEclosestgene_peakinfo_HDtoC)
write.csv(PND5CoDESEclosestgene_peakinfo_HDtoC, "PND5CoDESMeStroDAclosestgene_withpeakinfo_HDtoC.csv", row.names = TRUE)

#####differenatial activity between treatments----
PND5.CoDES.MeStro <- readRDS("PND5CoDESMeStroMeStro_ATAC_linked_peak.rds")
  
DefaultAssay(PND5.CoDES.MeStro) <- 'chromvar'
Idents(PND5.CoDES.MeStro) <- 'Treatment'

PND5CoDESEATAC.differential.activity <- FindMarkers(
  object = PND5.CoDES.MeStro,
  ident.1 = 'HighDES',
  ident.2 = 'Control',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = 'avg_diff',
  verbose = TRUE
)
  
MotifPlot(
  object = PND5.CoDES.MeStro,
  motifs = head(rownames(PND5CoDESEATAC.differential.activity)),
  assay = 'peaks'
)

write.csv(PND5CoDESEATAC.differential.activity, "PND5CoDESMeStroDiffAct_HDtoC.csv", row.names = TRUE)

head(PND5CoDESEATAC.differential.activity)
   
diffAct <- dplyr::filter(PND5CoDESEATAC.differential.activity, p_val_adj <= 0.05)
diffAct <- dplyr::filter(diffAct, avg_diff >= 1.5 | avg_diff <= -1.5)
diffAct <- diffAct[order(diffAct$avg_diff),]
head(diffAct,12)
  
#convert JASPAR IDs to names
PND5CoDESEATAC.differential.activity <- read.csv("PND5CoDESMeStroDiffAct_HDtoC.csv")
  
PND5CoDESEJASPARidsrownames <- PND5CoDESEATAC.differential.activity$X
JASPARnames.M <- ConvertMotifID(PND5.CoDES.MeStro$peaks@motifs, id = PND5CoDESEJASPARidsrownames)
  
M.IDName <- data.frame(PND5CoDESEJASPARidsrownames, JASPARnames.M)
  
PND5CoDESEATAC.differential.activity.name <- merge(PND5CoDESEATAC.differential.activity, M.IDName, by.x = 'X', by.y = 'PND5CoDESEJASPARidsrownames')
head(PND5CoDESEATAC.differential.activity.name,1)
  
write.csv(PND5CoDESEATAC.differential.activity.name, "PND5CoDESMeStroDiffAct_HDtoC.csv", row.names = TRUE)

#####volcano plot - differential activity ----
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
PND5CoDESEATAC.differential.activity.name.sig <- dplyr::filter(PND5CoDESEATAC.differential.activity.name, p_val_adj <= 0.05)

PND5CoDESEATAC.differential.activity.name.up <- dplyr::filter(PND5CoDESEATAC.differential.activity.name.sig, avg_diff >= 1)
head(PND5CoDESEATAC.differential.activity.name.up)

PND5CoDESEATAC.differential.activity.name.dn <- dplyr::filter(PND5CoDESEATAC.differential.activity.name.sig, avg_diff <= -1)
head(PND5CoDESEATAC.differential.activity.name.dn)
  
EnhancedVolcano(PND5CoDESEATAC.differential.activity.name, lab = PND5CoDESEATAC.differential.activity.name$JASPARnames.M,
                x = 'avg_diff', y = 'p_val_adj', title = 'Control vs DES',
                pCutoff = 0.05, pointSize = 1.0, labSize = 5.0,
                ylim = c(0, -log10(10e-300)), drawConnectors = FALSE,
                selectLab = c(PND5CoDESEATAC.differential.activity.name.up$JASPARnames.M, PND5CoDESEATAC.differential.activity.name.dn$JASPARnames.M))

PND5CoDES.MeStro.topDAct <- arrange(PND5CoDESEATAC.differential.activity.name, p_val_adj)
PND5CoDES.MeStro.topDAct <- head(PND5CoDES.MeStro.topDAct, 25) %>%
  arrange(, avg_diff)

#get cell ids for Browser Track output -----
PND5.CoDES.MeStro <- readRDS("Unintegrated/PND5CoDESMeStroMeStro_ATAC_linked_peak.rds")
DefaultAssay(PND5.CoDES.MeStro) <- 'ATAC'
Idents(PND5.CoDES.MeStro) <- 'Treatment'

PND5.CoDES.MeStro.C <- subset(PND5.CoDES.MeStro, idents='Control')
write(colnames(PND5.CoDES.MeStro.C), file="MeStroConCellIDlist_trim.txt", ncolumns=1);

PND5.CoDES.MeStro.HD <- subset(PND5.CoDES.MeStro, idents='HighDES')
write(colnames(PND5.CoDES.MeStro.HD), file="MeStroHDCellIDlist_trim.txt", ncolumns=1);

#find DAPs within 200kb of TSS of DEGs----
#use closest gene file and match to DEGs

PND5CoDESEclosestgene_peakinfo_HDtoC <- read.csv("PND5CoDESMeStroDAclosestgene_withpeakinfo_HDtoC.csv")
head(PND5CoDESEclosestgene_peakinfo_HDtoC,1)

PND5CoDESEclosestgene_peakinfo_HDtoC <- dplyr::filter(PND5CoDESEclosestgene_peakinfo_HDtoC, distance <= 200000) %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC <= 1.2 | avg_log2FC >= -1.2)

PND5.CODES.MeStro.markers.HDtoC <- read.csv("PND5CoDESMeStro.DEG.DEStoCon.csv")
head(PND5.CODES.MeStro.markers.HDtoC,1)

PND5.CODES.MeStro.markers.HDtoC <- dplyr::filter(PND5.CODES.MeStro.markers.HDtoC, p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC <= 1.2 | avg_log2FC >= -1.2)

PND5.CODES.MeStro.DEG.DAP.match <- inner_join(x = PND5.CODES.MeStro.markers.HDtoC, y = PND5CoDESEclosestgene_peakinfo_HDtoC, by = c('X' = 'gene_name'), suffix = c('_DEG', '_DAP'))
head(PND5.CODES.MeStro.DEG.DAP.match,1)

write.csv(PND5.CODES.MeStro.DEG.DAP.match, "PND5CoMeStro_DEGDAPmatch.csv", row.names = TRUE)

PND5.CoDES.MeStro <- readRDS("PND5CoDESMeStroMeStro_ATAC_linked_peak.rds")
Idents(PND5.CoDES.MeStro) <- 'Treatment'
DefaultAssay(PND5.CoDES.MeStro) <- 'ATAC'

###coverage plots of top 15 by DEG pvalue-----

#example
CoveragePlot(
  object = PND5.CoDES.MeStro,
  region = 'Tmem132c',
  extend.upstream = 200000,
  extend.downstream = 200000,
  features = 'Tmem132c',
  expression.assay = 'SCT'

```

Want to use previously generated bulk ChIP data for analysis with my snRNA/ATAC

```{r}
#bring in Wendy data for coordinates ----
library(readxl) #part of tidyverse but not loaded with core
##load ERa ChIP locations----
ERChIP_PND5_CO_table <- read_excel("PND5ER-ChIPseq-Wendy-Sylvia.xlsx", sheet = 'Wendy-HOMER-D5-CO')
ERChIP_PND5_HD_table <- read_excel("PND5ER-ChIPseq-Wendy-Sylvia.xlsx", sheet = 'Wendy-HOMER-D5-DES')
head(ERChIP_PND5_CO_table)
head(ERChIP_PND5_HD_table)

ERChIP_PND5_Co_GR <- makeGRangesFromDataFrame(ERChIP_PND5_CO_table, 
                                               seqnames.field = 'chr', 
                                               start.field = 'start', 
                                               end.field = 'end')
ERChIP_PND5_HD_GR <- makeGRangesFromDataFrame(ERChIP_PND5_HD_table, 
                                               seqnames.field = 'chr', 
                                               start.field = 'start', 
                                               end.field = 'end')

##load H3K27ac ChIP for enhancers----
H3K27ac_PND5_gain_table <- read_excel("PND5_H3K27ac-Enhancer-motif.xlsx", sheet = 'H3K27ac-Enhancer-Gain')
H3K27ac_PND5_loss_table <- read_excel("PND5_H3K27ac-Enhancer-motif.xlsx", sheet = 'H3K27ac-Enhancer-Loss')
head(H3K27ac_PND5_loss_table)

H3K27ac_PND5_loss_GR <- makeGRangesFromDataFrame(H3K27ac_PND5_loss_table, 
                                               seqnames.field = 'chr-peak', 
                                               start.field = 'start-peak', 
                                               end.field = 'end-peak')
H3K27ac_PND5_gain_GR <- makeGRangesFromDataFrame(H3K27ac_PND5_gain_table, 
                                               seqnames.field = 'chr-peak', 
                                               start.field = 'start-peak', 
                                               end.field = 'end-peak')

##determine overlaps between ER binding and H3K27ac marks----
#doesn't really give me what I want because H3K27ac lists are DES gain or loss lists not total lists
ERH3K27_overlap_PND5Co <- findOverlaps(ERChIP_PND5_Co_GR, H3K27ac_PND5_Co_GR, 
                                       minoverlap = 0L, select = 'all')
ERH3K27_overlap_PND5HD <- findOverlaps(ERChIP_PND5_HD_GR, H3K27ac_PND5_HD_GR, 
                                       minoverlap = 0L, select = 'all')

#probably better used instead to match to open chromatin?
#Epithelium
PND5CoDESE.DApeaks.HDtoC <- read.csv('PND5CoDESEpiMeStro_ATAC_DApeaks_HDtoC.csv')
PND5CoDESE.DApeaks.HDtoC <- separate(PND5CoDESE.DApeaks.HDtoC, X, sep = '-', into = c('chr', 'start', 'end'))
PND5CoDESE.DApeaks.HDtoC <- dplyr::filter(PND5CoDESE.DApeaks.HDtoC, p_val_adj <= 0.05)
head(PND5CoDESE.DApeaks.HDtoC)

PND5CoDESE.DApeaks.HDtoC.DESgain <- dplyr::filter(PND5CoDESE.DApeaks.HDtoC, avg_log2FC >= 0.5)
PND5CoDESE.DApeaks.HDtoC.DESloss <- dplyr::filter(PND5CoDESE.DApeaks.HDtoC, avg_log2FC <= -0.5)

head(PND5CoDESE.DApeaks.HDtoC.DESgain)
head(PND5CoDESE.DApeaks.HDtoC.DESloss)

DAPeaks_PND5_HDgain_GR <- makeGRangesFromDataFrame(PND5CoDESE.DApeaks.HDtoC.DESgain,
                                                   seqnames.field = 'chr', 
                                                   start.field = 'start', 
                                                   end.field = 'end')
DAPeaks_PND5_HDloss_GR <- makeGRangesFromDataFrame(PND5CoDESE.DApeaks.HDtoC.DESloss,
                                                   seqnames.field = 'chr', 
                                                   start.field = 'start', 
                                                   end.field = 'end')

##determine overlaps between ER binding or H3K27ac marks v DAPs----
overlapDAPgain_ERaCon_PND5 <- findOverlaps(DAPeaks_PND5_HDgain_GR, ERChIP_PND5_Co_GR, 
                                       minoverlap = 0L, select = 'all')
  #70 overlaps of 5514 accessible peaks gained in HD v 1884 ERa peaks in control
overlapDAPloss_ERaCon_PND5 <- findOverlaps(DAPeaks_PND5_HDloss_GR, ERChIP_PND5_Co_GR, 
                                       minoverlap = 0L, select = 'all')
  #10 overlaps of 902 accessible peaks lost in HD v 1884 ERa peaks in control
overlapDAPgain_ERaHD_PND5 <- findOverlaps(DAPeaks_PND5_HDgain_GR, ERChIP_PND5_HD_GR, 
                                       minoverlap = 0L, select = 'all')
  #647 overlaps of 5514 accessible peaks gained in HD v 5037 ERa peaks in HD
overlapDAPloss_ERaHD_PND5 <- findOverlaps(DAPeaks_PND5_HDloss_GR, ERChIP_PND5_HD_GR, 
                                       minoverlap = 0L, select = 'all')
  #19 overlaps of 902 accessible peaks lost in HD v 5037 ERa peaks in HD

overlapDAPgain_H3K27gain_PND5 <- findOverlaps(DAPeaks_PND5_HDgain_GR, H3K27ac_PND5_gain_GR, 
                                       minoverlap = 0L, select = 'all')
  #558 overlaps of 5514 accessible peaks gained in HD vs 8152 H3K27ac peaks gained in HD
overlapDAPloss_H3K27gain_PND5 <- findOverlaps(DAPeaks_PND5_HDloss_GR, H3K27ac_PND5_gain_GR, 
                                       minoverlap = 0L, select = 'all')
  #14 overlaps of 902 accessible peaks lost in HD vs 8152 H3K27ac peaks gained in HD
overlapDAPgain_H3K27loss_PND5 <- findOverlaps(DAPeaks_PND5_HDgain_GR, H3K27ac_PND5_loss_GR, 
                                       minoverlap = 0L, select = 'all')
  #4 overlaps of 5514 accessible peaks gained in HD vs 4070 H3K27ac peaks lost in HD
overlapDAPloss_H3K27loss_PND5 <- findOverlaps(DAPeaks_PND5_HDloss_GR, H3K27ac_PND5_loss_GR, 
                                       minoverlap = 0L, select = 'all')
  #55 overlaps of 902 accessible peaks lost in HD vs 4070 H3K27ac peaks lost in HD

#MeStro
PND5CoDESMS.DApeaks.HDtoC <- read.csv('PND5CoDESMeStro_ATAC_DApeaks_HDtoC.csv')
PND5CoDESMS.DApeaks.HDtoC <- separate(PND5CoDESMS.DApeaks.HDtoC, X, sep = '-', into = c('chr', 'start', 'end'))
PND5CoDESMS.DApeaks.HDtoC <- dplyr::filter(PND5CoDESMS.DApeaks.HDtoC, p_val_adj <= 0.05)
head(PND5CoDESMS.DApeaks.HDtoC)

PND5CoDESMS.DApeaks.HDtoC.DESgain <- dplyr::filter(PND5CoDESMS.DApeaks.HDtoC, avg_log2FC >= 0.5)
PND5CoDESMS.DApeaks.HDtoC.DESloss <- dplyr::filter(PND5CoDESMS.DApeaks.HDtoC, avg_log2FC <= -0.5)

head(PND5CoDESMS.DApeaks.HDtoC.DESgain)
head(PND5CoDESMS.DApeaks.HDtoC.DESloss)

DAPeaks_PND5MS_HDgain_GR <- makeGRangesFromDataFrame(PND5CoDESMS.DApeaks.HDtoC.DESgain,
                                                   seqnames.field = 'chr', 
                                                   start.field = 'start', 
                                                   end.field = 'end')
DAPeaks_PND5MS_HDloss_GR <- makeGRangesFromDataFrame(PND5CoDESMS.DApeaks.HDtoC.DESloss,
                                                   seqnames.field = 'chr', 
                                                   start.field = 'start', 
                                                   end.field = 'end')
```

Instead, Wendy thinks I should use heatmaps to look for patterns, rather than numbers
To do this, I need to prep and analyze the ChIP data

```{r}
#instead of doing overlap, make ChIP heatmaps at locations determined by ATAC----
#look for patterns rather than numbers

##ERa ChIP----
#the input files (bam-alignments)
bamfile_D5CO="ER-ChIPseq-PND5/BAM/1_030F_00D9NIEHS_d5-Co_ERalpha_mm10_i96.bam";
bamfile_D5DES="ER-ChIPseq-PND5/BAM/2_030G_00D9NIEHS_d5-DES_ERalpha_mm10_i68.bam";
bamfile_D5EpiCo="atac_EpiC.pseudobulkATAC.bam"
bamfile_D5EpiHD="atac_EpiHD.pseudobulkATAC.bam"
bamfile_D5MesCo="atac_MeStroC.pseudobulkATAC.bam"
bamfile_D5MesHD="atac_MeStroHD.pseudobulkATAC.bam"

#parameters
BSgenome="BSgenome.Mmusculus.UCSC.mm10"
uniq=TRUE
extend=0
shift=0
ws=250
standchr = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
             'chr16', 'chr17', 'chr18', 'chr18', 'chrX', 'chrY')

##PREPROCESSING: reading the input files----
D5CO_bulk <- MEDIPS.createSet(file = bamfile_D5CO,BSgenome = BSgenome, extend = extend,
                         shift = shift, window_size = ws, chr.select = standchr)
D5DES_bulk <- MEDIPS.createSet(file = bamfile_D5DES,BSgenome = BSgenome, extend = extend,
                          shift = shift, window_size = ws, chr.select = standchr)
D5CO_epi <- MEDIPS.createSet(file = bamfile_D5EpiCo,BSgenome = BSgenome, extend = extend,
                         shift = shift, window_size = ws, chr.select = standchr)
D5DES_epi <- MEDIPS.createSet(file = bamfile_D5EpiHD,BSgenome = BSgenome, extend = extend,
                          shift = shift, window_size = ws, chr.select = standchr)
D5CO_mes <- MEDIPS.createSet(file = bamfile_D5MesCo,BSgenome = BSgenome, extend = extend,
                         shift = shift, window_size = ws, chr.select = standchr)
D5DES_mes <- MEDIPS.createSet(file = bamfile_D5MesHD,BSgenome = BSgenome, extend = extend,
                          shift = shift, window_size = ws, chr.select = standchr)

#Create the result table and apply the statistical test for differential ERa binding
#bulk Co v DES
mr.edgeR = MEDIPS.meth(MSet1 = D5CO_bulk, MSet2 = D5DES_bulk, p.adj = "fdr", diff.method = "edgeR", MeDIP = F, CNV = F, minRowSum = 1)
head(mr.edgeR)
mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = 0.01, adj = F, ratio = 1.5, bg.counts = NULL, CNV = F)
head(mr.edgeR.s)
write.table(mr.edgeR.s, "ERaMEDIP-no-Input.tsv",sep="\t", row.names = FALSE)	#	Only significant

#epi Con vs bulk DES
epi.C.mr.edgeR = MEDIPS.meth(MSet1 = D5DES_bulk, MSet2 = D5CO_epi, p.adj = "fdr", diff.method = "edgeR", MeDIP = F, CNV = F, minRowSum = 1)
head(epi.C.mr.edgeR)
epi.C.mr.edgeR.s = MEDIPS.selectSig(results = epi.C.mr.edgeR, p.value = 0.01, adj = F, ratio = 1.5, bg.counts = NULL, CNV = F)
head(epi.C.mr.edgeR.s)
write.table(epi.C.mr.edgeR.s, "Epithelium/Unintegrated/ERaMEDIP-EpiConvsBulk.tsv",sep="\t", row.names = FALSE)	#	Only significant

#epi DES vs bulk DES
epi.HD.mr.edgeR = MEDIPS.meth(MSet1 = D5DES_bulk, MSet2 = D5DES_epi, p.adj = "fdr", diff.method = "edgeR", MeDIP = F, CNV = F, minRowSum = 1)
head(epi.HD.mr.edgeR)
epi.HD.mr.edgeR.s = MEDIPS.selectSig(results = epi.HD.mr.edgeR, p.value = 0.01, adj = F, ratio = 1.5, bg.counts = NULL, CNV = F)
head(epi.HD.mr.edgeR.s)
write.table(epi.HD.mr.edgeR.s, "Epithelium/Unintegrated/ERaMEDIP-EpiDESvsBulk.tsv",sep="\t", row.names = FALSE)	#	Only significant

#mes Con vs bulk DES
mes.C.mr.edgeR = MEDIPS.meth(MSet1 = D5DES_bulk, MSet2 = D5CO_mes, p.adj = "fdr", diff.method = "edgeR", MeDIP = F, CNV = F, minRowSum = 1)
head(mes.C.mr.edgeR)
mes.C.mr.edgeR.s = MEDIPS.selectSig(results = mes.C.mr.edgeR, p.value = 0.01, adj = F, ratio = 1.5, bg.counts = NULL, CNV = F)
head(mes.C.mr.edgeR.s)
write.table(mes.C.mr.edgeR.s, "MeStro/Unintegrated/ERaMEDIP-MesConvsBulk.tsv",sep="\t", row.names = FALSE)	#	Only significant

#mes DES vs bulk DES
mes.HD.mr.edgeR = MEDIPS.meth(MSet1 = D5DES_bulk, MSet2 = D5DES_mes, p.adj = "fdr", diff.method = "edgeR", MeDIP = F, CNV = F, minRowSum = 1)
head(mes.HD.mr.edgeR)
mes.HD.mr.edgeR.s = MEDIPS.selectSig(results = mes.HD.mr.edgeR, p.value = 0.01, adj = F, ratio = 1.5, bg.counts = NULL, CNV = F)
head(mes.HD.mr.edgeR.s)
write.table(mes.HD.mr.edgeR.s, "MeStro/Unintegrated/ERaMEDIP-MesDESvsBulk.tsv",sep="\t", row.names = FALSE)	#	Only significant

#merge sets to create a full peak set
Epiallpeaks <- MEDIPS.mergeSets(MSet1 = D5CO_epi, MSet2 = D5DES_epi, name = 'D5epi')
Mesallpeaks <- MEDIPS.mergeSets(MSet1 = D5CO_mes, MSet2 = D5DES_mes, name = 'D5mes')
allpeaks <- MEDIPS.mergeSets(MSet1 = Epiallpeaks, MSet2 = Mesallpeaks, name = 'D5all')

md.Epiallpeaks = MEDIPS.meth(MSet1 = Epiallpeaks, minRowSum = 1)
md.Mesallpeaks = MEDIPS.meth(MSet1 = Mesallpeaks, minRowSum = 1)
md.allpeaks = MEDIPS.meth(MSet1 = allpeaks, minRowSum = 1)

#export as bed
DESspecERaregions <- mr.edgeR.s[,c('chr', 'start', 'stop')]
DESspecERaregions <- rename(DESspecERaregions, 'chrom' = 'chr')
DESspecERaregions <- rename(DESspecERaregions, 'chromStart' = 'start')
DESspecERaregions <- rename(DESspecERaregions, 'chromEnd' = 'stop')
head(DESspecERaregions)
export(DESspecERaregions, "ERaMEDIP-DESspecregions.bed")

mes.C.ERaregions <- mes.C.mr.edgeR.s[,c('chr', 'start', 'stop')]
mes.C.ERaregions <- rename(mes.C.ERaregions, 'chrom' = 'chr')
mes.C.ERaregions <- rename(mes.C.ERaregions, 'chromStart' = 'start')
mes.C.ERaregions <- rename(mes.C.ERaregions, 'chromEnd' = 'stop')
head(mes.C.ERaregions)
export(mes.C.ERaregions, "MeStro/Unintegrated/meConERaRegions.bed")

mes.HD.ERaregions <- mes.HD.mr.edgeR.s[,c('chr', 'start', 'stop')]
mes.HD.ERaregions <- rename(mes.HD.ERaregions, 'chrom' = 'chr')
mes.HD.ERaregions <- rename(mes.HD.ERaregions, 'chromStart' = 'start')
mes.HD.ERaregions <- rename(mes.HD.ERaregions, 'chromEnd' = 'stop')
head(mes.HD.ERaregions)
export(mes.HD.ERaregions, "MeStro/Unintegrated/mesHDERaRegions.bed")

epi.C.ERaregions <- epi.C.mr.edgeR.s[,c('chr', 'start', 'stop')]
epi.C.ERaregions <- rename(epi.C.ERaregions, 'chrom' = 'chr')
epi.C.ERaregions <- rename(epi.C.ERaregions, 'chromStart' = 'start')
epi.C.ERaregions <- rename(epi.C.ERaregions, 'chromEnd' = 'stop')
head(epi.C.ERaregions)
export(epi.C.ERaregions, "Epithelium/Unintegrated/epiConERaRegions.bed")

epi.HD.ERaregions <- epi.HD.mr.edgeR.s[,c('chr', 'start', 'stop')]
epi.HD.ERaregions <- rename(epi.HD.ERaregions, 'chrom' = 'chr')
epi.HD.ERaregions <- rename(epi.HD.ERaregions, 'chromStart' = 'start')
epi.HD.ERaregions <- rename(epi.HD.ERaregions, 'chromEnd' = 'stop')
head(epi.HD.ERaregions)
export(epi.HD.ERaregions, "Epithelium/Unintegrated/epiHDERaRegions.bed")

md.Epiallpeaks <- md.Epiallpeaks[,c('chr', 'start', 'stop')]
md.Epiallpeaks <- rename(md.Epiallpeaks, 'chrom' = 'chr')
md.Epiallpeaks <- rename(md.Epiallpeaks, 'chromStart' = 'start')
md.Epiallpeaks <- rename(md.Epiallpeaks, 'chromEnd' = 'stop')
head(md.Epiallpeaks)
export(md.Epiallpeaks, "Epithelium/Unintegrated/epiallATACRegions.bed")

md.Mesallpeaks <- md.Mesallpeaks[,c('chr', 'start', 'stop')]
md.Mesallpeaks <- rename(md.Mesallpeaks, 'chrom' = 'chr')
md.Mesallpeaks <- rename(md.Mesallpeaks, 'chromStart' = 'start')
md.Mesallpeaks <- rename(md.Mesallpeaks, 'chromEnd' = 'stop')
head(md.Mesallpeaks)
export(md.Mesallpeaks, "MeStro/Unintegrated/MesallATACRegions.bed")

md.allpeaks <- md.allpeaks[,c('chr', 'start', 'stop')]
md.allpeaks <- rename(md.allpeaks, 'chrom' = 'chr')
md.allpeaks <- rename(md.allpeaks, 'chromStart' = 'start')
md.allpeaks <- rename(md.allpeaks, 'chromEnd' = 'stop')
head(md.allpeaks)
export(md.allpeaks, "EpiMeStro/allATACRegions.bed")

##H3K27ac ChIP----
#the input files (bam-alignments)
  #Wendy PND5 H3K27ac (from GO = GSE104402)
awk '{print $0"\t"NR}' H3K27PND5HDgain.BED > H3K27PND5HDgain.bed4

bedToBam -i H3K27PND5HDgain.bed4 -g mouse.mm10.genome > H3K27PND5HDgain.bam

samtools quickcheck H3K27PND5HDgain.bam

bamfile_D5DESdiff_Ac="H3K27PND5HDgain.bam";
  #my ATAC
bamfile_D5EpiCo="Epithelium/Unintegrated/atac_EpiC.pseudobulkATAC.bam"
bamfile_D5EpiHD="Epithelium/Unintegrated/atac_EpiHD.pseudobulkATAC.bam"
bamfile_D5MesCo="MeStro/Unintegrated/atac_MeStroC.pseudobulkATAC.bam"
bamfile_D5MesHD="MeStro/Unintegrated/atac_MeStroHD.pseudobulkATAC.bam"

#differential loop analysis-----
remotes::install_github("EricSDavis/hictoolsr")
install.packages("dbscan")
##read in files----
d5Co_loops <- "PND5COmerged_loops.bedpe"
d5DES_loops <- "PND5DESmerged_loops.bedpe"

d5CoDES_loops <- 
  mergeBedpe(bedpeFiles = c(d5Co_loops, d5DES_loops), res = 10e3) %>%
  as_ginteractions()

head(d5CoDES_loops)

```

Try doing multiomic UMAP with ATAC/RNA with Epithelium & MeStro separately to see if unique populations are revealed

```{r Multiomic UMAP of CoDES Epithelium}
#load in Epithelium ATAC file
PND5.CoDES.Epith <- readRDS("Epithelium/Unintegrated/PND5CoDESEpiMeStro_ATAC_linked_peak.rds")
PND5.CoDES.Epith

##joint neighbor graph with both assays----
PND5.CoDES.Epith <- FindMultiModalNeighbors(PND5.CoDES.Epith, reduction.list = list("pca", "lsi"),
                                  dims.list = list(3:30, 5:30),
                                  modality.weight.name = "SCT.weight",
                                  verbose = TRUE)

##joint UMAP visualization----
PND5.CoDES.Epith <- RunUMAP(PND5.CoDES.Epith, nn.name = 'weighted.nn', assay = 'SCT', verbose = TRUE)

pdf("Epithelium/Unintegrated/PND5CoDESEpith_JointRNAATAC.pdf")
DimPlot(PND5.CoDES.Epith, label = TRUE, repel = TRUE, reduction = 'umap')
dev.off()

pdf("Epithelium/Unintegrated/PND5CoDESEpith_JointRNAATAC_Treatment.pdf")
DimPlot(PND5.CoDES.Epith, label = TRUE, repel = TRUE, group.by ='Treatment', reduction = 'umap')
dev.off()

## Save Seurat of multimodal UMAP----
saveRDS(PND5.CoDES.Epith, file = "Epithelium/Unintegrated/PND5CoDESEpith_JointRNAATAC.Rds", compress = TRUE)

```

```{r Multiomic UMAP of CoDES MeStro}
#load in Epithelium ATAC file
PND5.CoDES.MeStro <- readRDS("MeStro/Unintegrated/PND5CoDESMeStroMeStro_ATAC_linked_peak.rds")
PND5.CoDES.MeStro

##joint neighbor graph with both assays----
PND5.CoDES.MeStro <- FindMultiModalNeighbors(PND5.CoDES.MeStro, reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:20, 3:30),
                                  modality.weight.name = "SCT.weight",
                                  verbose = TRUE)

##joint UMAP visualization----
PND5.CoDES.MeStro <- RunUMAP(PND5.CoDES.MeStro, nn.name = 'weighted.nn', assay = 'SCT', verbose = TRUE)

pdf("MeStro/Unintegrated/PND5CoDESMeStro_JointRNAATAC.pdf")
DimPlot(PND5.CoDES.MeStro, label = TRUE, repel = TRUE, reduction = 'umap')
dev.off()

pdf("MeStro/Unintegrated/PND5CoDESMeStro_JointRNAATAC_Treatment.pdf")
DimPlot(PND5.CoDES.MeStro, label = TRUE, repel = TRUE, group.by ='Treatment', reduction = 'umap')
dev.off()

## Save Seurat of multimodal UMAP----
saveRDS(PND5.CoDES.MeStro, file = "MeStro/Unintegrated/PND5CoDESMeStro_JointRNAATAC.Rds", compress = TRUE)

```
Not enough overlap between control and DES to find integrative-specific subclusters

```{r - shiny app}
#https://github.com/SGDDNB/ShinyCell?tab=readme-ov-file#quick-start-guide
#https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/4cloud.html
library(ShinyCell)

PND5.merge.Epith <- readRDS('Epithelium/Unintegrated/PND5mergedEpith.Rds')
PND5.merge.MeStro <- readRDS('MeStro/Unintegrated/PND5mergedMeStro.Rds')

#app of single 
scConf_MeStro <- createConfig(PND5.merge.MeStro)
makeShinyApp(PND5.merge.MeStro, scConf_MeStro, gene.mapping = TRUE, gex.assay = 'SCT',
             shiny.title = "PND5 Mesenchyme Control v DES")

scConf_Epi <- createConfig(PND5.merge.Epith)
makeShinyApp(PND5.merge.Epith, scConf_Epi, gene.mapping = TRUE, gex.assay = 'SCT',
             shiny.title = "PND5 Epithelium Control v DES")

```

