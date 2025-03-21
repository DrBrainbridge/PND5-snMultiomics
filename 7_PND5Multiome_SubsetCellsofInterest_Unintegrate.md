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
library(scplotter)


library('MEDIPS', lib = '~/R/x86_64-pc-linux-gnu-library/4.3/')
library(rtracklayer)
library(hictoolsr)
library(dbscan)
```
I will split the seruat by exposure and re-analyze 
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