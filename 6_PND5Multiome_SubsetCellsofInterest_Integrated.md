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
