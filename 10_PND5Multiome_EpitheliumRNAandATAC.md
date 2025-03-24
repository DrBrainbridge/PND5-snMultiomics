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
PND5.CoDES.Epith <- readRDS("/ddn/gs1/project/nextgen/post/williamsc5/NOV0595-Williams-Bainbridge/Rout_Con/ConDES_scratch/Epithelium/Unintegrated/PND5CoDESEpith.Rds")
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