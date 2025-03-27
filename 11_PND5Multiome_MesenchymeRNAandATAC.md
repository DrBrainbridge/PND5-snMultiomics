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