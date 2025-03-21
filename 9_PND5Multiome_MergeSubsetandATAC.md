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