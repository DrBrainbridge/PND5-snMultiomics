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