
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