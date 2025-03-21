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