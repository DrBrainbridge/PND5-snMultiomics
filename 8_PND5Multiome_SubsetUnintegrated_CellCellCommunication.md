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
Explore differences in possible cell-cell communication 

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