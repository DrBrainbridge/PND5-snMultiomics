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