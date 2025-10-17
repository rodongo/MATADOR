rm(list = ls())
setwd('../ROSMAP/')
# setwd('Z:/Other computers/ROdongo/ROSMAP')
load('ROSMAP-Neff2021.RData')

library(tidyverse)
library(openxlsx)
library(limma)
library(biomaRt)
library(viridis)
library(ggpubr)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load data
ROSMAP_ensmbl_unique <- openxlsx::read.xlsx("../ROSMAP_ensmbl_unique.xlsx") |>
  dplyr::select(-geneSymbol)|>
  column_to_rownames(var='EnsemblID')

metadatarosmap <- openxlsx::read.xlsx("../metadatarosmap.xlsx")|>
  dplyr::mutate(projid = sapply(projid, function(x){paste0('Proj',x,collapse = '')}))|>
  dplyr::filter(projid %in% colnames(ROSMAP_ensmbl_unique))

# Check the sample ids
length(intersect(colnames(ROSMAP_ensmbl_unique),metadatarosmap$projid))

# ROSMAP_ensmbl_unique <- ROSMAP_ensmbl_unique |> dplyr::select(contains(metadatarosmap$projid))
# Save expression data
openxlsx::write.xlsx(list(CONTROL = ROSMAP_ensmbl_unique |> 
                            dplyr::select(metadatarosmap$projid[which(metadatarosmap$condition=='control')])|>
                            rownames_to_column(var = 'GeneID'), 
                          AD = ROSMAP_ensmbl_unique |> 
                            dplyr::select(metadatarosmap$projid[which(metadatarosmap$condition=='AD')])|>
                            rownames_to_column(var = 'GeneID')),
                     './Expression Data/ROSMAP.xlsx')
# DEG analysis
degAnalysisLIMMA <- function(Group, df) {
  require(limma)
  design <- model.matrix(~0+Group)
  # ~0 + 
  colnames(design) <- paste0(c('AD','Control'))
  contrasts <- makeContrasts(AD - Control, levels = design)
  # Fit linear models
  # df <- voom(df, design, plot = F)
  fit <- lmFit(df, design)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit)
  return(topTable(fit, number = Inf)|> 
           as.data.frame()|>
           dplyr::select(logFC,	P.Value,	adj.P.Val)|>
           rownames_to_column(var = 'GeneID'))
}
# RNA-Seq
ROSMAP_ensmbl_unique <- ROSMAP_ensmbl_unique |> dplyr::select(contains(metadatarosmap$projid))
ROSMAP_ensmbl_unique.1 <- ROSMAP_ensmbl_unique + abs(min(ROSMAP_ensmbl_unique))+1
degRNASeq <- degAnalysisLIMMA(Group = metadatarosmap$condition, 
                              df = log2(ROSMAP_ensmbl_unique.1))
Control <- ROSMAP_ensmbl_unique |> dplyr::select(as.character(metadatarosmap$projid[which(metadatarosmap$condition=='control')]))
AD <- ROSMAP_ensmbl_unique |> dplyr::select(as.character(metadatarosmap$projid[which(metadatarosmap$condition=='AD')]))

degFC <- data.frame(FC = apply(AD,1,mean)/apply(Control,1,mean)) |>
  rownames_to_column(var = 'GeneID')

degRNASeq <- degRNASeq |>
  dplyr::left_join(degFC, by='GeneID')
openxlsx::write.xlsx(list(RNASeq = degRNASeq),
                     file = './DEG Data/ROSMAP_DEG.xlsx')
save.image('ROSMAP-Neff2021.RData')
