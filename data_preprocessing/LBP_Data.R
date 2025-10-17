rm(list = ls())
gc()
setwd('Z:/Other computers/ROdongo/PhD Thesis/Data/SYNAPSE - Human/')
# load('FreshMicro.RData')
# load('fm_egdeR.RData')
# source('geTMM.R')
# source("D:/Projects/Common/Differential Gene Expression Analysis/PPDE_Calculator.R")
library(tidyverse)
# library(biomaRt)
library(openxlsx)
library(limma)
library(biomaRt)
library(viridis)
library(ggpubr)
library(AnnotationDbi)
library(org.Hs.eg.db)
pcaPlot <- function(dfnew,Key,title){
  require(ggplot2)
  as.matrix(dfnew)
  pcaANAL <- prcomp(t(dfnew), scale. = T)
  pcaPlotDF <- data.frame(pcaANAL$x[,1:2]) 
  
  ggplot(data = pcaPlotDF, aes(x = PC1, y = PC2, 
                               color = Key))+
    # label = row.names(pcaPlotDF),
    geom_point(aes(size = 0.2))+
    
    scale_color_viridis(discrete = TRUE)+
    # geom_text(aes(label = row.names(pcaPlotDF)),hjust=0.1, vjust=(-0.5), size = 5)+
    # geom_text_repel() +
    labs(title = paste0(title))+
    theme_pubr(base_size = 10, legend = "right")+
    labs_pubr(base_size = 10, base_family = "")+ 
    theme_bw()
}
TMM <- function(
    df,#Raw counts data matrix with Length of each gene in the raw count matrix
    sampleGrps#A factor object with sample-group information
){
  require(edgeR)
  
  countsMat <- as.matrix(df)
  myCPM <- cpm(countsMat, log = F)
  thresh <- myCPM > 0.1
  keep <- rowSums(thresh) >= min(table(sampleGrps))
  # Subset the rows of countdata to keep the more highly expressed genes
  countsMat <- countsMat[keep,]
  return(cpm(calcNormFactors(DGEList(counts=countsMat, group=sampleGrps), method = 'TMM'), log = F)|>
           as.data.frame())
}
degAnalysisLIMMA <- function(Group, df) {
  require(limma)
  design <- model.matrix(~0+Group)
  colnames(design) <- paste0(c('Control','PD'))
  # contrasts <- makeContrasts(Cancer - Control, levels = design)
  # Fit linear models
  # df <- voom(df, design, plot = F)
  fit <- lmFit(log2(df+1), design)
  contrasts <- makeContrasts(PD - Control, levels = colnames(coef(fit)))
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit)
  return(topTable(fit, number = Inf)|> 
           as.data.frame()|>
           dplyr::select(logFC,	P.Value,	adj.P.Val)|>
           rownames_to_column(var = 'GeneID')|>
           dplyr::mutate(FC = (2^logFC)))
}
geTMM <- function(
    df,#Raw counts data matrix with Length of each gene in the raw count matrix
    sampleGrps#A factor object with sample-group information
){
  require(edgeR)

  countsMat <- as.matrix(df |> dplyr::select(-Length))
  myCPM <- cpm(countsMat, log = F)
  thresh <- myCPM > 0.1
  keep <- rowSums(thresh) >= min(table(sampleGrps))
  # Subset the rows of countdata to keep the more highly expressed genes
  countsMat <- countsMat[keep,]
  geneLength <- as.numeric(df$Length[keep])
  rpk <- countsMat/(geneLength)
  return(cpm(calcNormFactors(DGEList(counts=rpk, group=sampleGrps), method = 'TMM'), log = F)|>
           as.data.frame())
}
geneLen.ROSMAP <- openxlsx::read.xlsx("Z:/Other computers/ROdongo/Benchmarking Study/ExonLengths_GRCh38.v37.xlsx")|>
  dplyr::rename(GeneID = 1, Length = 2)
ensembl <- read_csv("gProfiler_hsapiens_convert.csv")|>
  dplyr::select(initial_alias,converted_alias)|>
  dplyr::filter(converted_alias != 'None')|>
  distinct()
tx2gene <- read_delim("Z:/Other computers/ROdongo/LUHMES Study/Kallisto v2/transcripts_to_genes.txt",
                      delim = "\t", escape_double = FALSE,
                      col_names = FALSE, trim_ws = TRUE)|>
  dplyr::rename(TxID = 1, GeneID = 2, SYMBOL=3)|>
  # dplyr::left_join(ensembl,by=c('SYMBOL'='initial_alias'))|>
  distinct()|>
  # dplyr::filter(converted_alias != 'none')|>
  dplyr::mutate(ENSEMBL = sapply(GeneID, function(x){str_split_i(x,'\\.',1)}))

# openxlsx::write.xlsx(tx2gene,'tx2gene.xlsx')
# Living Brain Project

# LBP_assay_RNAseq_metadata <- read_csv("D:/Projects/PhD Thesis/Data/SYNAPSE - Human/syn26337520/LBP_assay_RNAseq_metadata.csv")
# LBP_biospecimen_metadata <- read_csv("D:/Projects/PhD Thesis/Data/SYNAPSE - Human/syn26337520/LBP_biospecimen_metadata.csv")
covAdj_v2 <- function(metaData, exprDATA) {
  metaData <- metaData|>
    dplyr::filter(specimenID %in% colnames(exprDATA))
  exprDATA <- exprDATA |>
    as.data.frame()|>
    # column_to_rownames(var = 'GeneID')|>
    dplyr::select(metaData$specimenID)|>
    as.matrix()
  # Calculate regression coefficients
  lmCoeff <- lm(t(exprDATA) ~ metaData$CONDITION + metaData$BATCH + metaData$PM + 
                  metaData$SEX + metaData$AGE+ metaData$ETHNICITY + metaData$TISSUE)$coefficients |>
    t() |>
    as.data.frame() |>
    dplyr::rename(INTERCEPT = 1,CONDITION=2,BATCH=3,PM=4,SEX=5,AGE = 6,ETHNICITY=7,TISSUE=8)|>
    dplyr::select(BATCH,PM,SEX,AGE,ETHNICITY,TISSUE)|>
    as.matrix()
  fitsv <- lm.fit(lmCoeff, exprDATA)
  
  fitsv <- fitsv$residuals
  # return(df1 + abs(min(df1)) + 1)
  return(fitsv+abs(min(fitsv))+1)
}
gtf2 <- rtracklayer::readGFF('Z:/Other computers/ROdongo/PhD Thesis/Homo_sapiens.GRCh38.114.gtf.gz')|>
  group_by(gene_id)|>
  dplyr::summarise(gene_biotype = paste0(unique(gene_biotype)),
                   gene_name = paste0(unique(gene_name)))|>
  dplyr::filter(gene_biotype=='protein_coding')
LBP_FlagshipPaper_Metadata <- read_csv("syn26337520/LBP_FlagshipPaper_Metadata.csv")
# |>
#   dplyr::select(specimenID,mymet_postmortem,mymet_tissue,mymet_phe,mymet_age,mymet_sex,mymet_ethnicity,Bank)|>
#   dplyr::rename(condition = 4)|>
#   dplyr::mutate(AGE = as.numeric(mymet_age),
#                 ETHNICITY = as.numeric(as.factor(mymet_ethnicity)),
#                 PM = as.numeric(as.factor(mymet_postmortem)),
#                 SEX = as.numeric(as.factor(mymet_sex)),
#                 TISSUE = as.numeric(as.factor(mymet_tissue)),
#                 BATCH = as.numeric(as.factor(Bank)),
#                 DIAGNOSIS = as.numeric(as.factor(condition)))
# table(LBP_FlagshipPaper_Metadata$DIAGNOSIS)
# table(LBP_FlagshipPaper_Metadata$ETHNICITY)
# table(LBP_FlagshipPaper_Metadata$PM)

individualID <- LBP_FlagshipPaper_Metadata|>
  dplyr::group_by(individualID)|>
  dplyr::summarise(Bank = paste0(unique(Bank),collapse = ','),
                   Brain.Biobank = paste0(unique(mymet_bank),collapse = ','),
                   Brain.Axis = paste0(unique(mymet_tissue),collapse = ','),
                   DIAGNOSIS = paste0(unique(mymet_phe),collapse = ','),
                   ETHNICITY = paste0(unique(mymet_ethnicity),collapse = ','),
                   SEX = paste0(unique(mymet_sex),collapse = ','),
                   BATCH = paste0(unique(mymet_extractionbatch),collapse = ','),
                   AGE = paste0(unique(mymet_age),collapse = ','))
idx <- sapply(unique(LBP_FlagshipPaper_Metadata$individualID), function(x){
  df <- LBP_FlagshipPaper_Metadata|>
    dplyr::filter(individualID==x)
  if(('R_Brain'%in%df$mymet_tissue) & ('L_Brain'%in%df$mymet_tissue)){
    idx1 <- df$specimenID[df$mymet_tissue=='L_Brain']
  }
  else{
    idx1 <- df$specimenID
  }
  return(idx1)
})

LBP_FlagshipPaper_Metadata.1 <- LBP_FlagshipPaper_Metadata|>
  dplyr::filter(specimenID %in% idx)|>
  dplyr::mutate(Bank = Bank,
                Brain.Biobank = mymet_bank,
                TISSUE = as.numeric(as.factor(mymet_tissue)),
                DIAGNOSIS = as.numeric(as.factor(mymet_phe)),
                PM = as.numeric(mymet_postmortem),
                ETHNICITY = as.numeric(as.factor(mymet_ethnicity)),
                SEX = as.numeric(as.factor(mymet_sex)),
                BATCH = as.numeric(as.factor(mymet_extractionbatch)),
                AGE = as.numeric(mymet_age))

# openxlsx::write.xlsx(list(Summary.Phenotype = table(individualID$Phenotype) |> as.data.frame()|>dplyr::rename(Condition=1,Participants=2),
#                           # Summary.Living = individualID |> dplyr::group_by(Bank)|>dplyr::summarise(Phenotype_n=n_distinct(Phenotype),Sex_n = n_distinct(Sex)),
#                           Metadata = individualID),
#                      file = 'syn26337520/LBP_Summarised.xlsx',asTable = T)
individualID$Brain.Axis <- LBP_FlagshipPaper_Metadata.1$mymet_tissue[match(LBP_FlagshipPaper_Metadata.1$individualID,individualID$individualID)]
individualID$specimenID <- LBP_FlagshipPaper_Metadata.1$specimenID[match(LBP_FlagshipPaper_Metadata.1$individualID,individualID$individualID)]
individualID$BATCH <- LBP_FlagshipPaper_Metadata.1$BATCH[match(LBP_FlagshipPaper_Metadata.1$individualID,individualID$individualID)]
individualID$PM <- LBP_FlagshipPaper_Metadata.1$PM[match(LBP_FlagshipPaper_Metadata.1$individualID,individualID$individualID)]
# LBP_individual_metadata <- read_csv("D:/Projects/PhD Thesis/Data/SYNAPSE - Human/syn26337520/LBP_individual_metadata.csv")

LBP_FlagshipPaper_featureCounts <- read_csv("./syn26337520/LBP_FlagshipPaper_featureCounts.csv")|>
  dplyr::rename(GeneID = 1)|>
  dplyr::select(GeneID, individualID$specimenID)|>
  dplyr::mutate(GeneID = sapply(GeneID, function(x){str_split_i(x,'\\.',1)}))|>
  # dplyr::left_join(geneLen.ROSMAP, by = 'GeneID')|>
  distinct()|>
  drop_na(GeneID)|>
  dplyr::filter(GeneID %in% gtf2$gene_id)|>
  dplyr::select(GeneID,individualID$specimenID)
# gtf <- rtracklayer::readGFF("D:/Projects/PhD Thesis/Homo_sapiens.GRCh38.114.gtf.gz")

# gtf2 <- gtf |> dplyr::mutate(GeneID = paste0(gene_id,'.',gene_version))|>
#   dplyr::filter(gene_biotype=='protein_coding')
# ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# proteinCoding <- getBM(attributes=c('ensembl_gene_id', 'ensembl_gene_id_version','hgnc_symbol','gene_biotype'), 
#                        filters = 'ensembl_gene_id_version', values =LBP_FlagshipPaper_featureCounts$GeneID, mart = ensembl)|>
#   as.data.frame()|>
#   dplyr::filter(gene_biotype=='protein_coding')

dim(LBP_FlagshipPaper_featureCounts)


variance_cols <- colnames(LBP_FlagshipPaper_featureCounts)[-1]
LBP_FlagshipPaper_featureCounts <- LBP_FlagshipPaper_featureCounts %>%
  mutate(row_variance = apply(.[, variance_cols], 1, var)) %>%
  group_by(GeneID) %>%
  slice_max(order_by = row_variance, n = 1, with_ties = F) %>%
  ungroup()|>
  dplyr::select(-row_variance)|>
  tidyr::drop_na()|>
  dplyr::mutate(Length = geneLen.ROSMAP$Length[match(GeneID,geneLen.ROSMAP$GeneID)]) |>
  remove_rownames()|>column_to_rownames(var = 'GeneID')|>
  geTMM(sampleGrps = individualID$DIAGNOSIS)

openxlsx::write.xlsx(LBP_FlagshipPaper_featureCounts|>rownames_to_column(var = 'GeneID'),'LBP_geTMM_Normalized.Xlsx')
dim(LBP_FlagshipPaper_featureCounts)
KEY <- individualID$DIAGNOSIS
print(pcaPlot(LBP_FlagshipPaper_featureCounts,KEY,'Living Brain Project - RNA Seq (before covariate adjustment'))

ggsave(paste0('./syn26337520/PCA Plots/LBP_RNASeq_beforeCovAdj_geTMM.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')

individualID <- individualID|>
  dplyr::mutate(TISSUE=as.numeric(as.factor(Brain.Axis)),
                CONDITION=as.numeric(as.factor(DIAGNOSIS)),
                ETHNICITY=as.numeric(as.factor(ETHNICITY)),
                SEX=as.numeric(as.factor(SEX)),
                AGE=as.numeric(AGE))
LBP_FlagshipPaper_featureCounts.covAdj <- covAdj_v2(metaData = individualID, exprDATA = LBP_FlagshipPaper_featureCounts)

max(LBP_FlagshipPaper_featureCounts.covAdj);min(LBP_FlagshipPaper_featureCounts.covAdj)
# exprDATA_CovAdj <- exprDATA_CovAdj + abs(exprDATA_CovAdj)+ 1
# max(exprDATA_CovAdj);min(exprDATA_CovAdj)
print(pcaPlot(LBP_FlagshipPaper_featureCounts.covAdj,KEY,'Living Brain Project - RNA Seq (after covariate adjustment'))
ggsave(paste0('./syn26337520/PCA Plots/LBP_RNASeq_afterCovAdj_geTMM.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')

openxlsx::write.xlsx(list(CONTROL = LBP_FlagshipPaper_featureCounts.covAdj |> 
                            as.data.frame()|>
                            dplyr::select(individualID$specimenID[which(individualID$DIAGNOSIS=='Control')])|>
                            rownames_to_column(var = 'GeneID'), 
                          AD = LBP_FlagshipPaper_featureCounts.covAdj |> 
                            as.data.frame()|>
                            dplyr::select(individualID$specimenID[which(individualID$DIAGNOSIS=='PD')])|>
                            rownames_to_column(var = 'GeneID')),
                     './syn26337520/Expression Data/LBP_geTMM_final.xlsx')

degLBP <- degAnalysisLIMMA(Group = individualID$DIAGNOSIS, 
                            df = LBP_FlagshipPaper_featureCounts.covAdj|> as.data.frame()|>dplyr::select(individualID$specimenID))

openxlsx::write.xlsx(degLBP, './syn26337520/DEG Data/LBP_DEG_geTMM_final.xlsx')


# # Covariate adjustment
# # ATAC-Seq
# covAdj <- function(metaData, df) {
#   metaData <- metaData|>
#     dplyr::filter(specimenID %in% colnames(df))
#   df <- df |>
#     # remove_rownames()|>
#     # column_to_rownames(var = 'GeneID')|>
#     dplyr::select(metaData$specimenID)
#   # Calculate regression coefficients
#   lmCoeff <- lm(t(df) ~ metaData$CONDITION + metaData$BATCH + metaData$PM + 
#                   metaData$SEX + metaData$AGE+ metaData$ETHNICITY + metaData$TISSUE)$coefficients |>
#     t() |>
#     as.data.frame() |>
#     dplyr::rename(INTERCEPT = 1,CONDITION=2,BATCH=3,PM=4,SEX=5,AGE = 6,ETHNICITY=7,TISSUE=8)
#   df1 <- df
#   # Subtract the effect of each covariate
#   for(i in 1:nrow(metaData)){
#     # i <- 1
#     adj.Samp <- metaData$specimenID[i]
#     # Remove condition!
#     adj.Vals <- rowSums(cbind(lmCoeff$BATCH*metaData$BATCH[i],lmCoeff$PM*metaData$PM[i],
#                               lmCoeff$SEX*metaData$SEX[i],lmCoeff$AGE*metaData$AGE[i],
#                               lmCoeff$ETHNICITY*metaData$ETHNICITY[i],lmCoeff$TISSUE*metaData$TISSUE[i]))
#     adj.Vec <- which(colnames(df1)==adj.Samp)
#     df1[,adj.Vec] <- df1[,adj.Vec]-adj.Vals
#   }
#   # return(df1 + abs(min(df1)) + 1)
#   return(df1 + abs(min(df1)) + 1)
# }
# 
# LBP_FlagshipPaper_featureCounts_CovAdj <- covAdj_v2(individualID,LBP_FlagshipPaper_featureCounts)
# KEY <- individualID$DIAGNOSIS
# print(pcaPlot(LBP_FlagshipPaper_featureCounts_CovAdj,KEY,'Living Brain Project - RNA Seq (after covariate adjustment)'))
# 
# ggsave(paste0('LBP_RNASeq_afterCovAdj.tiff'),
#        units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
# # LBP_FlagshipPaper_featureCounts_CovAdj <- LBP_FlagshipPaper_featureCounts_CovAdj + abs(LBP_FlagshipPaper_featureCounts_CovAdj)+1
# min(LBP_FlagshipPaper_featureCounts_CovAdj)
# 
# 
# 
# # Mitochondrial genes
# # Human_MitoCarta2_0 <- openxlsx::read.xlsx("D:/Projects/LUHMES Study/Human MitoCarta2.0.xlsx")
# # tx2gene1 <- tx2gene[!duplicated(tx2gene$GeneID),]
# # ddsDEG_1 <- do.call(cbind,ddsDEG)|> as.data.frame() |> 
# #   dplyr::mutate(GeneID = ddsDEG@rownames,
# #                 GeneID = sapply(GeneID, function(x){str_split_i(x,'_',1)}),
# #                 FC = 2^as.numeric(log2FoldChange))|>
# #   dplyr::select(GeneID, FC,log2FoldChange,pvalue,padj)|>
# #   dplyr::rename(GeneID = 1, logFC = 3,PValue=4, FDR = 5)|>
# #   drop_na(FDR)|>
# #   left_join(tx2gene1[,2:3], by = 'GeneID')|>
# #   # dplyr::select(-GeneID)|>
# #   dplyr::rename(Symbol = 6)|>
# #   dplyr::select(GeneID,Symbol,logFC,FC,PValue,FDR)
# 
# # Create a design matrix (replace 'Group' with your group labels)
# degAnalysisLIMMA <- function(Group, df) {
#   require(limma)
#   design <- model.matrix(~0+Group)
#   # ~0 + 
#   colnames(design) <- paste0(c('PD','Control'))
#   contrasts <- makeContrasts(PD - Control, levels = design)
#   # Fit linear models
#   # df <- voom(df, design, plot = F)
#   fit <- lmFit(log2(df+1), design)
#   fit <- contrasts.fit(fit, contrasts)
#   fit <- eBayes(fit)
#   return(topTable(fit, number = Inf)|> 
#            as.data.frame()|>
#            dplyr::select(logFC,	P.Value,	adj.P.Val)|>
#            rownames_to_column(var = 'GeneID'))
# }
# # RNA-Seq
# 
# control <- LBP_FlagshipPaper_featureCounts_CovAdj|>
#   dplyr::select(individualID$specimenID[which(individualID$DIAGNOSIS=='Control')])
# control$AverageCOntrol <- apply(control,1,mean)
# PD <- LBP_FlagshipPaper_featureCounts_CovAdj|>
#   dplyr::select(individualID$specimenID[which(individualID$DIAGNOSIS=='PD')])
# PD$AveragePD <- apply(PD,1,mean)
# degData <- cbind(control,PD)|>
#   dplyr::mutate(FC = AveragePD/AverageCOntrol)|>
#   rownames_to_column(var = 'GeneID')|>
#   dplyr::select(GeneID,FC)
# degRNASeq <- degAnalysisLIMMA(Group = individualID$DIAGNOSIS, 
#                               df = LBP_FlagshipPaper_featureCounts_CovAdj)|>
#   left_join(degData,by='GeneID')
# # |>
# #   dplyr::mutate(FC=apply(LBP_FlagshipPaper_featureCounts_CovAdj[,individualID$specimenID[which(individualID$DIAGNOSIS=='PD')]],1,mean)/apply(LBP_FlagshipPaper_featureCounts_CovAdj[,individualID$specimenID[which(individualID$DIAGNOSIS=='Control')]],1,mean),
# #                 logFC2 = log2(FC))
# openxlsx::write.xlsx(list(RNASeq = degRNASeq),
#                      file = './syn26337520/DEG Data/LBP_DEG_TMM.xlsx')
rm(list = ls())
gc()
