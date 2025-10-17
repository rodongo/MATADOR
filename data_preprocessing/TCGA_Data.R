
setwd('../Cancer Data/')
rm(list = ls())
# BiocManager::install('TCGAbiolinks')
# BiocManager::install('TCGAWorkflow')
# BiocManager::install('TCGAWorkflowData')
# BiocManager::install('SummarizedExperiment')
# library(TCGAbiolinks)
# library(TCGAWorkflow)
# library(TCGAWorkflowData)
# library(DT)
library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(viridis)
library(ggpubr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
# ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# TCGAbiolinks::getGDCprojects()
# Create a design matrix (replace 'Group' with your group labels)
degAnalysisLIMMA <- function(Group, df) {
  require(limma)
  design <- model.matrix(~0+Group)
  colnames(design) <- paste0(c('Cancer','Control'))
  # contrasts <- makeContrasts(Cancer - Control, levels = design)
  # Fit linear models
  # df <- voom(df, design, plot = F)
  fit <- lmFit(log2(df+1), design)
  contrasts <- makeContrasts(Cancer - Control, levels = colnames(coef(fit)))
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit)
  return(topTable(fit, number = Inf)|> 
           as.data.frame()|>
           dplyr::select(logFC,	P.Value,	adj.P.Val)|>
           rownames_to_column(var = 'GeneID')|>
           dplyr::mutate(FC = (2^logFC)))
}
covAdj <- function(metaData, df) {
  metaData <- metaData|>
    dplyr::filter(specimenID %in% colnames(df))
  df <- df |>
    remove_rownames()|>
    column_to_rownames(var = 'GeneID')|>
    dplyr::select(metaData$specimenID)
  # Calculate regression coefficients
  lmCoeff <- lm(t(df) ~ metaData$DIAGNOSIS + metaData$AGE)$coefficients |>
    t() |>
    as.data.frame() |>
    dplyr::rename(INTERCEPT = 1,DIAGNOSIS=2,AGE = 3)
  df1 <- df
  # Subtract the effect of each covariate
  for(i in 1:nrow(metaData)){
    # i <- 1
    adj.Samp <- metaData$specimenID[i]
    # Remove condition!
    # adj.Vals <- rowSums(cbind(lmCoeff$AGE*metaData$AGE[i]))
    adj.Vec <- which(colnames(df1)==adj.Samp)
    df1[,adj.Vec] <- df1[,adj.Vec]-lmCoeff$AGE*metaData$AGE[i]
  }
  return(df1 + abs(min(df1)) + 1)
  # return(df1)
}
covAdj_v2 <- function(metaData, exprDATA) {
  metaData <- metaData|>
    dplyr::filter(specimenID %in% colnames(exprDATA))
  exprDATA <- exprDATA |>
    as.data.frame()|>
    # column_to_rownames(var = 'GeneID')|>
    dplyr::select(metaData$specimenID)|>
    as.matrix()
  # Calculate regression coefficients
  lmCoeff <- lm(t(exprDATA) ~ metaData$DIAGNOSIS + metaData$AGE)$coefficients |>
    t() |>
    as.data.frame() |>
    dplyr::rename(INTERCEPT = 1,DIAGNOSIS=2,AGE = 3)|>
    dplyr::select(AGE)|>
    as.matrix()
  fitsv <- lm.fit(lmCoeff, exprDATA)
  fitsv <- fitsv$residuals
  return(fitsv+abs(min(fitsv))+1)
}
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

# Load GTF file to extract protein coding gene information
gtf2 <- rtracklayer::readGFF('Z:/Other computers/ROdongo/PhD Thesis/Homo_sapiens.GRCh38.114.gtf.gz')|>
  group_by(gene_id)|>
  dplyr::summarise(gene_biotype = paste0(unique(gene_biotype)),
                   gene_name = paste0(unique(gene_name)))|>
  dplyr::filter(gene_biotype=='protein_coding')

geTMM <- function(
    df,#Raw counts data matrix with Length of each gene in the raw count matrix
    sampleGrps#A factor object with sample-group information
){
  require(edgeR)
  countsMat <- as.matrix(df |> dplyr::select(-LENGTH))
  myCPM <- cpm(countsMat, log = F)
  thresh <- myCPM > 0.1
  keep <- rowSums(thresh) >= min(table(sampleGrps))
  # Subset the rows of countdata to keep the more highly expressed genes
  countsMat <- countsMat[keep,]
  geneLength <- as.numeric(df$LENGTH[keep])
  rpk <- countsMat/(geneLength/1000)
  return(cpm(calcNormFactors(DGEList(counts=rpk, group=sampleGrps), method = 'TMM'), log = T)|>
           as.data.frame())
}
ExonLengths_GRCh38_v37 <- openxlsx::read.xlsx("../ExonLengths_GRCh38.v37.xlsx")|>
  dplyr::rename(GeneID = 1, LENGTH = 2)
dim(ExonLengths_GRCh38_v37)
ExonLengths_GRCh38_v37 <- ExonLengths_GRCh38_v37[!duplicated(ExonLengths_GRCh38_v37$GeneID),]
dim(ExonLengths_GRCh38_v37)

# Access, download, normalize and adjust data for covariates
# Lung Adenocarcinoma (LUAD)
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  # platform = "Illumina HiSeq", 
  workflow.type = "STAR - Counts",
  # file.type = "results",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)
gc()
LUAD.data <- GDCprepare(query,
                        save = TRUE,
                        save.filename = "TCGA_LUAD.rda",
                        remove.files.prepared = TRUE
)
tcga_LUAD_metadata <- data.frame(LUAD.data$patient,LUAD.data$sample,
                                 LUAD.data$tumor_descriptor,LUAD.data$definition,
                                 LUAD.data$barcode, LUAD.data$days_to_diagnosis,
                                 LUAD.data$days_to_last_follow_up,LUAD.data$age_at_diagnosis,
                                 LUAD.data$primary_diagnosis,LUAD.data$diagnosis_id)

data2 <- LUAD.data@assays@data@listData$unstranded %>% 
  as.data.frame()|>
  dplyr::rename_with(~LUAD.data@colData@rownames)|>
  dplyr::mutate(GeneID = LUAD.data@rowRanges@elementMetadata$gene_id,
                GeneName = LUAD.data@rowRanges@elementMetadata$gene_name,
                HGNC = LUAD.data@rowRanges@elementMetadata$hgnc_id)
saveRDS(data2, 'LUAD_TCGA.RDS')
gc()
tcga_LUAD_metadata <- tcga_LUAD_metadata |>
  as.data.frame()|>
  dplyr::rename(SOURCE = 4, ID = 5)|>
  dplyr::mutate(DIAGNOSIS = case_when(SOURCE == 'Solid Tissue Normal' ~ 'Control',
                                      SOURCE == 'Primary solid Tumor' ~ 'Cancer',
                                      .default = SOURCE))|>
  dplyr::filter(DIAGNOSIS == 'Control'|DIAGNOSIS=='Cancer')|>
  dplyr::rename(specimenID = 5,
                AGE = 8)
gc()
saveRDS(tcga_LUAD_metadata,'tcga_LUAD_metadata.RDS')

tcga_LUAD_metadata <- readRDS('tcga_LUAD_metadata.RDS')|>
  dplyr::rename(specimenID = 5,
                AGE = 8)|>
  drop_na(AGE)
exprDATA <- readRDS('LUAD_TCGA.RDS') |>
  dplyr::select(-c(HGNC,GeneName))|>
  drop_na(GeneID)|>
  distinct()|>
  dplyr::mutate(GeneID2 = sapply(GeneID, function(x){str_split_i(x,'\\.',1)}))|>
  dplyr::filter(GeneID2 %in% gtf2$gene_id)|>
  dplyr::select(-GeneID)|>
  dplyr::rename(GeneID = GeneID2)|>
  dplyr::select(GeneID,tcga_LUAD_metadata$specimenID)

dim(exprDATA) # [1] 32610   601
variance_cols <- colnames(exprDATA)[-1]
exprDATA <- exprDATA %>%
  mutate(row_variance = apply(.[, variance_cols], 1, var)) %>%
  group_by(GeneID) %>%
  slice_max(order_by = row_variance, n = 1, with_ties = F) %>%
  ungroup()|>
  dplyr::select(-row_variance)|>
  tidyr::drop_na()|>
  remove_rownames()|>column_to_rownames(var = 'GeneID')|>
  dplyr::select(tcga_LUAD_metadata$specimenID)|>
  TMM(tcga_LUAD_metadata$DIAGNOSIS)

dim(exprDATA) # [1] 20976   560
KEY <- tcga_LUAD_metadata$DIAGNOSIS[match(colnames(exprDATA),tcga_LUAD_metadata$specimenID)] 
# Cancer Control 
# 508      52
print(table(KEY))
print(pcaPlot(exprDATA,KEY,'LUAD'))
ggsave(paste0('./TMM/PCA Plots/LUAD_beforeCovAdj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
# exprDATA_CovAdj <- covAdj(tcga_LUAD_metadata,exprDATA)
exprDATA_CovAdj <- covAdj_v2(metaData = tcga_LUAD_metadata,exprDATA = exprDATA)

max(exprDATA_CovAdj);min(exprDATA_CovAdj)
# exprDATA_CovAdj <- exprDATA_CovAdj + abs(exprDATA_CovAdj)+ 1
# max(exprDATA_CovAdj);min(exprDATA_CovAdj)
print(pcaPlot(exprDATA_CovAdj,KEY,'LUAD'))
ggsave(paste0('./TMM/PCA Plots/LUAD_afterCovAdj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
# DEG

degLUAD <- degAnalysisLIMMA(Group = tcga_LUAD_metadata$DIAGNOSIS, 
                            df = exprDATA_CovAdj|> as.data.frame()|>dplyr::select(tcga_LUAD_metadata$specimenID))

openxlsx::write.xlsx(degLUAD, '../DEG Data/LUAD_DEG.xlsx')
df <- list()
df[['Control']] <- exprDATA_CovAdj|> as.data.frame()|> dplyr::select(tcga_LUAD_metadata$specimenID[which(tcga_LUAD_metadata$DIAGNOSIS == 'Control')])|> rownames_to_column(var = 'GeneID')
df[['Cancer']] <- exprDATA_CovAdj|> as.data.frame()|> dplyr::select(tcga_LUAD_metadata$specimenID[which(tcga_LUAD_metadata$DIAGNOSIS == 'Cancer')])|> rownames_to_column(var = 'GeneID')
openxlsx::write.xlsx(df, '../Expression Data/LUAD.xlsx')
gc()
# 
# BRCA
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)
BRCA.data <- GDCprepare(query,
                        save = TRUE,
                        save.filename = "TCGA_BRCA.rda",
                        remove.files.prepared = TRUE
)
tcga_BRCA_metadata <- data.frame(BRCA.data$patient,BRCA.data$sample,
                                  BRCA.data$tumor_descriptor,BRCA.data$definition,
                                  BRCA.data$barcode, BRCA.data$days_to_diagnosis,
                                  BRCA.data$days_to_last_follow_up,BRCA.data$age_at_diagnosis,
                                  BRCA.data$primary_diagnosis,BRCA.data$diagnosis_id)

saveRDS(tcga_BRCA_metadata,'tcga_BRCA_metadata.RDS')
data2 <- BRCA.data@assays@data@listData$unstranded |>
  as.data.frame()|>
  dplyr::rename_with(~BRCA.data@colData@rownames)|>
  dplyr::mutate(GeneID = BRCA.data@rowRanges@elementMetadata$gene_id,
                GeneName = BRCA.data@rowRanges@elementMetadata$gene_name,
                HGNC = BRCA.data@rowRanges@elementMetadata$hgnc_id)
# colnames(data2)[1:431] <- paste0(BRCA.data@colData@rownames)
saveRDS(data2, 'TCGA_BRCA.RDS')
rm(list = ls())

tcga_BRCA_metadata <- readRDS('tcga_BRCA_metadata.RDS')|>
  as.data.frame()|>
  as.data.frame()|>
  dplyr::rename(SOURCE = 4, ID = 5)|>
  dplyr::mutate(DIAGNOSIS = case_when(SOURCE == 'Solid Tissue Normal' ~ 'Control',
                                      SOURCE == 'Primary solid Tumor' ~ 'Cancer',
                                      .default = SOURCE))|>
  dplyr::filter(DIAGNOSIS == 'Control'|DIAGNOSIS=='Cancer')|>
  dplyr::rename(specimenID = 5,
                AGE = 8)|>
  drop_na(AGE)
gc()
tcga_BRCA_metadata <- tcga_BRCA_metadata |>
  as.data.frame()|>
  dplyr::rename(SOURCE = 4, ID = 5)|>
  dplyr::mutate(DIAGNOSIS = case_when(SOURCE == 'Solid Tissue Normal' ~ 'Control',
                                      SOURCE == 'Primary solid Tumor' ~ 'Cancer',
                                      .default = SOURCE))|>
  dplyr::filter(DIAGNOSIS == 'Control'|DIAGNOSIS=='Cancer')
gc()
exprDATA <- readRDS('TCGA_BRCA.RDS') |>
  dplyr::select(-c(HGNC,GeneName))|>
  drop_na(GeneID)|>
  distinct()|>
  dplyr::mutate(GeneID2 = sapply(GeneID, function(x){str_split_i(x,'\\.',1)}))|>
  dplyr::filter(GeneID2 %in% gtf2$gene_id)|>
  dplyr::select(-GeneID)|>
  dplyr::rename(GeneID = GeneID2)|>
  dplyr::select(GeneID,tcga_BRCA_metadata$specimenID)

dim(exprDATA) # [1] 32610   601
variance_cols <- colnames(exprDATA)[-1]
exprDATA <- exprDATA %>%
  mutate(row_variance = apply(.[, variance_cols], 1, var)) %>%
  group_by(GeneID) %>%
  slice_max(order_by = row_variance, n = 1, with_ties = F) %>%
  ungroup()|>
  dplyr::select(-row_variance)|>
  tidyr::drop_na()|>
  remove_rownames()|>column_to_rownames(var = 'GeneID')|>
  dplyr::select(tcga_BRCA_metadata$specimenID)|>
  TMM(tcga_BRCA_metadata$DIAGNOSIS)
dim(exprDATA) # 14451  1233
KEY <- tcga_BRCA_metadata$DIAGNOSIS[match(colnames(exprDATA),tcga_BRCA_metadata$specimenID)] 
# Cancer Control 
# 1095     112 
print(table(KEY))
# exprDATA <- exprDATA[-which(apply(exprDATA,1,var)==0),]
# exprDATA1 <- matrix(data = sapply(exprDATA,as.numeric),nrow = nrow(exprDATA),ncol = ncol(exprDATA))
print(pcaPlot(exprDATA,KEY,'BRCA'))
ggsave(paste0('./TMM/PCA Plots/BRCA_beforeCovAdj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
exprDATA_CovAdj <- covAdj_v2(tcga_BRCA_metadata,exprDATA)
max(exprDATA_CovAdj);min(exprDATA_CovAdj)
# exprDATA_CovAdj <- exprDATA_CovAdj + abs(exprDATA_CovAdj)+ 1
# max(exprDATA_CovAdj);min(exprDATA_CovAdj)
print(pcaPlot(exprDATA_CovAdj,KEY,'BRCA'))
ggsave(paste0('./TMM/PCA Plots/BRCA_afterCovAdj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
# DEG

degBRCA <- degAnalysisLIMMA(Group = tcga_BRCA_metadata$DIAGNOSIS, 
                            df = exprDATA_CovAdj|> as.data.frame()|>dplyr::select(tcga_BRCA_metadata$specimenID))

openxlsx::write.xlsx(degBRCA, './TMM/DEG Data/BRCA_DEG.xlsx')
df <- list()
df[['Control']] <- exprDATA_CovAdj|> as.data.frame()|> dplyr::select(tcga_BRCA_metadata$specimenID[which(tcga_BRCA_metadata$DIAGNOSIS == 'Control')])|> rownames_to_column(var = 'GeneID')
df[['Cancer']] <- exprDATA_CovAdj|> as.data.frame()|> dplyr::select(tcga_BRCA_metadata$specimenID[which(tcga_BRCA_metadata$DIAGNOSIS == 'Cancer')])|> rownames_to_column(var = 'GeneID')
openxlsx::write.xlsx(df, './TMM/Expression Data/BRCA.xlsx')
gc()
# Colorectal Cancer (COAD)
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)
COAD.data <- GDCprepare(query,
                       save = TRUE,
                       save.filename = "TCGA_COAD.rda",
                       remove.files.prepared = TRUE
)

# Metadata
tcga_COAD_metadata <- data.frame(COAD.data$patient,COAD.data$sample,
                                 COAD.data$tumor_descriptor,COAD.data$definition,
                                 COAD.data$barcode, COAD.data$days_to_diagnosis,
                                 COAD.data$days_to_last_follow_up,COAD.data$age_at_diagnosis,
                                 COAD.data$primary_diagnosis,COAD.data$diagnosis_id)

saveRDS(tcga_COAD_metadata,'tcga_COAD_metadata.RDS')
data2 <- COAD.data@assays@data@listData$unstranded |>
  as.data.frame()|>
  dplyr::rename_with(~COAD.data@colData@rownames)|>
  dplyr::mutate(GeneID = COAD.data@rowRanges@elementMetadata$gene_id,
                GeneName = COAD.data@rowRanges@elementMetadata$gene_name,
                HGNC = COAD.data@rowRanges@elementMetadata$hgnc_id)
saveRDS(data2, 'TCGA_COAD.RDS')

tcga_COAD_metadata <- readRDS('tcga_COAD_metadata.RDS')|>
  as.data.frame()|>
  dplyr::rename(SOURCE = 4, ID = 5)|>
  dplyr::mutate(DIAGNOSIS = case_when(SOURCE == 'Solid Tissue Normal' ~ 'Control',
                                      SOURCE == 'Primary solid Tumor' ~ 'Cancer',
                                      .default = SOURCE))|>
  dplyr::filter(DIAGNOSIS == 'Control'|DIAGNOSIS=='Cancer')|>
  dplyr::rename(specimenID = 5,
                AGE = 8)|>
  drop_na(AGE)
gc()

exprDATA <- readRDS('TCGA_COAD.RDS') |>
  dplyr::select(-c(HGNC,GeneName))|>
  drop_na(GeneID)|>
  distinct()|>
  dplyr::mutate(GeneID2 = sapply(GeneID, function(x){str_split_i(x,'\\.',1)}))|>
  dplyr::filter(GeneID2 %in% gtf2$gene_id)|>
  dplyr::select(-GeneID)|>
  dplyr::rename(GeneID = GeneID2)|>
  dplyr::select(GeneID,tcga_COAD_metadata$specimenID)

dim(exprDATA) # [1] 32610   601
variance_cols <- colnames(exprDATA)[-1]
exprDATA <- exprDATA %>%
  mutate(row_variance = apply(.[, variance_cols], 1, var)) %>%
  group_by(GeneID) %>%
  slice_max(order_by = row_variance, n = 1, with_ties = F) %>%
  ungroup()|>
  dplyr::select(-row_variance)|>
  tidyr::drop_na()|>
  remove_rownames()|>column_to_rownames(var = 'GeneID')|>
  dplyr::select(tcga_COAD_metadata$specimenID)|>
  TMM(tcga_COAD_metadata$DIAGNOSIS)
dim(exprDATA) 
KEY <- tcga_COAD_metadata$DIAGNOSIS[match(colnames(exprDATA),tcga_COAD_metadata$specimenID)] 
# Cancer Control 
# 477      41
print(table(KEY))
# exprDATA <- exprDATA[-which(apply(exprDATA,1,var)==0),]
# exprDATA1 <- matrix(data = sapply(exprDATA,as.numeric),nrow = nrow(exprDATA),ncol = ncol(exprDATA))
print(pcaPlot(exprDATA,KEY,'COAD'))
ggsave(paste0('./TMM/PCA Plots/COAD_beforeCovAdj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
exprDATA_CovAdj <- covAdj_v2(tcga_COAD_metadata,exprDATA)
max(exprDATA_CovAdj);min(exprDATA_CovAdj)

print(pcaPlot(exprDATA_CovAdj,KEY,'COAD'))
ggsave(paste0('./TMM/PCA Plots/COAD_afterCov_Adj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
# DEG

degCOAD <- degAnalysisLIMMA(Group = tcga_COAD_metadata$DIAGNOSIS, 
                            df = exprDATA_CovAdj|>as.data.frame() |>dplyr::select(tcga_COAD_metadata$specimenID))
openxlsx::write.xlsx(degCOAD, './TMM/DEG Data/COAD_DEG.xlsx')
df <- list()
df[['Control']] <- exprDATA_CovAdj|>as.data.frame()|> dplyr::select(tcga_COAD_metadata$specimenID[which(tcga_COAD_metadata$DIAGNOSIS == 'Control')])|> rownames_to_column(var = 'GeneID')
df[['Cancer']] <- exprDATA_CovAdj|>as.data.frame()|> dplyr::select(tcga_COAD_metadata$specimenID[which(tcga_COAD_metadata$DIAGNOSIS == 'Cancer')])|> rownames_to_column(var = 'GeneID')
openxlsx::write.xlsx(df, './TMM/Expression Data/COAD.xlsx')
gc()


# Liver Adernocarcinoma (LIHC)
query <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)
GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)
LIHC.data <- GDCprepare(query,
                       save = TRUE,
                       save.filename = "TCGA_LIHC.rda",
                       remove.files.prepared = TRUE
)
# Metadata
tcga_LIHC.metadata <- data.frame(LIHC.data$patient,LIHC.data$sample,
                                 LIHC.data$tumor_descriptor,LIHC.data$definition,
                                 LIHC.data$barcode, LIHC.data$days_to_diagnosis,
                                 LIHC.data$days_to_last_follow_up,LIHC.data$age_at_diagnosis,
                                 LIHC.data$primary_diagnosis,LIHC.data$diagnosis_id)

saveRDS(tcga_LIHC.metadata,'tcga_LIHC.metadata.RDS')
data2 <- LIHC.data@assays@data@listData$unstranded |>
  as.data.frame()|>
  dplyr::rename_with(~LIHC.data@colData@rownames)|>
  dplyr::mutate(GeneID = LIHC.data@rowRanges@elementMetadata$gene_id,
                GeneName = LIHC.data@rowRanges@elementMetadata$gene_name,
                HGNC = LIHC.data@rowRanges@elementMetadata$hgnc_id)
saveRDS(data2, 'TCGA_LIHC.RDS')
tcga_LIHC.metadata <- readRDS('tcga_LIHC.metadata.RDS')  |>
  as.data.frame()|>
  dplyr::rename(SOURCE = 4, ID = 5)|>
  dplyr::mutate(DIAGNOSIS = case_when(SOURCE == 'Solid Tissue Normal' ~ 'Control',
                                      SOURCE == 'Primary solid Tumor' ~ 'Cancer',
                                      .default = SOURCE))|>
  dplyr::filter(DIAGNOSIS == 'Control'|DIAGNOSIS=='Cancer')|>
  dplyr::rename(specimenID = 5,
                AGE = 8)|>
  drop_na(AGE)
gc()

exprDATA <- readRDS('TCGA_LIHC.RDS') |>
  dplyr::select(-c(HGNC,GeneName))|>
  drop_na(GeneID)|>
  distinct()|>
  dplyr::mutate(GeneID2 = sapply(GeneID, function(x){str_split_i(x,'\\.',1)}))|>
  dplyr::filter(GeneID2 %in% gtf2$gene_id)|>
  dplyr::select(-GeneID)|>
  dplyr::rename(GeneID = GeneID2)|>
  dplyr::select(GeneID,tcga_LIHC.metadata$specimenID)

dim(exprDATA) # [1] 32610   601
variance_cols <- colnames(exprDATA)[-1]
exprDATA <- exprDATA %>%
  mutate(row_variance = apply(.[, variance_cols], 1, var)) %>%
  group_by(GeneID) %>%
  slice_max(order_by = row_variance, n = 1, with_ties = F) %>%
  ungroup()|>
  dplyr::select(-row_variance)|>
  tidyr::drop_na()|>
  remove_rownames()|>column_to_rownames(var = 'GeneID')|>
  dplyr::select(tcga_LIHC.metadata$specimenID)|>
  TMM(tcga_LIHC.metadata$DIAGNOSIS)
dim(exprDATA)  
KEY <- tcga_LIHC.metadata$DIAGNOSIS[match(colnames(exprDATA),tcga_LIHC.metadata$specimenID)] 
# Cancer Control 
# 367      48
print(table(KEY))
# exprDATA <- exprDATA[-which(apply(exprDATA,1,var)==0),]
# exprDATA1 <- matrix(data = sapply(exprDATA,as.numeric),nrow = nrow(exprDATA),ncol = ncol(exprDATA))
print(pcaPlot(exprDATA,KEY,'LIHC'))
ggsave(paste0('./TMM/PCA Plots/LIHC_beforeCovAdj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
exprDATA_CovAdj <- covAdj_v2(tcga_LIHC.metadata,exprDATA)
max(exprDATA_CovAdj);min(exprDATA_CovAdj)
# exprDATA_CovAdj <- exprDATA_CovAdj + abs(exprDATA_CovAdj)+ 1
# max(exprDATA_CovAdj);min(exprDATA_CovAdj)
print(pcaPlot(exprDATA_CovAdj,KEY,'LIHC'))
ggsave(paste0('./TMM/PCA Plots/LIHC_afterCovAdj.tiff'), 
       units="in", width=6.5, height=7.5, dpi=600, compression = 'lzw')
# DEG

degLIHC <- degAnalysisLIMMA(Group = tcga_LIHC.metadata$DIAGNOSIS, 
                            df = exprDATA_CovAdj|> as.data.frame()|>dplyr::select(tcga_LIHC.metadata$specimenID))
openxlsx::write.xlsx(degLIHC, './TMM/DEG Data/LIHC_DEG.xlsx')
df <- list()
df[['Control']] <- exprDATA_CovAdj|> as.data.frame()|> dplyr::select(tcga_LIHC.metadata$specimenID[which(tcga_LIHC.metadata$DIAGNOSIS == 'Control')])|> rownames_to_column(var = 'GeneID')
df[['Cancer']] <- exprDATA_CovAdj|> as.data.frame()|> dplyr::select(tcga_LIHC.metadata$specimenID[which(tcga_LIHC.metadata$DIAGNOSIS == 'Cancer')])|> rownames_to_column(var = 'GeneID')
openxlsx::write.xlsx(df, './TMM/Expression Data/LIHC.xlsx')
gc()

rm(list = ls())
gc()
