# Prepare metabolite data from the hmdb and perform Fisher's-based disease over-representation analysis
# BiocManager::install('OmnipathR')
library(OmnipathR)
library(data.table)
library(tidyverse)
setwd('..')

name_map_KEGGID <- readr::read_csv("./name_map_KEGGID.csv")
met_disease_map <- fread("./hmdb_disease_data.tsv")|>
  dplyr::mutate(KEGGID = name_map_KEGGID$Query[match(HMDB_ID,name_map_KEGGID$HMDB)])

# Sanity checks for some diseases in the study
ad.METS <- met_disease_map |> dplyr::filter(Disease_Name=="Alzheimer's disease")
pd.METS <- met_disease_map |> dplyr::filter(Disease_Name=="Parkinson's disease")


met_disease_map.1 <- met_disease_map |>
  dplyr::rename(Metabolite=2,Disease=3)
# Export data for annotation using the "annotate_hmdb_metabolites.R" script
openxlsx::write.xlsx(met_disease_map.1,'../hmdb_disease_data_v2.xlsx')

# Load annotated metabolites
ref_met2dis <- openxlsx::read.xlsx('E:/Metabolic Modeling Tools/HMDB Data/met2dis_Ref.xlsx')

# Load iMAT derived personalized models summarized by Fisher's Exact Test
sheetnames <- openxlsx::getSheetNames('../iMAT Mapping/2025-07-09_AMS_iMATFisher.xlsx')
metList <- lapply(sheetnames, function(x){
  metDrugTarget2(openxlsx::read.xlsx('../iMAT Mapping/2025-07-09_AMS_iMATFisher.xlsx',
                                     sheet = paste0(x))|>
                   dplyr::mutate(RxnActivity = ifelse(OR==65535 & (ActiveAD>ActiveCon), (ActiveAD+1)/(ActiveCon+1), OR),
                                 # RxnActivity = ifelse(is.na(RxnActivity), 0, RxnActivity),
                                 FDR = p.adjust(PValue, method = 'fdr'))|> 
                   dplyr::filter(PValue < 0.05 & (OR > 1 | RxnActivity>1)),
                 HumanGEMMetabolite,
                 Human_GEM,
                 CurrencyMets,
                 metaboliteGEM,
                 metaboliteGEM.POOL,
                 met2RxnDATA,
                 topNPercent=.95)$FullMetaboliteList
})|>
  magrittr::set_names(paste0(sheetnames))
# Fisher's based disease over-representation analysis

overrepAnalysis <- lapply(metList, function(x){
  tarMets <- unique(x$KEGGID)
  tarMets <- tarMets[!is.na(tarMets)]
  
  # Remove duplicate associations
  disAssoc <- unique(ref_met2dis)
  
  # Define the "universe" of all unique metabolites in your background set
  allMets <- unique(disAssoc$KEGGID); allMets <- allMets[!is.na(allMets)]
  N <- length(allMets)

  tarMets <- intersect(tarMets, allMets)
  k <- length(tarMets)
  
  # all unique diseases to test
  allDis <- unique(disAssoc$Disease)
  
  results <- list()
  # Loop through each disease to perform the test
  for (disName in allDis) {
    
    # Get all metabolites associated with the current disease
    met2dis <- subset(disAssoc, Disease == disName)$KEGGID
    
    # --- Calculate the 4 numbers for the hypergeometric test ---
    
    # q: The number of metabolites in your list that ARE associated with the disease.
    q <- length(intersect(tarMets, met2dis))
    
    # m: The total number of metabolites in the universe that ARE associated with the disease.
    m <- length(met2dis)
    
    # n: The total number of metabolites in the universe that ARE NOT associated with the disease.
    n <- N - m
    
    # k: The number of metabolites in your list of interest. (defined earlier)
    
    # Perform the hypergeometric test
    pValue <- phyper(q - 1, m, n, k, lower.tail = FALSE)
    results[[disName]] <- data.frame(
      Disease = disName,
      PValue = pValue,
      `number of Metabolites Altered in disease` = q,
      `Number of Altered Metabolites` = k,
      `number of Metabolites in Disease` = m,
      `number of All Metabolites` = N
    )
  }
  dfResults <- do.call(rbind, results)
  
  # Adjust p-values for multiple testing using the Benjamini-Hochberg method
  dfResults$FDR <- p.adjust(dfResults$p_value, method = "BH")
  dfResults <- dfResults[order(dfResults$p_adjusted), ]|>
    remove_rownames()
  return(dfResults)
})


overrepAnalysis_filtered <- lapply(overrepAnalysis, function(x){
  x |> dplyr::filter(p_adjusted < 0.05)
})

openxlsx::write.xlsx(overrepAnalysis,'../Disease Ontology/met2Dis_OverrepAnalysis.xlsx')
