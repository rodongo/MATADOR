# =============================================================================
# Antimetabolite Target Prioritisation from iMAT Reaction Enrichment Results
# =============================================================================
#
# DESCRIPTION:
#   Integrates Human-GEM metabolic model annotations with iMAT-derived
#   reaction enrichment statistics (odds ratios and p-values) to identify
#   and prioritise candidate antimetabolite drug targets. For each comparison
#   group in the input enrichment file, the script:
#     1. Filters reactions reaching nominal significance (p < 0.05).
#     2. Applies FDR correction across all tested reactions.
#     3. Identifies reactions enriched in the disease state (OR > 1).
#     4. Maps significant reactions to their associated metabolites via
#        metDrugTarget2, excluding currency metabolites (e.g., ATP, NADH).
#     5. Exports per-comparison metabolite target lists to Excel.
#
# INPUT FILES:
#   antimetaboliteTargetPrioritizer.R              – helper functions including
#                                                    metDrugTarget2()
#   Human-GEM.xlsx                                 – Human-GEM reaction table
#                                                    (default sheet and METS sheet)
#   genes.tsv                                      – Human-GEM gene annotations
#   metabolites.tsv                                – Human-GEM metabolite table
#                                                    with KEGG IDs
#   HumanGEMMetabolite.xlsx                        – curated metabolite metadata
#   CurrencyMets.tsv                               – list of currency metabolite
#                                                    names to exclude
#   met2Rxn_v2.RDS                                 – metabolite-to-reaction
#                                                    mapping object
#   ./iMAT Mapping/2024-11-14_AMS_iMATFisher_edgeR.xlsx
#                                                  – multi-sheet workbook of
#                                                    iMAT reaction enrichment
#                                                    results; one sheet per
#                                                    comparison group
#
# OUTPUT FILES:
#   ./Metabolite Targets/<ComparisonName>.xlsx      – one file per comparison
#                                                    group containing the
#                                                    deduplicated metabolite
#                                                    target list (element [[6]]
#                                                    of metDrugTarget2 output)
#
# KEY PARAMETERS:
#   p-value threshold   : 0.05  (nominal, pre-FDR filter)
#   FDR method          : Benjamini-Hochberg ('fdr')
#   OR enrichment filter: OR > 1 (disease-enriched reactions only)
#   Connectivity cutoff : 0.95  (passed to metDrugTarget2; excludes the top
#                                5% most-connected metabolites as hubs)
#
# NOTES:
#   - OR == 65535 is a sentinel value (overflow/missing); replaced with 0.
#   - All file paths are hard-coded; update before deployment on a new system.
#   - metDrugTarget2 returns a list; element [[6]] contains the final
#     metabolite target table. Inspect other elements for intermediate outputs.
#
# DEPENDENCIES:
#   openxlsx, readr, dplyr, magrittr, tidyverse
#   metDrugTarget2 (sourced from antimetaboliteTargetPrioritizer.R)
#
# AUTHORS / VERSION:
#   MATADOR Toolbox pipeline  |  Antimetabolite target prioritisation module
# =============================================================================

# -----------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP AND HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Source helper script containing metDrugTarget2() and related utilities

source('antimetaboliteTargetPrioritizer.R')

library(tidyverse)
library(openxlsx)
library(readr)
library(magrittr)

# -----------------------------------------------------------------------------
# 1. LOAD HUMAN-GEM REFERENCE DATA
# -----------------------------------------------------------------------------

# Full Human-GEM reaction table (used for reaction ID filtering and annotation)
Human_GEM <- openxlsx::read.xlsx(
  "Human-GEM.xlsx"
)

# Gene-to-reaction mapping table from Human-GEM
metabolicGenes <- readr::read_delim(
  "genes.tsv",
  delim        = "\t",
  escape_double = FALSE,
  trim_ws       = TRUE
)

# Curated metabolite metadata (used by metDrugTarget2 for name resolution)
HumanGEMMetabolite <- openxlsx::read.xlsx(
  "HumanGEMMetabolite.xlsx"
)

# Currency metabolites to exclude (e.g., ATP, NADH, H2O).
# These are highly connected hub metabolites that would otherwise dominate
# the target list without reflecting meaningful drug target biology.
CurrencyMets <- readr::read_csv(
  "CurrencyMets.tsv",
  col_names = FALSE
)$X1

# -----------------------------------------------------------------------------
# 2. BUILD METABOLITE-TO-REACTION MAPPING TABLES
# -----------------------------------------------------------------------------

# Load metabolite TSV and retain only the metabolite ID and KEGG ID columns.
# These are used to annotate metabolite targets with external database IDs.
metabolites <- readr::read_delim(
  "metabolites.tsv",
  delim         = "\t",
  escape_double = FALSE,
  trim_ws       = TRUE
) |>
  dplyr::rename(METABOLITE.ID = 1, KEGGID = 4) |>
  dplyr::select(METABOLITE.ID, KEGGID)

# Build a metabolite name → KEGG ID lookup by joining the METS sheet of
# Human-GEM with the metabolites TSV on the shared metabolite ID field.
metaboliteGEM <- openxlsx::read.xlsx(
  "Human-GEM.xlsx",
  sheet = "METS"
) |>
  dplyr::rename(METABOLITE = 3, METABOLITE.ID = 9) |>
  dplyr::select(METABOLITE, METABOLITE.ID) |>
  dplyr::inner_join(metabolites, by = 'METABOLITE.ID') |>
  dplyr::select(-METABOLITE.ID)   # ID no longer needed after join

# Identify 'pool' metabolites (aggregated species, not individual metabolites).
# These are excluded downstream to prevent spurious target assignments.
metaboliteGEM.POOL <- metaboliteGEM$METABOLITE[
  grepl('pool', metaboliteGEM$METABOLITE, ignore.case = TRUE)
]

# Clean up the intermediate metabolites table; no longer needed
rm(metabolites)

# Pre-computed metabolite-to-reaction mapping object (serialised RDS).
# Generated upstream; contains reaction membership per metabolite.
met2RxnDATA <- readRDS(
  'met2Rxn_v2.RDS'
)

# -----------------------------------------------------------------------------
# 3. LOAD AND FILTER iMAT REACTION ENRICHMENT RESULTS
# -----------------------------------------------------------------------------

filesR <- "./iMAT Mapping/2026-03-09_AMS_iMATFisher.xlsx"

# Restrict Human-GEM to reactions present in the enrichment file.
# This reduces subsequent join sizes and confirms data consistency.
Human_GEM2 <- Human_GEM |>
  dplyr::filter(ID %in% openxlsx::read.xlsx(filesR)$Reaction)

# Read all comparison sheets from the enrichment workbook.
# Each sheet contains one comparison group's Fisher/edgeR enrichment results.
sheetNames <- openxlsx::getSheetNames(filesR)

compDATA <- lapply(sheetNames, function(x) {
  openxlsx::read.xlsx(filesR, sheet = x) |>
    dplyr::rename(PValue = 2) |>
    dplyr::mutate(
      # Replace sentinel overflow value (65535) and NA odds ratios with 0
      OR  = ifelse(OR == 65535, 0, OR),
      OR  = ifelse(is.na(OR),   0, OR),
      # Compute FDR-corrected p-values within each comparison group
      FDR = p.adjust(PValue, method = 'fdr')
    ) |>
    dplyr::filter(PValue < 0.05)   # retain nominally significant reactions
}) |>
  magrittr::set_names(sheetNames)

# -----------------------------------------------------------------------------
# 4. FILTER FOR DISEASE-ENRICHED REACTIONS (OR > 1)
# -----------------------------------------------------------------------------
# Reactions with OR > 1 are more active in the disease state and represent
# potential antimetabolite targets (blocking their associated metabolites
# may shift the metabolic phenotype back towards the healthy state).

compDATA_OR <- lapply(compDATA, function(x) {
  x |> dplyr::filter(OR > 1)
})

# Report the number of disease-enriched reactions per comparison group
message("Disease-enriched reactions (OR > 1) per comparison group:")
lapply(compDATA_OR, nrow) |> print()

# -----------------------------------------------------------------------------
# 5. MAP REACTIONS TO METABOLITE DRUG TARGETS
# -----------------------------------------------------------------------------
# metDrugTarget2 maps enriched reactions to their associated metabolites,
# filters out currency metabolites and pool species, and ranks candidates
# by connectivity. The 0.95 connectivity cutoff removes the top 5% most
# connected metabolites to avoid promiscuous hub assignments.

# -- Disease-enriched reactions only (OR > 1) --------------------------------
compDATA_Mets_OR <- lapply(compDATA_OR, function(x) {
  metDrugTarget2(
    x,
    HumanGEMMetabolite,
    Human_GEM,
    CurrencyMets,
    metaboliteGEM,
    metaboliteGEM.POOL,
    met2RxnDATA,
    .95
  )
})

# -- All significant reactions (p < 0.05, any OR) ----------------------------
# Retained for exploratory analysis; full target space before OR filtering.
compDATA_ALL_Mets <- lapply(compDATA, function(x) {
  metDrugTarget2(
    x,
    HumanGEMMetabolite,
    Human_GEM,
    CurrencyMets,
    metaboliteGEM,
    metaboliteGEM.POOL,
    met2RxnDATA,
    .95
  )
})

# -----------------------------------------------------------------------------
# 6. EXPORT METABOLITE TARGET LISTS
# -----------------------------------------------------------------------------
# Write the deduplicated metabolite target table (element [[6]] of the
# metDrugTarget2 output list) to a separate Excel file per comparison group.

lapply(seq_along(compDATA_Mets_OR), function(x) {
  ifelse(dir.exists('Metabolite Targets'),print('Created already'),dir.create('Metabolite Targets'))
  openxlsx::write.xlsx(
    compDATA_Mets_OR[[x]][[6]] |> dplyr::distinct(),
    file = paste0('Metabolite Targets/', names(compDATA_Mets_OR)[x], '.xlsx')
  )
})
message("Target prioritisation complete. Results written to ./Metabolite Targets/")