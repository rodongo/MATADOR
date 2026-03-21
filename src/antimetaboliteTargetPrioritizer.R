################################################################################
# metDrugTarget.R
#
# Functions for metabolic network construction and drug target prioritisation
# from iMAT flux results. Metabolites are ranked by betweenness and closeness
# centrality in weighted and unweighted substrate-product graphs.
#
# Two analysis strategies are provided:
#   - metDrugTarget()  : module-aware (decomposes graph into strongly connected
#                        components; top-N% of metabolites per module selected)
#   - metDrugTarget2() : whole-graph strategy (all metabolites ranked globally)
#
# Dependencies: tidyverse, igraph, BioNet, openxlsx, readr, stringr, stringi
# Regan Odongo, 2024-06-20
################################################################################

library(tidyverse)
library(igraph)
library(BioNet)
library(openxlsx)
library(readr)
library(stringr)
library(stringi)


# ------------------------------------------------------------------------------
#' Compute betweenness and closeness centrality per strongly-connected module
#'
#' Decomposes a graph into strongly connected components (≥ 4 vertices) and
#' returns centrality scores for each component separately.
#'
#' @param graph   An igraph object. Edge weights (if present) are used when
#'                \code{directed = TRUE}.
#' @param directed Logical. If \code{TRUE}, uses edge weights and directional
#'                 closeness (\code{mode = "out"}); otherwise uses unweighted
#'                 undirected centrality.
#'
#' @return A named list of data frames (one per module), each with columns:
#'   \describe{
#'     \item{Metabolite}{Vertex name}
#'     \item{Betweenness}{Betweenness centrality score}
#'     \item{Closeness}{Closeness centrality score}
#'   }
#'   List elements are named \code{Module_1}, \code{Module_2}, …
# ------------------------------------------------------------------------------
subnetTopographyModular <- function(graph, directed = TRUE) {
  subnetworks <- decompose(graph, mode = "strong", max.comps = NA, min.vertices = 4)
  netAnalysis <- lapply(seq_along(subnetworks), function(i) {
    net <- subnetworks[[i]]
    if (directed) {
      Betweenness <- igraph::betweenness(net, v = V(net), weights = E(net)$weight)
      Closeness   <- igraph::closeness(net, vids = V(net), weights = E(net)$weight, mode = "out")
    } else {
      Betweenness <- igraph::betweenness(net, v = V(net), weights = NULL)
      Closeness   <- igraph::closeness(net, vids = V(net), weights = NULL, mode = "all")
    }
    data.frame(
      Metabolite  = V(net)$name,
      Betweenness = Betweenness,
      Closeness   = Closeness
    )
  })
  names(netAnalysis) <- paste0("Module_", seq_along(subnetworks))
  netAnalysis
}


# ------------------------------------------------------------------------------
#' Compute betweenness and closeness centrality on the full graph
#'
#' Unlike \code{subnetTopographyModular}, this function operates on the entire
#' graph rather than decomposed modules, giving a global centrality ranking.
#'
#' @param graph   An igraph object.
#' @param directed Logical. See \code{subnetTopographyModular}.
#'
#' @return A single data frame with columns \code{Metabolite},
#'   \code{Betweenness}, and \code{Closeness}.
# ------------------------------------------------------------------------------
subnetTopographyGlobal <- function(graph, directed = TRUE) {
  if (directed) {
    Betweenness <- igraph::betweenness(graph, v = V(graph), weights = E(graph)$weight, directed = TRUE)
    Closeness   <- igraph::closeness(graph, vids = V(graph), weights = E(graph)$weight, mode = "out")
  } else {
    Betweenness <- igraph::betweenness(graph, v = V(graph), weights = NULL, directed = TRUE)
    Closeness   <- igraph::closeness(graph, vids = V(graph), weights = NULL, mode = "all")
  }
  data.frame(
    Metabolite  = V(graph)$name,
    Betweenness = Betweenness,
    Closeness   = Closeness
  )
}


# ------------------------------------------------------------------------------
#' Extract reaction context for a set of metabolites
#'
#' For each metabolite, retrieves associated reactions (all, exchange-only,
#' producing, consuming) from a substrate-product edge table.
#'
#' @param mets    Character vector of metabolite names to query.
#' @param cenData Data frame with columns \code{SubstrateMet}, \code{ProductMet},
#'                \code{AssoRxn}, \code{ExRxn}, and \code{AssociatedConsuming}.
#'
#' @return A data frame (one row per unique \code{AssociatedConsuming} entry)
#'   with columns: \code{Metabolite}, \code{AssociatedRxn}, \code{ExchangeRxn},
#'   \code{AssociatedConsuming}, \code{ProducingRxn}, \code{ConsumingRxn}.
# ------------------------------------------------------------------------------
extractMetCent <- function(mets, cenData) {
  betMat <- data.frame(
    Metabolite         = mets,
    AssociatedRxn      = NA_character_,
    ExchangeRxn        = NA_character_,
    AssociatedConsuming = NA_character_,
    ProducingRxn       = NA_character_,
    ConsumingRxn       = NA_character_,
    stringsAsFactors   = FALSE
  )

  for (a in seq_along(mets)) {
    met <- mets[a]
    is_sub <- cenData$SubstrateMet == met
    is_pro <- cenData$ProductMet   == met

    Rxn_all <- cenData$AssoRxn[is_sub | is_pro]
    Rxn_ex  <- cenData$ExRxn[is_sub | is_pro]
    Rxn_ex  <- Rxn_ex[Rxn_ex != "Not Exchange"]

    betMat$AssociatedRxn[a] <- paste(unique(Rxn_all), collapse = ", ")
    betMat$ExchangeRxn[a]   <- paste(unique(Rxn_ex),  collapse = ", ")

    consuming_idx <- is_sub & !is_pro
    producing_idx <- !is_sub & is_pro

    betMat$ProducingRxn[a]       <- paste(unique(cenData$AssoRxn[producing_idx]), collapse = ", ")
    betMat$ConsumingRxn[a]       <- paste(unique(cenData$AssoRxn[consuming_idx]), collapse = ", ")
    betMat$AssociatedConsuming[a] <- paste(unique(cenData$AssociatedConsuming[consuming_idx]), collapse = ", ")
  }

  betMat[!duplicated(betMat$AssociatedConsuming), ]
}


# ------------------------------------------------------------------------------
#' Build a substrate-product edge table from iMAT results
#'
#' Internal helper: filters significant reactions, parses metabolic equations,
#' expands substrate-product pairs, and removes currency/pool metabolites.
#'
#' @param df           Data frame of iMAT results with columns \code{Reaction}
#'                     and \code{PValue}.
#' @param rxnMetEqn    Data frame with columns \code{Reactions} and \code{Equations}.
#' @param rxnPATH      Data frame with columns \code{ID} and \code{SUBSYSTEM}.
#' @param curMETS      Character vector of currency metabolite names to exclude.
#' @param poolMETS     Character vector of pool metabolite names to exclude.
#'
#' @return A data frame with columns \code{SubstrateMet}, \code{ProductMet},
#'   \code{RxnWeight}, \code{AssoRxn}, \code{ExRxn}.
# ------------------------------------------------------------------------------
.buildEdgeTable <- function(df, rxnMetEqn, rxnPATH, curMETS, poolMETS) {
  df <- df |>
    dplyr::filter(PValue < 0.05) |>
    dplyr::mutate(
      Equation  = rxnMetEqn$Equations[match(Reaction, rxnMetEqn$Reactions)],
      Subsystem = rxnPATH$SUBSYSTEM[match(Reaction, rxnPATH$ID)]
    ) |>
    dplyr::select(Reaction, Subsystem, Equation, PValue)

  met2RxnMat <- data.frame(
    SubstrateMet = character(),
    ProductMet   = character(),
    RxnWeight    = numeric(),
    AssoRxn      = character(),
    ExRxn        = character(),
    stringsAsFactors = FALSE
  )

  for (l in seq_len(nrow(df))) {
    eqn      <- df$Equation[l]
    weight   <- -log10(as.numeric(df$PValue[l]))
    is_exchange <- "Exchange/demand reactions" %in% df$Subsystem[l]

    # Split on => or ->
    parts <- str_split(eqn, "=>", n = 2)[[1]]
    if (length(parts) < 2) parts <- str_split(eqn, "->", n = 2)[[1]]

    subs_str <- gsub("<", "", parts[1])
    prod_str <- gsub(">", "", parts[2])

    # Exchange reactions: product side may be empty
    if (trimws(prod_str) == "") prod_str <- subs_str

    parse_side <- function(s) {
      x <- str_split(s, " \\+ ")[[1]]
      x <- str_trim(x, side = "both")
      # Remove stoichiometric coefficients (integer or decimal at start)
      x <- sub("^[0-9]*\\.?[0-9]+[[:space:]]+", "", x)
      x
    }

    subs <- parse_side(subs_str)
    prod <- parse_side(prod_str)

    pairs <- expand.grid(SubstrateMet = subs, ProductMet = prod, stringsAsFactors = FALSE)
    pairs$RxnWeight <- weight
    pairs$AssoRxn   <- df$Reaction[l]
    pairs$ExRxn     <- ifelse(is_exchange, df$Reaction[l], "Not Exchange")

    met2RxnMat <- rbind(met2RxnMat, pairs)
  }

  met2RxnMat |>
    as.data.frame() |>
    tidyr::drop_na(SubstrateMet) |>
    dplyr::filter(
      !SubstrateMet %in% curMETS,
      !ProductMet   %in% curMETS,
      !SubstrateMet %in% poolMETS,
      !ProductMet   %in% poolMETS,
      SubstrateMet  != "",
      ProductMet    != ""
    )
}


# ------------------------------------------------------------------------------
#' Build and clean the metabolite-metabolite igraph objects
#'
#' Internal helper: constructs weighted and unweighted undirected graphs from
#' an edge table, removes isolated vertices and self-loops.
#'
#' @param met2RxnMat Edge table from \code{.buildEdgeTable}.
#'
#' @return A list with elements \code{Weighted} and \code{Unweighted},
#'   each an igraph object.
# ------------------------------------------------------------------------------
.buildMetGraphs <- function(met2RxnMat) {
  edge_df <- data.frame(
    Substrate = met2RxnMat$SubstrateMet,
    Product   = met2RxnMat$ProductMet,
    Weight    = met2RxnMat$RxnWeight
  )

  g_raw <- graph_from_data_frame(edge_df, directed = FALSE)
  g_raw <- delete.vertices(g_raw, which(igraph::degree(g_raw, mode = "out") == 0))

  clean_graph <- function(g, add_weights = FALSE) {
    if (add_weights) E(g)$weight <- edge_df$Weight
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    delete.vertices(g, igraph::degree(g) == 0)
  }

  list(
    Weighted   = clean_graph(g_raw, add_weights = TRUE),
    Unweighted = clean_graph(g_raw, add_weights = FALSE)
  )
}


# ------------------------------------------------------------------------------
#' Annotate edge table with consuming reactions from a reference dataset
#'
#' @param met2RxnMat   Edge table (output of \code{.buildEdgeTable}).
#' @param metAssRxnDATA Data frame with columns \code{Substrate} and
#'                     \code{AssociatedReaction}.
#'
#' @return \code{met2RxnMat} with an additional \code{AssociatedConsuming} column.
# ------------------------------------------------------------------------------
.annotateConsuming <- function(met2RxnMat, metAssRxnDATA) {
  met2RxnMat$AssociatedConsuming <- vapply(
    met2RxnMat$SubstrateMet,
    function(sub) {
      rxns <- metAssRxnDATA$AssociatedReaction[metAssRxnDATA$Substrate == sub]
      paste(unique(rxns), collapse = ", ")
    },
    character(1)
  )
  met2RxnMat
}


# ------------------------------------------------------------------------------
#' Post-filter centrality result tables
# ------------------------------------------------------------------------------
.filterCentralityTable <- function(df) {
  df |>
    as.data.frame() |>
    dplyr::filter(!grepl("pool", Metabolite, ignore.case = TRUE)) |>
    dplyr::filter(AssociatedConsuming != "") |>
    tidyr::drop_na(Metabolite) |>
    dplyr::distinct(Metabolite, .keep_all = TRUE)
}


# ------------------------------------------------------------------------------
#' Module-aware metabolic drug target prioritisation
#'
#' Builds a metabolite-metabolite interaction network from iMAT flux results,
#' decomposes it into strongly connected modules, and ranks metabolites by
#' betweenness and closeness centrality within each module. The top
#' \code{topNPercent} of metabolites per module are returned as candidate
#' drug targets.
#'
#' @param df            Data frame of iMAT results with columns \code{Reaction}
#'                      and \code{PValue}.
#' @param rxnMetEqn     Data frame with columns \code{Reactions} and
#'                      \code{Equations}.
#' @param rxnPATH       Data frame with columns \code{ID} and \code{SUBSYSTEM}.
#' @param curMETS       Character vector of currency metabolite names to exclude
#'                      (e.g. ATP, NADH).
#' @param metKEGG       Data frame with columns \code{METABOLITE} and \code{KEGGID}.
#' @param poolMETS      Character vector of pool metabolite names to exclude.
#' @param metAssRxnDATA Data frame with columns \code{Substrate} and
#'                      \code{AssociatedReaction}.
#' @param topNPercent   Numeric in (0, 1]. Fraction of top-ranked metabolites
#'                      to retain per module.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{FullMetaboliteList}{All metabolites in the network with KEGG IDs}
#'     \item{Network}{The full igraph network object}
#'     \item{Subnetwork}{List with \code{Weighted} and \code{Unweighted} module centrality tables}
#'     \item{met2RxnMat}{Annotated substrate-product edge table}
#'     \item{BetweennessTop}{Top metabolites ranked by betweenness centrality}
#'     \item{ClosenessTop}{Top metabolites ranked by closeness centrality}
#'   }
#'
#' @seealso \code{\link{metDrugTarget2}} for a whole-graph alternative.
#' @export
# ------------------------------------------------------------------------------
metDrugTarget <- function(df,
                          rxnMetEqn,
                          rxnPATH,
                          curMETS,
                          metKEGG,
                          poolMETS,
                          metAssRxnDATA,
                          topNPercent) {
  netRes <- list()

  met2RxnMat <- .buildEdgeTable(df, rxnMetEqn, rxnPATH, curMETS, poolMETS)
  graphs      <- .buildMetGraphs(met2RxnMat)

  # Full metabolite list with KEGG annotation
  all_mets <- unique(c(met2RxnMat$SubstrateMet, met2RxnMat$ProductMet))
  netRes[["FullMetaboliteList"]] <- data.frame(
    MetaboliteName = all_mets,
    KEGGID         = metKEGG$KEGGID[match(all_mets, metKEGG$METABOLITE)],
    stringsAsFactors = FALSE
  )
  netRes[["Network"]] <- graphs$Weighted

  # Module-level centrality
  topo_W <- subnetTopographyModular(graphs$Weighted,   directed = TRUE)
  topo_U <- subnetTopographyModular(graphs$Unweighted, directed = FALSE)
  netRes[["Subnetwork"]] <- list(Weighted = topo_W, Unweighted = topo_U)

  # Annotate consuming reactions
  met2RxnMat <- .annotateConsuming(met2RxnMat, metAssRxnDATA)
  netRes[["met2RxnMat"]] <- met2RxnMat

  # Helper: select top-N% from each module and extract reaction context
  topFromModules <- function(topo, centrality_col, net_label) {
    mets <- do.call(c, lapply(topo, function(mod) {
      mod <- dplyr::arrange(mod, dplyr::desc(.data[[centrality_col]]))
      mod$Metabolite[seq_len(round(topNPercent * nrow(mod)))]
    }))
    extractMetCent(mets = mets, cenData = met2RxnMat) |>
      dplyr::mutate(NetType = net_label)
  }

  netRes[["BetweennessTop"]] <- .filterCentralityTable(rbind(
    topFromModules(topo_W, "Betweenness", "Weighted"),
    topFromModules(topo_U, "Betweenness", "Unweighted")
  ))

  netRes[["ClosenessTop"]] <- .filterCentralityTable(rbind(
    topFromModules(topo_W, "Closeness", "Weighted"),
    topFromModules(topo_U, "Closeness", "Unweighted")
  ))

  netRes
}


# ------------------------------------------------------------------------------
#' Whole-graph metabolic drug target prioritisation
#'
#' Identical workflow to \code{metDrugTarget} but centrality is computed on the
#' full network rather than per module, and all metabolites are returned ranked
#' (no \code{topNPercent} truncation). Use this strategy when the network is
#' densely connected or module decomposition yields very few components.
#'
#' @inheritParams metDrugTarget
#'
#' @return A named list with the same structure as \code{metDrugTarget}, except
#'   that \code{Subnetwork} contains single data frames (not lists of modules).
#'   \code{BetweennessTop} and \code{ClosenessTop} contain all metabolites,
#'   globally ranked.
#'
#' @seealso \code{\link{metDrugTarget}} for the module-aware alternative.
#' @export
# ------------------------------------------------------------------------------
metDrugTarget2 <- function(df,
                           rxnMetEqn,
                           rxnPATH,
                           curMETS,
                           metKEGG,
                           poolMETS,
                           metAssRxnDATA,
                           topNPercent = NULL) {
  netRes <- list()

  met2RxnMat <- .buildEdgeTable(df, rxnMetEqn, rxnPATH, curMETS, poolMETS)
  graphs      <- .buildMetGraphs(met2RxnMat)

  # Full metabolite list with KEGG annotation
  all_mets <- unique(c(met2RxnMat$SubstrateMet, met2RxnMat$ProductMet))
  full_met_df <- data.frame(
    MetaboliteName = all_mets,
    KEGGID         = metKEGG$KEGGID[match(all_mets, metKEGG$METABOLITE)],
    stringsAsFactors = FALSE
  )
  netRes[["FullMetaboliteList"]] <- full_met_df[!duplicated(full_met_df), ]
  netRes[["Network"]] <- graphs$Weighted

  # Whole-graph centrality
  topo_W <- subnetTopographyGlobal(graphs$Weighted,   directed = TRUE)
  topo_U <- subnetTopographyGlobal(graphs$Unweighted, directed = FALSE)
  netRes[["Subnetwork"]] <- list(Weighted = topo_W, Unweighted = topo_U)

  # Annotate consuming reactions
  met2RxnMat <- .annotateConsuming(met2RxnMat, metAssRxnDATA)
  netRes[["met2RxnMat"]] <- met2RxnMat[!duplicated(met2RxnMat), ]

  # Extract centrality-ordered reaction context for all metabolites
  rankAndExtract <- function(topo, centrality_col, net_label) {
    ranked <- dplyr::arrange(topo, dplyr::desc(.data[[centrality_col]]))
    extractMetCent(mets = ranked$Metabolite, cenData = met2RxnMat) |>
      dplyr::mutate(NetType = net_label)
  }

  netRes[["BetweennessTop"]] <- .filterCentralityTable(rbind(
    rankAndExtract(topo_W, "Betweenness", "Weighted"),
    rankAndExtract(topo_U, "Betweenness", "Unweighted")
  ))

  netRes[["ClosenessTop"]] <- .filterCentralityTable(rbind(
    rankAndExtract(topo_W, "Closeness", "Weighted"),
    rankAndExtract(topo_U, "Closeness", "Unweighted")
  ))

  netRes
}