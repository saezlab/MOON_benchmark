# CytoSig MOON benchmark
#
# Manuscript section 2.2 / Figure 2.
#
# Default behavior intentionally uses the original MOON cache that was used for
# the submitted manuscript figures. The 2024-03-19 rerun is kept as a transparent
# network-correction comparison, but it is not mixed into the frozen paper
# figures unless explicitly selected.
#
# Fast paper reproduction:
#   Rscript scripts/citosig_moon.R
#
# Useful environment switches:
#   CITOSIG_SNAPSHOT=paper_original              # default
#   CITOSIG_SNAPSHOT=corrected_20240319          # use corrected rerun for plots
#   CITOSIG_WRITE_FIGURES=false                  # build objects without writing PDFs
#   CITOSIG_COMPARE_CORRECTED=true               # summarize original vs corrected
#   CITOSIG_EXPORT_PANEL_D_NETWORKS=true         # export SIF/ATT files for examples
#   CITOSIG_RECOMPUTE_TF=true                    # expensive, off by default
#   CITOSIG_RECOMPUTE_MOON=true                  # very expensive, off by default

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(reshape2)
  library(ggplot2)
  library(ggrepel)
  library(ggExtra)
  library(cosmosR)
  library(decoupleR)
})

source("scripts/support_functions.R")


# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)

  if (is.na(value) || value == "") {
    return(default)
  }

  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

load_single_rdata_object <- function(path) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)

  if (length(object_names) != 1) {
    stop("Expected exactly one object in ", path, ", found: ",
         paste(object_names, collapse = ", "))
  }

  env[[object_names]]
}

SNAPSHOTS <- list(
  paper_original = list(
    label = "paper_original",
    description = "Frozen manuscript snapshot using the original COSMOS PKN.",
    moon_cache = "results/citosig_moon_activties.RData",
    ligand_table = "results/moon_res_df_long_target_only_GPTcleanedup_withMeanScore.csv",
    ligand_summary_col = "median_score_by_target",
    meta_network = "support/neta_network.RData",
    scoring_network_characteristics = "results/moon_scoring_network_characteristics.RData"
  ),
  corrected_20240319 = list(
    label = "corrected_20240319",
    description = "Transparency rerun using the 2024-03-19 COSMOS/Omnipath PKN.",
    moon_cache = "results/citosig_moon_activties_20240319.RData",
    ligand_table = "results/moon_res_df_long_target_only_20240319.csv",
    ligand_summary_col = "mean_ligand_score",
    meta_network = "support/meta_network_20240319.RData",
    scoring_network_characteristics = NA_character_
  )
)

SNAPSHOT_NAME <- Sys.getenv("CITOSIG_SNAPSHOT", unset = "paper_original")
if (!SNAPSHOT_NAME %in% names(SNAPSHOTS)) {
  stop("Unknown CITOSIG_SNAPSHOT='", SNAPSHOT_NAME, "'. Available snapshots: ",
       paste(names(SNAPSHOTS), collapse = ", "))
}

SNAPSHOT <- SNAPSHOTS[[SNAPSHOT_NAME]]
WRITE_FIGURES <- env_flag("CITOSIG_WRITE_FIGURES", default = TRUE)
COMPARE_CORRECTED <- env_flag("CITOSIG_COMPARE_CORRECTED", default = FALSE)
EXPORT_PANEL_D_NETWORKS <- env_flag("CITOSIG_EXPORT_PANEL_D_NETWORKS", default = FALSE)
RECOMPUTE_TF <- env_flag("CITOSIG_RECOMPUTE_TF", default = FALSE)
RECOMPUTE_MOON <- env_flag("CITOSIG_RECOMPUTE_MOON", default = FALSE)
RECOMPUTE_SCORING_NETWORKS <- env_flag("CITOSIG_RECOMPUTE_SCORING_NETWORKS", default = FALSE)

FIGURE_DIR <- if (SNAPSHOT_NAME == "paper_original") {
  "results/figures"
} else {
  file.path("results/figures", SNAPSHOT$label)
}
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

message("Using CytoSig MOON snapshot: ", SNAPSHOT$label)
message(SNAPSHOT$description)
message("Full MOON recomputation enabled: ", RECOMPUTE_MOON)


# ------------------------------------------------------------------------------
# Manuscript figure map
# ------------------------------------------------------------------------------
#
# Figure 2A:
#   MOON score of each applied ligand in its matching CytoSig experiment.
#   Output: results/figures/moon_scores_ranked.pdf
#
# Figure 2B:
#   Mean within-experiment quantile and mean raw score of each ligand in
#   experiments where it was applied versus experiments where it was not applied.
#   Outputs:
#     results/figures/ligand_quantiles_long.pdf
#     results/figures/ligand_score_long.pdf
#
# Figure 2C:
#   Number of edges, nodes, and TFs in the MOON scoring networks.
#   Outputs:
#     results/figures/moon_scoring_network_characteristics.pdf
#     results/figures/moon_scoring_network_characteristics_scatter.pdf
#
# Figure 2D:
#   Network schematics assembled outside this script from exported SIF/ATT files
#   and the Cytoscape session. Set CITOSIG_EXPORT_PANEL_D_NETWORKS=true to export
#   example scoring networks without rerunning MOON.


# ------------------------------------------------------------------------------
# Shared data loaders
# ------------------------------------------------------------------------------

read_cytosig_metadata <- function(path = "data/cytosig/zscore_meta_clean_renamed.csv") {
  metadata <- as.data.frame(readr::read_csv(path, show_col_types = FALSE))
  metadata$experiment_target <- sub("[@&].*", "", metadata$id)
  metadata
}

read_zscore_matrix <- function(path = "data/cytosig/zscore_final_clean_filtered.csv") {
  zscores <- as.data.frame(readr::read_csv(path, show_col_types = FALSE))
  row.names(zscores) <- zscores$gene
  zscores[, setdiff(names(zscores), "gene"), drop = FALSE]
}

load_collectri <- function(path = "support/collectri.RData") {
  load_single_rdata_object(path)
}

load_meta_network <- function(path = SNAPSHOT$meta_network) {
  meta_network <- load_single_rdata_object(path)
  names(meta_network) <- c("source", "interaction", "target")
  meta_network
}

load_moon_cache <- function(path = SNAPSHOT$moon_cache) {
  moon_res_list <- load_single_rdata_object(path)

  if (!is.list(moon_res_list)) {
    stop("Expected moon cache to contain a list: ", path)
  }

  moon_res_list
}

moon_experiment_ids <- function(moon_res_list) {
  vapply(moon_res_list, function(x) names(x)[2], character(1))
}

save_panel <- function(plot, filename, width, height) {
  if (WRITE_FIGURES) {
    ggplot2::ggsave(
      filename = file.path(FIGURE_DIR, filename),
      plot = plot,
      width = width,
      height = height,
      units = "in"
    )
  }

  invisible(plot)
}


# ------------------------------------------------------------------------------
# Figure 2A: applied-ligand MOON scores
# ------------------------------------------------------------------------------

read_applied_ligand_table <- function(snapshot = SNAPSHOT) {
  ligand_table <- as.data.frame(
    readr::read_csv(snapshot$ligand_table, show_col_types = FALSE)
  )

  if (!snapshot$ligand_summary_col %in% names(ligand_table)) {
    stop("Missing expected summary column '", snapshot$ligand_summary_col,
         "' in ", snapshot$ligand_table)
  }

  ligand_table$summary_score <- ligand_table[[snapshot$ligand_summary_col]]

  ligand_table <- ligand_table[!is.na(ligand_table$id), , drop = FALSE]
  ligand_table$id <- as.character(ligand_table$id)
  ligand_table$source <- as.character(ligand_table$source)
  ligand_table
}

lookup_level_from_cache <- function(ligand_table, moon_res_list) {
  if ("level" %in% names(ligand_table)) {
    return(ligand_table)
  }

  ids <- moon_experiment_ids(moon_res_list)
  cache_index <- setNames(seq_along(ids), ids)

  ligand_table$level <- mapply(
    FUN = function(experiment_id, ligand) {
      cache_position <- unname(cache_index[experiment_id])

      if (length(cache_position) == 0 || is.na(cache_position)) {
        return(NA_integer_)
      }

      moon_df <- moon_res_list[[cache_position]]
      rows <- which(moon_df$source_original == ligand & moon_df$level != 0)

      if (length(rows) == 0) {
        return(NA_integer_)
      }

      as.integer(moon_df$level[rows[1]])
    },
    experiment_id = ligand_table$id,
    ligand = ligand_table$source
  )

  ligand_table
}

prepare_panel_2a_data <- function(ligand_table, moon_res_list) {
  panel_data <- lookup_level_from_cache(ligand_table, moon_res_list)
  panel_data <- panel_data[!is.na(panel_data$score), , drop = FALSE]

  panel_data <- panel_data %>%
    group_by(.data$source) %>%
    mutate(
      mean_ligand_score = mean(.data$score, na.rm = TRUE),
      sd_ligand_score = sd(.data$score, na.rm = TRUE)
    ) %>%
    ungroup()

  panel_data$quality_colour <- ifelse(
    panel_data$summary_score > 1.7,
    "lightgreen",
    ifelse(panel_data$summary_score > 0.3, "orange", "red")
  )

  panel_data <- panel_data[order(panel_data$summary_score, decreasing = TRUE), ]
  panel_data$source <- factor(panel_data$source, levels = unique(panel_data$source))
  panel_data$level <- as.factor(panel_data$level)
  panel_data
}

plot_panel_2a <- function(panel_data) {
  ggplot(panel_data, aes(x = .data$source, y = .data$score, group = .data$source)) +
    geom_boxplot(
      aes(fill = .data$quality_colour),
      coef = 6,
      color = "black",
      alpha = 1,
      outlier.shape = NA
    ) +
    geom_jitter(aes(color = .data$level), width = 0.2, height = 0, alpha = 0.8) +
    geom_hline(yintercept = 0) +
    scale_fill_identity() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "MOON score", color = "Level")
}


# ------------------------------------------------------------------------------
# Figure 2B: true-treatment versus background quantile/score
# ------------------------------------------------------------------------------

build_score_matrix <- function(moon_res_list) {
  ids <- moon_experiment_ids(moon_res_list)
  sources <- unique(unlist(lapply(moon_res_list, function(moon_df) {
    moon_df$source_original[moon_df$level != 0]
  }), use.names = FALSE))

  score_matrix <- matrix(
    NA_real_,
    nrow = length(sources),
    ncol = length(ids),
    dimnames = list(sources, ids)
  )

  for (i in seq_along(moon_res_list)) {
    moon_df <- moon_res_list[[i]]
    keep <- moon_df$level != 0
    score_matrix[moon_df$source_original[keep], i] <- as.numeric(moon_df[[2]][keep])
  }

  score_matrix
}

rank_by_experiment <- function(score_matrix) {
  apply(score_matrix, 2, function(scores) {
    rank(scores, na.last = "keep") / sum(!is.na(scores))
  })
}

build_panel_2b_summary <- function(score_matrix, applied_ligand_table) {
  rank_matrix <- rank_by_experiment(score_matrix)
  ligands <- unique(as.character(applied_ligand_table$source))

  ligand_summary <- lapply(ligands, function(ligand) {
    true_experiments <- unique(applied_ligand_table$id[applied_ligand_table$source == ligand])
    true_experiments <- intersect(true_experiments, colnames(score_matrix))
    false_experiments <- setdiff(colnames(score_matrix), true_experiments)

    if (!ligand %in% rownames(score_matrix)) {
      return(data.frame(
        ligand_to_rank = ligand,
        mean_quantile_in_trueExp = NA_real_,
        mean_quantile_in_falseExp = NA_real_,
        mean_score_in_trueExp = NA_real_,
        mean_score_in_falseExp = NA_real_
      ))
    }

    data.frame(
      ligand_to_rank = ligand,
      mean_quantile_in_trueExp = mean(rank_matrix[ligand, true_experiments], na.rm = TRUE),
      mean_quantile_in_falseExp = mean(rank_matrix[ligand, false_experiments], na.rm = TRUE),
      mean_score_in_trueExp = mean(score_matrix[ligand, true_experiments], na.rm = TRUE),
      mean_score_in_falseExp = mean(score_matrix[ligand, false_experiments], na.rm = TRUE)
    )
  })

  bind_rows(ligand_summary)
}

plot_panel_2b_quantiles <- function(panel_2b_summary) {
  plot_data <- reshape2::melt(
    panel_2b_summary[order(panel_2b_summary$mean_score_in_trueExp, decreasing = FALSE),
                     c("ligand_to_rank",
                       "mean_quantile_in_trueExp",
                       "mean_quantile_in_falseExp")],
    id.vars = "ligand_to_rank"
  )
  plot_data$ligand_to_rank <- factor(plot_data$ligand_to_rank,
                                     levels = unique(plot_data$ligand_to_rank))

  ggplot(plot_data, aes(x = .data$value, y = .data$ligand_to_rank, color = .data$variable)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(x = "Mean score quantile", y = NULL, color = NULL)
}

plot_panel_2b_scores <- function(panel_2b_summary) {
  plot_data <- reshape2::melt(
    panel_2b_summary[order(panel_2b_summary$mean_score_in_trueExp, decreasing = FALSE),
                     c("ligand_to_rank",
                       "mean_score_in_trueExp",
                       "mean_score_in_falseExp")],
    id.vars = "ligand_to_rank"
  )
  plot_data$ligand_to_rank <- factor(plot_data$ligand_to_rank,
                                     levels = unique(plot_data$ligand_to_rank))

  ggplot(plot_data, aes(x = .data$value, y = .data$ligand_to_rank, color = .data$variable)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(x = "Mean MOON score", y = NULL, color = NULL)
}


# ------------------------------------------------------------------------------
# Figure 2C: scoring-network characteristics
# ------------------------------------------------------------------------------

load_scoring_network_characteristics <- function(
    path = SNAPSHOT$scoring_network_characteristics) {
  if (is.na(path)) {
    return(NULL)
  }

  load_single_rdata_object(path)
}

summarise_scoring_network_characteristics <- function(characteristics) {
  as.data.frame(sapply(characteristics[c("nedges", "nnodes", "nTFs")], unlist))
}

plot_panel_2c_boxplots <- function(characteristics_df) {
  plot_data <- reshape2::melt(
    characteristics_df,
    measure.vars = names(characteristics_df),
    variable.name = "variable",
    value.name = "value"
  )

  ggplot(plot_data, aes(x = .data$variable, y = .data$value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
    facet_wrap(~variable, scales = "free", nrow = 1) +
    theme_minimal() +
    labs(
      title = "Moon scoring network characteristics",
      x = NULL,
      y = "Count"
    )
}

plot_panel_2c_scatter <- function(characteristics_df, include_marginals = WRITE_FIGURES) {
  base_plot <- ggplot(characteristics_df, aes(x = .data$nedges,
                                             y = .data$nnodes,
                                             size = .data$nTFs / 10)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(
      title = "Moon scoring network characteristics",
      x = "Number of edges",
      y = "Number of nodes",
      size = "TFs / 10"
    )

  if (!include_marginals) {
    return(base_plot)
  }

  ggExtra::ggMarginal(base_plot, type = "density", fill = "lightgrey")
}


# ------------------------------------------------------------------------------
# Figure 2D / Supplementary S2 network exports
# ------------------------------------------------------------------------------

export_moon_scoring_network <- function(experiment_id,
                                        ligand,
                                        zscores,
                                        meta_network,
                                        moon_res_list,
                                        output_prefix,
                                        keep_upstream_node_peers = FALSE) {
  ids <- moon_experiment_ids(moon_res_list)
  moon_res_experiment <- moon_res_list[[match(experiment_id, ids)]]

  if (is.null(moon_res_experiment)) {
    stop("Experiment not found in MOON cache: ", experiment_id)
  }

  expressed_genes <- zscores[, experiment_id, drop = FALSE]
  expressed_genes <- expressed_genes[complete.cases(expressed_genes), , drop = FALSE]
  expressed_genes <- setNames(expressed_genes[, 1], nm = row.names(expressed_genes))

  moon_res_experiment <- moon_res_experiment[, c(4, 2, 3)]
  names(moon_res_experiment) <- c("source", "score", "level")

  meta_network_filtered <- cosmosR:::filter_pkn_expressed_genes(
    expressed_genes_entrez = names(expressed_genes),
    meta_network
  )

  moon_scoring_network <- get_moon_scoring_network(
    upstream_node = ligand,
    meta_network = meta_network_filtered,
    moon_scores = moon_res_experiment,
    keep_upstream_node_peers = keep_upstream_node_peers
  )

  names(moon_scoring_network$SIF)[2] <- "sign"
  names(moon_scoring_network$ATT)[2] <- "moon_score"

  readr::write_csv(moon_scoring_network$SIF, paste0(output_prefix, "_SIF.csv"))
  readr::write_csv(moon_scoring_network$ATT, paste0(output_prefix, "_ATT.csv"))

  invisible(moon_scoring_network)
}

export_panel_2d_networks <- function(zscores, meta_network, moon_res_list) {
  examples <- data.frame(
    panel = c("Figure2D", "Figure2D", "Figure2D", "SupplementaryS2", "SupplementaryS2"),
    ligand = c("OSM", "IL6", "IL13", "TGFB1", "EGF"),
    experiment_id = c(
      "OSM@Condition:HCT116@GSE53295.MicroArray.GPL6480",
      "IL6@Condition:HuH-7@E-MTAB-4570.MicroArray.HTA-2_0",
      "IL13@Condition:Blood (PBMC)&Duration:24h@GSE79027.RNASeq.SRP071332_GRCh38",
      "TGFB1@Condition:HaCat&Duration:2h@E-MTAB-265.MicroArray.log_rsn_normalized",
      "EGF&Duration:60min@Condition:HCT116@GSE94374.RNASeq.SRP098688_GRCh38"
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(examples))) {
    output_prefix <- file.path(
      "results",
      paste(examples$panel[i], SNAPSHOT$label, examples$ligand[i], sep = "_")
    )

    export_moon_scoring_network(
      experiment_id = examples$experiment_id[i],
      ligand = examples$ligand[i],
      zscores = zscores,
      meta_network = meta_network,
      moon_res_list = moon_res_list,
      output_prefix = output_prefix,
      keep_upstream_node_peers = FALSE
    )
  }
}


# ------------------------------------------------------------------------------
# Optional original-vs-corrected comparison
# ------------------------------------------------------------------------------

compare_snapshot_ligand_tables <- function() {
  original <- as.data.frame(
    readr::read_csv(SNAPSHOTS$paper_original$ligand_table, show_col_types = FALSE)
  )
  corrected <- as.data.frame(
    readr::read_csv(SNAPSHOTS$corrected_20240319$ligand_table, show_col_types = FALSE)
  )

  original_summary <- original %>%
    group_by(.data$source) %>%
    summarise(mean_score_original = mean(.data$score, na.rm = TRUE), .groups = "drop")

  corrected_summary <- corrected %>%
    group_by(.data$source) %>%
    summarise(mean_score_corrected_20240319 = mean(.data$score, na.rm = TRUE),
              .groups = "drop")

  comparison <- full_join(original_summary, corrected_summary, by = "source") %>%
    mutate(delta_corrected_minus_original =
             .data$mean_score_corrected_20240319 - .data$mean_score_original) %>%
    arrange(desc(abs(.data$delta_corrected_minus_original)))

  print(head(comparison, 20))
  invisible(comparison)
}


# ------------------------------------------------------------------------------
# Expensive recomputation path
# ------------------------------------------------------------------------------

recompute_tf_activities <- function(zscores, collectri,
                                    output = "results/citosig_TF_activties.RData") {
  TF_activities <- apply(zscores, 2, function(x, collectri) {
    TF_act <- decoupleR::run_ulm(na.omit(x), network = collectri)
    TF_act[, c(2, 4)]
  }, collectri = collectri)

  for (i in seq_along(TF_activities)) {
    names(TF_activities[[i]])[2] <- names(TF_activities)[i]
  }

  save(TF_activities, file = output)
  invisible(TF_activities)
}

load_or_recompute_tf_activities <- function(zscores, collectri) {
  if (RECOMPUTE_TF) {
    return(recompute_tf_activities(zscores, collectri))
  }

  load_single_rdata_object("results/citosig_TF_activties.RData")
}

run_moon_for_experiment <- function(experiment_id,
                                    RNA_input,
                                    TF_input,
                                    meta_network,
                                    collectri,
                                    n_steps = 10) {
  RNA_input <- RNA_input[!is.na(RNA_input)]

  meta_network_filtered <- meta_network_cleanup(meta_network)
  meta_network_filtered <- cosmosR:::filter_pkn_expressed_genes(
    names(RNA_input),
    meta_pkn = meta_network_filtered
  )

  TF_input <- TF_input[!is.na(TF_input)]
  TF_input_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(
    TF_input,
    meta_network_filtered
  )
  meta_network_filtered <- cosmosR:::keep_observable_neighbours(
    meta_network_filtered,
    n_steps = n_steps,
    observed_nodes = names(TF_input_filtered)
  )

  meta_network_compressed_list <- compress_same_children(
    meta_network_filtered,
    sig_input = c(0),
    metab_input = TF_input_filtered
  )

  meta_network_compressed <- meta_network_cleanup(
    meta_network_compressed_list$compressed_network
  )
  meta_network_compressed_to_run <- meta_network_compressed

  before <- 1
  after <- 0
  iteration <- 1

  while (before != after && iteration < 10) {
    before <- nrow(meta_network_compressed_to_run)
    moon_res <- cosmosR::moon(
      downstream_input = TF_input_filtered,
      meta_network = meta_network_compressed_to_run,
      n_layers = n_steps,
      statistic = "ulm"
    )

    meta_network_compressed_to_run <- filter_incohrent_TF_target(
      moon_res,
      collectri,
      meta_network_compressed_to_run,
      RNA_input
    )

    after <- nrow(meta_network_compressed_to_run)
    iteration <- iteration + 1
  }

  message(
    experiment_id,
    ": ",
    ifelse(iteration < 10,
           paste0("converged after ", iteration - 1, " iterations"),
           paste0("interrupted after ", iteration, " iterations"))
  )

  decompress_moon_result(
    moon_res,
    meta_network_compressed_list,
    meta_network_compressed_to_run
  )
}

recompute_moon_cache <- function(zscores,
                                 TF_activities,
                                 meta_network,
                                 collectri,
                                 output = SNAPSHOT$moon_cache,
                                 n_steps = 10) {
  TF_activities_df <- TF_activities %>%
    reduce(full_join, by = "source") %>%
    as.data.frame()

  moon_res_list <- vector("list", ncol(zscores))

  for (experiment_number in seq_len(ncol(zscores))) {
    experiment_id <- names(zscores)[experiment_number]
    message("Running MOON for experiment ", experiment_number, "/", ncol(zscores),
            ": ", experiment_id)

    RNA_input <- zscores[, experiment_number]
    names(RNA_input) <- row.names(zscores)

    TF_input <- TF_activities_df[, experiment_number + 1]
    names(TF_input) <- TF_activities_df[, 1]

    moon_res <- run_moon_for_experiment(
      experiment_id = experiment_id,
      RNA_input = RNA_input,
      TF_input = TF_input,
      meta_network = meta_network,
      collectri = collectri,
      n_steps = n_steps
    )
    names(moon_res)[2] <- experiment_id
    moon_res_list[[experiment_number]] <- moon_res
  }

  save(moon_res_list, file = output)
  invisible(moon_res_list)
}


# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

moon_res_list <- load_moon_cache()
applied_ligands <- read_applied_ligand_table()
panel_2a_data <- prepare_panel_2a_data(applied_ligands, moon_res_list)

panel_2a_plot <- plot_panel_2a(panel_2a_data)
save_panel(panel_2a_plot, "moon_scores_ranked.pdf", width = 12, height = 6)

score_matrix <- build_score_matrix(moon_res_list)
panel_2b_summary <- build_panel_2b_summary(score_matrix, panel_2a_data)

panel_2b_quantiles_plot <- plot_panel_2b_quantiles(panel_2b_summary)
save_panel(panel_2b_quantiles_plot, "ligand_quantiles_long.pdf", width = 7, height = 9)

panel_2b_scores_plot <- plot_panel_2b_scores(panel_2b_summary)
save_panel(panel_2b_scores_plot, "ligand_score_long.pdf", width = 7, height = 9)

moon_scoring_network_characteristics <- load_scoring_network_characteristics()
if (!is.null(moon_scoring_network_characteristics)) {
  panel_2c_data <- summarise_scoring_network_characteristics(
    moon_scoring_network_characteristics
  )

  panel_2c_boxplot <- plot_panel_2c_boxplots(panel_2c_data)
  save_panel(panel_2c_boxplot, "moon_scoring_network_characteristics.pdf",
             width = 7, height = 4)

  panel_2c_scatter <- plot_panel_2c_scatter(panel_2c_data)
  save_panel(panel_2c_scatter, "moon_scoring_network_characteristics_scatter.pdf",
             width = 6, height = 5)
} else {
  message("Skipping Figure 2C for snapshot ", SNAPSHOT$label,
          ": no matching scoring-network characteristic cache is configured.")
}

message("Figure 2A ligands: ", length(unique(panel_2a_data$source)))
message("Figure 2A ligand-score rows: ", nrow(panel_2a_data))
message("Figure 2B ligands with higher true-experiment quantile: ",
        sum(panel_2b_summary$mean_quantile_in_trueExp >
              panel_2b_summary$mean_quantile_in_falseExp, na.rm = TRUE),
        "/", nrow(panel_2b_summary))
if (exists("panel_2c_data")) {
  message("Figure 2C median edges/nodes/TFs: ",
          paste(
            c(
              median(panel_2c_data$nedges, na.rm = TRUE),
              median(panel_2c_data$nnodes, na.rm = TRUE),
              median(panel_2c_data$nTFs, na.rm = TRUE)
            ),
            collapse = "/"
          ))
}

if (COMPARE_CORRECTED) {
  snapshot_comparison <- compare_snapshot_ligand_tables()
}

if (EXPORT_PANEL_D_NETWORKS) {
  zscores <- read_zscore_matrix()
  meta_network <- load_meta_network()
  export_panel_2d_networks(zscores, meta_network, moon_res_list)
}

if (RECOMPUTE_TF || RECOMPUTE_MOON) {
  zscores <- read_zscore_matrix()
  collectri <- load_collectri()
  meta_network <- load_meta_network()

  network_nodes <- unique(c(meta_network$source, meta_network$target))
  collectri <- collectri[collectri$source %in% network_nodes, ]
  collectri <- collectri[collectri$target %in% row.names(zscores), ]

  TF_activities <- load_or_recompute_tf_activities(zscores, collectri)

  if (RECOMPUTE_MOON) {
    moon_res_list <- recompute_moon_cache(
      zscores = zscores,
      TF_activities = TF_activities,
      meta_network = meta_network,
      collectri = collectri
    )
  }
}

if (RECOMPUTE_SCORING_NETWORKS) {
  stop("RECOMPUTE_SCORING_NETWORKS is intentionally not automated here yet. ",
       "Use export_panel_2d_networks() for paper examples, or port the legacy ",
       "all-ligand scoring-network loop once the desired exemplar set is fixed.")
}
