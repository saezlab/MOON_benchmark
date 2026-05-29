# Compare Figure 2A applied-ligand MOON scores between two CytoSig snapshots.
#
# The default comparison is the frozen paper snapshot versus the 2024-03-19
# corrected-network rerun. The plot mirrors the Figure 2A ligand-ranked
# box/jitter layout, including the same legacy CytoSig target mapping used by
# scripts/citosig_moon.R for moon_scores_ranked.pdf, but the y-axis is the
# performance delta:
#
#   delta = sign * (updated_score - original_score)
#
# For activating perturbations, this is simply updated - original. Positive
# values mean the updated snapshot moved the applied ligand in the expected
# direction; negative values mean it moved away from the expected direction.
#
# Usage:
#   Rscript scripts/plot_citosig_snapshot_delta.R
#   Rscript scripts/plot_citosig_snapshot_delta.R \
#     results/citosig_moon_activties.RData \
#     results/citosig_moon_activties_20240319.RData \
#     results/figures/moon_score_delta_original_vs_20240319.pdf

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

ORIGINAL_CACHE <- if (length(args) >= 1) {
  args[[1]]
} else {
  "results/citosig_moon_activties.RData"
}
UPDATED_CACHE <- if (length(args) >= 2) {
  args[[2]]
} else {
  "results/citosig_moon_activties_20240319.RData"
}
OUTPUT_PDF <- if (length(args) >= 3) {
  args[[3]]
} else {
  "results/figures/moon_score_delta_original_vs_20240319.pdf"
}
METADATA_PATH <- if (length(args) >= 4) {
  args[[4]]
} else {
  "data/cytosig/zscore_meta_clean_renamed.csv"
}
OUTPUT_CSV <- if (length(args) >= 5) {
  args[[5]]
} else {
  sub("\\.pdf$", ".csv", OUTPUT_PDF)
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

moon_experiment_ids <- function(moon_res_list) {
  vapply(moon_res_list, function(x) names(x)[2], character(1))
}

map_legacy_figure_target <- function(experiment_ids, metadata) {
  metadata$cytosig_target <- sub("[@&].*", "", metadata$id)
  target_mapping <- unique(metadata[, c("cytosig_target", "treatment")])
  target_mapping <- target_mapping[!is.na(target_mapping$cytosig_target), ]
  target_mapping <- setNames(target_mapping$treatment, target_mapping$cytosig_target)

  # This intentionally mirrors the original Figure 2A script. Experiments whose
  # IDs contain "&Dose" before the first "@" do not map to a target here, because
  # that is how the submitted figure's input table was produced.
  experiment_label <- sub("@.*", "", experiment_ids)
  unname(target_mapping[experiment_label])
}

extract_applied_ligand_scores <- function(moon_res_list, metadata, snapshot_label) {
  ids <- moon_experiment_ids(moon_res_list)
  metadata <- metadata[match(ids, metadata$id), , drop = FALSE]

  if (any(is.na(metadata$id))) {
    missing_metadata <- ids[is.na(metadata$id)]
    stop("Missing metadata for ", length(missing_metadata), " experiment IDs. ",
         "First missing ID: ", missing_metadata[[1]])
  }

  figure_targets <- map_legacy_figure_target(ids, metadata)
  scores <- vector("list", length(moon_res_list))

  for (i in seq_along(moon_res_list)) {
    moon_df <- moon_res_list[[i]]
    ligand <- figure_targets[[i]]

    if (is.na(ligand)) {
      scores[[i]] <- data.frame(
        id = ids[[i]],
        ligand = NA_character_,
        sign = metadata$sign[[i]],
        score = NA_real_,
        level = NA_integer_,
        snapshot = snapshot_label,
        figure_target_mapped = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }

    rows <- which(moon_df$source_original == ligand & moon_df$level != 0)

    if (length(rows) == 0) {
      scores[[i]] <- data.frame(
        id = ids[[i]],
        ligand = ligand,
        sign = metadata$sign[[i]],
        score = NA_real_,
        level = NA_integer_,
        snapshot = snapshot_label,
        figure_target_mapped = TRUE,
        stringsAsFactors = FALSE
      )
      next
    }

    scores[[i]] <- data.frame(
      id = ids[[i]],
      ligand = ligand,
      sign = metadata$sign[[i]],
      score = as.numeric(moon_df[[2]][rows[[1]]]),
      level = as.integer(moon_df$level[rows[[1]]]),
      snapshot = snapshot_label,
      figure_target_mapped = TRUE,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(scores)
}

build_delta_table <- function(original_scores, updated_scores) {
  compared <- full_join(
    original_scores,
    updated_scores,
    by = c("id", "ligand", "sign", "figure_target_mapped"),
    suffix = c("_original", "_updated")
  )

  compared$sign[is.na(compared$sign)] <- 1
  compared$performance_original <- compared$sign * compared$score_original
  compared$performance_updated <- compared$sign * compared$score_updated
  compared$performance_delta <- (
    compared$performance_updated - compared$performance_original
  )
  compared$level_delta <- compared$level_updated - compared$level_original
  compared$comparison_status <- dplyr::case_when(
    !compared$figure_target_mapped ~ "not_in_figure_2a_mapping",
    is.na(compared$score_original) & is.na(compared$score_updated) ~ "missing_both",
    is.na(compared$score_original) ~ "updated_only",
    is.na(compared$score_updated) ~ "original_only",
    TRUE ~ "shared"
  )

  compared
}

plot_delta_by_ligand <- function(delta_table) {
  plot_data <- delta_table %>%
    filter(.data$comparison_status == "shared") %>%
    group_by(.data$ligand) %>%
    mutate(
      mean_delta_by_ligand = mean(.data$performance_delta, na.rm = TRUE),
      median_delta_by_ligand = median(.data$performance_delta, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    arrange(desc(.data$mean_delta_by_ligand))

  plot_data$ligand <- factor(plot_data$ligand, levels = unique(plot_data$ligand))
  plot_data$mean_direction <- ifelse(
    plot_data$mean_delta_by_ligand >= 0,
    "improved",
    "degraded"
  )
  plot_data$level_delta <- factor(plot_data$level_delta)

  ggplot(plot_data, aes(x = .data$ligand,
                        y = .data$performance_delta,
                        group = .data$ligand)) +
    geom_boxplot(
      aes(fill = .data$mean_direction),
      coef = 6,
      color = "black",
      alpha = 0.75,
      outlier.shape = NA
    ) +
    geom_jitter(aes(color = .data$level_delta),
                width = 0.2,
                height = 0,
                alpha = 0.8,
                size = 1.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c(improved = "lightgreen", degraded = "tomato")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = NULL,
      y = "Performance delta (updated - original)",
      fill = "Mean delta",
      color = "Level delta",
      title = "Change in applied-ligand MOON score between snapshots"
    )
}

message("Original cache: ", ORIGINAL_CACHE)
message("Updated cache: ", UPDATED_CACHE)

metadata <- as.data.frame(readr::read_csv(METADATA_PATH, show_col_types = FALSE))
original_moon <- load_single_rdata_object(ORIGINAL_CACHE)
updated_moon <- load_single_rdata_object(UPDATED_CACHE)

original_scores <- extract_applied_ligand_scores(
  original_moon,
  metadata,
  snapshot_label = "original"
)
updated_scores <- extract_applied_ligand_scores(
  updated_moon,
  metadata,
  snapshot_label = "updated"
)

delta_table <- build_delta_table(original_scores, updated_scores)
plot_data <- delta_table %>% filter(.data$comparison_status == "shared")

message("Shared applied-ligand scores: ", nrow(plot_data))
message("Original-only applied-ligand scores: ",
        sum(delta_table$comparison_status == "original_only"))
message("Updated-only applied-ligand scores: ",
        sum(delta_table$comparison_status == "updated_only"))
message("Ligands in shared comparison: ", length(unique(plot_data$ligand)))
message("Ligands with positive mean delta: ",
        sum(tapply(plot_data$performance_delta, plot_data$ligand, mean,
                   na.rm = TRUE) > 0))

dir.create(dirname(OUTPUT_CSV), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(delta_table, OUTPUT_CSV)

delta_plot <- plot_delta_by_ligand(delta_table)
dir.create(dirname(OUTPUT_PDF), recursive = TRUE, showWarnings = FALSE)
ggsave(OUTPUT_PDF, delta_plot, width = 12, height = 6, units = "in")

message("Wrote delta table: ", OUTPUT_CSV)
message("Wrote delta plot: ", OUTPUT_PDF)
