#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(igraph)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
set.seed(20260802)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
script_dir <- dirname(script_path)
package_root <- dirname(script_dir)
workspace_root <- dirname(dirname(package_root))
platform_root <- file.path(workspace_root, "RareMol-AI_platform")
derived_dir <- file.path(package_root, "derived")
table_dir <- file.path(package_root, "tables")
private_verification_dir <- file.path(package_root, "verification", "private")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(private_verification_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  predictions = file.path(
    package_root, "source_data", "v1_6", "V1_6_PAIR_LEVEL_MODEL_PREDICTIONS.tsv"
  ),
  v16_metrics = file.path(
    package_root, "source_data", "v1_6", "V1_6_MODEL_STRATUM_METRICS.tsv"
  ),
  chembl_metadata = file.path(
    platform_root, "validation", "target_aware_phase16_chembl37_blind_20260729",
    "CHEMBL37_BLIND_PAIRS_WITH_IDENTITY_CATEGORIES.csv"
  ),
  chembl_clusters = file.path(
    platform_root, "validation", "target_aware_phase16_chembl37_blind_20260729",
    "PHASE16_COMBINED_PROTEIN_IDENTITY_CLUSTERS.csv"
  ),
  biolip_pairs = file.path(
    platform_root, "validation", "science_advances_v1_5_restricted_representation_20260801",
    "V1_5_RESTRICTED_PAIRS_WITH_COLD_CATEGORIES.tsv"
  ),
  biolip_outcomes = file.path(
    platform_root, "validation", "science_advances_v1_5_restricted_phase_b_20260801",
    "evaluation", "V1_5_PAIR_LEVEL_PERFORMANCE_LEDGER.tsv"
  ),
  biolip_clusters = file.path(
    platform_root, "validation", "science_advances_v1_5_restricted_representation_20260801",
    "V1_5_COMBINED_PROTEIN_IDENTITY_CLUSTERS.csv"
  ),
  development = file.path(
    platform_root, "data", "curated", "bindingdb_snapshot_temporal_202401_202607",
    "bindingdb_articles_202401_kd_qualifier_aware.csv"
  ),
  v15_metrics = file.path(
    platform_root, "validation", "science_advances_v1_5_restricted_phase_b_20260801",
    "evaluation", "V1_5_STRATUM_METRICS.tsv"
  ),
  v15_permutations = file.path(
    platform_root, "validation", "science_advances_v1_5_restricted_phase_b_20260801",
    "evaluation", "V1_5_PERMUTATION_NULL_METRICS.tsv"
  )
)

missing_paths <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0L) {
  stop("Missing required inputs: ", paste(missing_paths, collapse = ", "))
}

model_ids <- c(
  "hgb_molecule",
  "hgb_protein",
  "hgb_fusion_additive",
  "hgb_fusion_interaction",
  "phase15_esm2_aft_frozen"
)
model_labels <- c(
  hgb_molecule = "HGB ligand only",
  hgb_protein = "HGB protein only",
  hgb_fusion_additive = "HGB additive fusion",
  hgb_fusion_interaction = "HGB interaction fusion",
  phase15_esm2_aft_frozen = "Frozen ESM2-AFT"
)
model_families <- c(
  hgb_molecule = "Histogram gradient boosting",
  hgb_protein = "Histogram gradient boosting",
  hgb_fusion_additive = "Histogram gradient boosting",
  hgb_fusion_interaction = "Histogram gradient boosting",
  phase15_esm2_aft_frozen = "XGBoost accelerated failure time"
)
strata <- c(
  "overall",
  "drug_cold",
  "protein_identity_cold",
  "strict_scaffold_protein_identity_double_cold"
)
stratum_labels <- c(
  overall = "Overall",
  drug_cold = "Ligand cold",
  protein_identity_cold = "Target cold",
  strict_scaffold_protein_identity_double_cold = "Double cold"
)

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[keep], y[keep], method = method))
}

effective_clusters <- function(values) {
  proportions <- as.numeric(table(values)) / length(values)
  1 / sum(proportions^2)
}

first_token <- function(values) {
  values <- as.character(values)
  tokens <- strsplit(values, "[|;]")
  vapply(tokens, function(x) {
    x <- trimws(x[nzchar(trimws(x))])
    if (length(x) == 0L) NA_character_ else x[[1L]]
  }, character(1))
}

document_components <- function(document_values, prefix) {
  lists <- lapply(strsplit(as.character(document_values), "[|;]"), function(values) {
    values <- sort(unique(trimws(values[nzchar(trimws(values))])))
    if (length(values) == 0L) NA_character_ else values
  })
  empty <- which(vapply(lists, function(x) all(is.na(x)), logical(1)))
  if (length(empty) > 0L) {
    for (index in empty) lists[[index]] <- paste0(prefix, "_NO_DOCUMENT_", index)
  }
  vertices <- unique(unlist(lists, use.names = FALSE))
  edges <- rbindlist(lapply(lists, function(values) {
    if (length(values) <= 1L) return(NULL)
    data.table(from = values[[1L]], to = values[-1L])
  }), fill = TRUE)
  if (nrow(edges) > 0L) {
    graph <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = vertices))
    membership <- components(graph)$membership
  } else {
    membership <- setNames(seq_along(vertices), vertices)
  }
  primary <- vapply(lists, `[[`, character(1), 1L)
  data.table(
    primary_document_id = primary,
    document_component_id = paste0(prefix, "_DOC_COMPONENT_", unname(membership[primary])),
    document_count = lengths(lists)
  )
}

predictions <- fread(paths$predictions)
if (anyDuplicated(predictions$pair_id)) stop("v1.6 pair identifiers are not unique.")

chembl_meta <- fread(paths$chembl_metadata)
setnames(chembl_meta, "row_id", "pair_id")
chembl_doc <- document_components(chembl_meta$document_chembl_ids, "CHEMBL37")
chembl_meta <- cbind(chembl_meta, chembl_doc)
chembl_meta[, source := "ChEMBL37"]
chembl_meta[, assay_cluster_id := first_token(assay_chembl_ids)]
chembl_meta[, protein_length_derived := as.integer(protein_length)]
chembl_meta[, provenance_candidate_count := as.integer(measurement_count)]
chembl_meta[, provenance_label_range := as.numeric(pkd_max - pkd_min)]

biolip_pairs <- fread(paths$biolip_pairs)
biolip_outcomes <- fread(paths$biolip_outcomes)
biolip_meta <- merge(
  biolip_pairs,
  biolip_outcomes[, .(
    frozen_pair_id, pkd_outcome = pkd, pkd_min, pkd_max, pkd_range,
    exact_kd_measurement_count, contributing_candidate_count,
    document_block_id, raw_affinity_text_retained
  )],
  by = "frozen_pair_id", all.x = TRUE, sort = FALSE
)
setnames(biolip_meta, "frozen_pair_id", "pair_id")
biolip_meta[, source := "BioLiP/MOAD"]
biolip_meta[, primary_document_id := as.character(document_block_id)]
biolip_meta[, document_component_id := as.character(document_block_id)]
biolip_meta[, document_count := 1L]
biolip_meta[, assay_cluster_id := first_token(pdb_ids)]
biolip_meta[, protein_length_derived := nchar(protein_sequence)]
biolip_meta[, provenance_candidate_count := as.integer(contributing_candidate_count)]
biolip_meta[, provenance_label_range := as.numeric(pkd_range)]

metadata_columns <- c(
  "pair_id", "source", "frozen_scaffold_id", "protein_identity_cluster_id",
  "primary_document_id", "document_component_id", "document_count",
  "assay_cluster_id", "protein_length_derived", "provenance_candidate_count",
  "provenance_label_range", "pkd_min", "pkd_max"
)
chembl_link <- chembl_meta[, ..metadata_columns]
biolip_link <- biolip_meta[, ..metadata_columns]
metadata <- rbindlist(list(chembl_link, biolip_link), use.names = TRUE, fill = TRUE)
analysis <- merge(predictions, metadata, by = c("pair_id", "source"), all.x = TRUE, sort = FALSE)
if (nrow(analysis) != nrow(predictions)) stop("Metadata merge changed pair membership.")
required_linkage <- c("frozen_scaffold_id", "protein_identity_cluster_id", "document_component_id")
if (analysis[, anyNA(.SD), .SDcols = required_linkage]) stop("Incomplete dependence linkage.")
analysis[, exact_pair_key := paste(canonical_smiles, protein_sha256, sep = "|")]
analysis[, ligand_identity := canonical_smiles]
analysis[, document_component_size := .N, by = .(source, document_component_id)]
analysis[, protein_identity_size := .N, by = .(source, protein_sha256)]
analysis[, protein_cluster_size := .N, by = .(source, protein_identity_cluster_id)]
analysis[, scaffold_size := .N, by = .(source, frozen_scaffold_id)]
analysis[, affinity_regime := cut(
  pkd,
  breaks = c(-Inf, 5, 6, 7, 8, Inf),
  labels = c("<5", "5-<6", "6-<7", "7-<8", ">=8"),
  right = FALSE
)]

subset_stratum <- function(frame, stratum) {
  if (stratum == "overall") frame else frame[external_cold_category == stratum]
}

cluster_vcov <- function(fit, frame, cluster_cols) {
  if (length(cluster_cols) == 0L) return(vcovHC(fit, type = "HC1"))
  vcovCL(
    fit,
    cluster = lapply(cluster_cols, function(column) frame[[column]]),
    type = "HC1",
    fix = TRUE,
    multi0 = TRUE
  )
}

rank_interval <- function(frame, prediction_column, cluster_cols, z_value = 1.96) {
  work <- data.frame(
    truth_rank = as.numeric(scale(rank(frame$pkd, ties.method = "average"))),
    prediction_rank = as.numeric(scale(rank(frame[[prediction_column]], ties.method = "average")))
  )
  fit <- lm(truth_rank ~ prediction_rank, data = work)
  covariance <- cluster_vcov(fit, frame, cluster_cols)
  estimate <- unname(coef(fit)[["prediction_rank"]])
  standard_error <- sqrt(max(0, diag(covariance)[["prediction_rank"]]))
  c(
    estimate = estimate,
    se = standard_error,
    lower = max(-1, estimate - z_value * standard_error),
    upper = min(1, estimate + z_value * standard_error)
  )
}

residual_rank_correlation <- function(frame, prediction_column, group_column) {
  groups <- frame[[group_column]]
  group_sizes <- ave(rep(1, nrow(frame)), groups, FUN = sum)
  keep <- group_sizes >= 2L
  if (sum(keep) < 3L) return(c(rho = NA_real_, n = 0, groups = 0))
  truth_rank <- rank(frame$pkd, ties.method = "average")
  prediction_rank <- rank(frame[[prediction_column]], ties.method = "average")
  truth_residual <- truth_rank - ave(truth_rank, groups, FUN = mean)
  prediction_residual <- prediction_rank - ave(prediction_rank, groups, FUN = mean)
  c(
    rho = safe_cor(truth_residual[keep], prediction_residual[keep], method = "pearson"),
    n = sum(keep),
    groups = uniqueN(groups[keep])
  )
}

rank_covariance_decomposition <- function(frame, prediction_column, group_column) {
  truth <- as.numeric(scale(rank(frame$pkd, ties.method = "average"), center = TRUE, scale = FALSE))
  prediction <- as.numeric(scale(
    rank(frame[[prediction_column]], ties.method = "average"), center = TRUE, scale = FALSE
  ))
  groups <- frame[[group_column]]
  truth_group_mean <- ave(truth, groups, FUN = mean)
  prediction_group_mean <- ave(prediction, groups, FUN = mean)
  total_covariance <- mean(truth * prediction)
  between_covariance <- mean(truth_group_mean * prediction_group_mean)
  within_covariance <- mean(
    (truth - truth_group_mean) * (prediction - prediction_group_mean)
  )
  c(
    total_covariance = total_covariance,
    between_covariance = between_covariance,
    within_covariance = within_covariance,
    between_fraction = ifelse(
      abs(total_covariance) < 1e-12, NA_real_, between_covariance / total_covariance
    )
  )
}

macro_within_group_rank <- function(frame, prediction_column, group_column) {
  group_sizes <- frame[, .N, by = group_column]
  eligible <- group_sizes[N >= 3L][[group_column]]
  if (length(eligible) == 0L) {
    return(c(groups = 0, n = 0, macro_rho = NA_real_, weighted_rho = NA_real_))
  }
  summaries <- frame[get(group_column) %in% eligible, .(
    n = .N,
    rho = safe_cor(pkd, get(prediction_column))
  ), by = group_column]
  summaries <- summaries[is.finite(rho)]
  c(
    groups = nrow(summaries),
    n = sum(summaries$n),
    macro_rho = mean(summaries$rho),
    weighted_rho = weighted.mean(summaries$rho, summaries$n)
  )
}

rank_rows <- list()
row_index <- 0L
for (source_name in sort(unique(analysis$source))) {
  source_frame <- analysis[source == source_name]
  for (stratum_name in strata) {
    frame <- subset_stratum(source_frame, stratum_name)
    if (nrow(frame) < 20L) next
    for (model_id in model_ids) {
      prediction_column <- paste0("prediction__", model_id)
      pair_ci <- rank_interval(frame, prediction_column, character())
      document_ci <- rank_interval(frame, prediction_column, "document_component_id")
      scaffold_ci <- rank_interval(frame, prediction_column, "frozen_scaffold_id")
      protein_cluster_ci <- rank_interval(
        frame, prediction_column, "protein_identity_cluster_id"
      )
      multiway_ci <- rank_interval(
        frame,
        prediction_column,
        c("frozen_scaffold_id", "protein_identity_cluster_id", "document_component_id")
      )
      dependence_envelope_se <- max(c(
        document_ci[["se"]], scaffold_ci[["se"]],
        protein_cluster_ci[["se"]], multiway_ci[["se"]]
      ), na.rm = TRUE)
      dependence_envelope_low <- max(
        -1, pair_ci[["estimate"]] - 1.96 * dependence_envelope_se
      )
      dependence_envelope_high <- min(
        1, pair_ci[["estimate"]] + 1.96 * dependence_envelope_se
      )
      document_residual <- residual_rank_correlation(
        frame, prediction_column, "document_component_id"
      )
      protein_residual <- residual_rank_correlation(frame, prediction_column, "protein_sha256")
      protein_cluster_residual <- residual_rank_correlation(
        frame, prediction_column, "protein_identity_cluster_id"
      )
      scaffold_residual <- residual_rank_correlation(
        frame, prediction_column, "frozen_scaffold_id"
      )
      protein_decomposition <- rank_covariance_decomposition(
        frame, prediction_column, "protein_sha256"
      )
      document_decomposition <- rank_covariance_decomposition(
        frame, prediction_column, "document_component_id"
      )
      protein_macro <- macro_within_group_rank(frame, prediction_column, "protein_sha256")
      row_index <- row_index + 1L
      rank_rows[[row_index]] <- data.table(
        source = source_name,
        stratum = stratum_name,
        stratum_label = stratum_labels[[stratum_name]],
        model_id = model_id,
        model_label = model_labels[[model_id]],
        algorithm_family = model_families[[model_id]],
        n = nrow(frame),
        spearman = pair_ci[["estimate"]],
        pair_hc1_ci95_low = pair_ci[["lower"]],
        pair_hc1_ci95_high = pair_ci[["upper"]],
        document_cluster_se = document_ci[["se"]],
        document_cluster_ci95_low = document_ci[["lower"]],
        document_cluster_ci95_high = document_ci[["upper"]],
        scaffold_cluster_se = scaffold_ci[["se"]],
        scaffold_cluster_ci95_low = scaffold_ci[["lower"]],
        scaffold_cluster_ci95_high = scaffold_ci[["upper"]],
        protein_cluster_se = protein_cluster_ci[["se"]],
        protein_cluster_ci95_low = protein_cluster_ci[["lower"]],
        protein_cluster_ci95_high = protein_cluster_ci[["upper"]],
        three_way_cluster_se = multiway_ci[["se"]],
        three_way_cluster_ci95_low = multiway_ci[["lower"]],
        three_way_cluster_ci95_high = multiway_ci[["upper"]],
        dependence_envelope_se = dependence_envelope_se,
        dependence_envelope_ci95_low = dependence_envelope_low,
        dependence_envelope_ci95_high = dependence_envelope_high,
        document_residual_rank_rho = document_residual[["rho"]],
        document_residual_eligible_n = as.integer(document_residual[["n"]]),
        protein_residual_rank_rho = protein_residual[["rho"]],
        protein_residual_eligible_n = as.integer(protein_residual[["n"]]),
        protein_cluster_residual_rank_rho = protein_cluster_residual[["rho"]],
        scaffold_residual_rank_rho = scaffold_residual[["rho"]],
        exact_protein_between_covariance_fraction = protein_decomposition[["between_fraction"]],
        document_between_covariance_fraction = document_decomposition[["between_fraction"]],
        exact_protein_macro_groups = as.integer(protein_macro[["groups"]]),
        exact_protein_macro_n = as.integer(protein_macro[["n"]]),
        exact_protein_macro_spearman = protein_macro[["macro_rho"]],
        exact_protein_n_weighted_spearman = protein_macro[["weighted_rho"]],
        interpretation = "Post-unblinding decomposition; rank CIs use rank-regression sandwich inference"
      )
    }
  }
}
rank_decomposition <- rbindlist(rank_rows, use.names = TRUE, fill = TRUE)
fwrite(
  rank_decomposition,
  file.path(derived_dir, "V1_7_RANK_SIGNAL_DECOMPOSITION.tsv"),
  sep = "\t"
)

bootstrap_overall_rank <- function(frame, replicates = 1000L, seed = 20260802L) {
  set.seed(seed)
  grouped_rows <- split(seq_len(nrow(frame)), frame$document_component_id)
  group_names <- names(grouped_rows)
  observed <- vapply(model_ids, function(model_id) {
    safe_cor(frame$pkd, frame[[paste0("prediction__", model_id)]])
  }, numeric(1))
  bootstrap <- matrix(NA_real_, nrow = replicates, ncol = length(model_ids))
  colnames(bootstrap) <- model_ids
  for (replicate in seq_len(replicates)) {
    sampled_groups <- sample(group_names, length(group_names), replace = TRUE)
    sampled_rows <- unlist(grouped_rows[sampled_groups], use.names = FALSE)
    truth <- frame$pkd[sampled_rows]
    for (model_index in seq_along(model_ids)) {
      prediction <- frame[[paste0("prediction__", model_ids[[model_index]])]][sampled_rows]
      bootstrap[replicate, model_index] <- safe_cor(truth, prediction)
    }
  }
  standard_errors <- apply(bootstrap, 2, sd, na.rm = TRUE)
  t_statistics <- sweep(bootstrap, 2, observed, "-")
  t_statistics <- sweep(t_statistics, 2, standard_errors, "/")
  simultaneous_quantile <- as.numeric(quantile(
    apply(abs(t_statistics), 1, max, na.rm = TRUE), 0.95, na.rm = TRUE
  ))
  rbindlist(lapply(seq_along(model_ids), function(index) {
    values <- bootstrap[, index]
    data.table(
      model_id = model_ids[[index]],
      model_label = model_labels[[model_ids[[index]]]],
      estimate = observed[[index]],
      bootstrap_se = standard_errors[[index]],
      percentile_ci95_low = as.numeric(quantile(values, 0.025, na.rm = TRUE)),
      percentile_ci95_high = as.numeric(quantile(values, 0.975, na.rm = TRUE)),
      simultaneous_ci95_low = max(-1, observed[[index]] - simultaneous_quantile * standard_errors[[index]]),
      simultaneous_ci95_high = min(1, observed[[index]] + simultaneous_quantile * standard_errors[[index]]),
      simultaneous_max_t_quantile = simultaneous_quantile,
      replicates = replicates,
      bootstrap_unit = "connected_document_component",
      familywise_scope = "five_nonreference_specifications_within_source"
    )
  }))
}

overall_bootstrap <- rbindlist(lapply(sort(unique(analysis$source)), function(source_name) {
  result <- bootstrap_overall_rank(
    analysis[source == source_name],
    replicates = 1000L,
    seed = 20260802L + match(source_name, sort(unique(analysis$source)))
  )
  result[, source := source_name]
  result
}), use.names = TRUE, fill = TRUE)
setcolorder(overall_bootstrap, c("source", setdiff(names(overall_bootstrap), "source")))
fwrite(
  overall_bootstrap,
  file.path(derived_dir, "V1_7_OVERALL_DOCUMENT_BOOTSTRAP_RANK.tsv"),
  sep = "\t"
)

metric_interval <- function(frame, outcome, prediction_column = NULL, cluster_cols = character()) {
  if (is.null(prediction_column)) {
    fit <- lm(frame[[outcome]] ~ 1)
    coefficient <- "(Intercept)"
  } else {
    fit_frame <- data.frame(outcome = frame[[outcome]], prediction = frame[[prediction_column]])
    fit <- lm(outcome ~ prediction, data = fit_frame)
    coefficient <- "prediction"
  }
  covariance <- cluster_vcov(fit, frame, cluster_cols)
  estimate <- unname(coef(fit)[[coefficient]])
  standard_error <- sqrt(max(0, diag(covariance)[[coefficient]]))
  c(
    estimate = estimate,
    se = standard_error,
    lower = estimate - 1.96 * standard_error,
    upper = estimate + 1.96 * standard_error
  )
}

utility_rows <- list()
utility_index <- 0L
for (source_name in sort(unique(analysis$source))) {
  source_frame <- analysis[source == source_name]
  for (stratum_name in strata) {
    frame <- subset_stratum(source_frame, stratum_name)
    if (nrow(frame) < 20L) next
    for (model_id in model_ids) {
      prediction_column <- paste0("prediction__", model_id)
      work <- copy(frame)
      work[, prediction := get(prediction_column)]
      work[, delta_mae := abs(pkd - prediction) - abs(pkd - prediction__development_median)]
      work[, signed_error := prediction - pkd]
      delta_document <- metric_interval(work, "delta_mae", cluster_cols = "document_component_id")
      delta_scaffold <- metric_interval(work, "delta_mae", cluster_cols = "frozen_scaffold_id")
      delta_protein <- metric_interval(
        work, "delta_mae", cluster_cols = "protein_identity_cluster_id"
      )
      delta_multiway <- metric_interval(
        work, "delta_mae",
        cluster_cols = c("frozen_scaffold_id", "protein_identity_cluster_id", "document_component_id")
      )
      delta_envelope_se <- max(c(
        delta_document[["se"]], delta_scaffold[["se"]],
        delta_protein[["se"]], delta_multiway[["se"]]
      ), na.rm = TRUE)
      bias_multiway <- metric_interval(
        work, "signed_error",
        cluster_cols = c("frozen_scaffold_id", "protein_identity_cluster_id", "document_component_id")
      )
      calibration_multiway <- metric_interval(
        work, "pkd", prediction_column = "prediction",
        cluster_cols = c("frozen_scaffold_id", "protein_identity_cluster_id", "document_component_id")
      )
      coverage <- coverage_multiway <- c(estimate = NA, se = NA, lower = NA, upper = NA)
      coverage_scaffold <- coverage_protein <- c(
        estimate = NA, se = NA, lower = NA, upper = NA
      )
      coverage_envelope_se <- NA_real_
      if (model_id == "phase15_esm2_aft_frozen") {
        work[, interval_covered := as.numeric(abs(pkd - prediction) <= phase15_interval_half_width_pkd)]
        coverage <- metric_interval(work, "interval_covered", cluster_cols = "document_component_id")
        coverage_scaffold <- metric_interval(
          work, "interval_covered", cluster_cols = "frozen_scaffold_id"
        )
        coverage_protein <- metric_interval(
          work, "interval_covered", cluster_cols = "protein_identity_cluster_id"
        )
        coverage_multiway <- metric_interval(
          work, "interval_covered",
          cluster_cols = c("frozen_scaffold_id", "protein_identity_cluster_id", "document_component_id")
        )
        coverage_envelope_se <- max(c(
          coverage[["se"]], coverage_scaffold[["se"]],
          coverage_protein[["se"]], coverage_multiway[["se"]]
        ), na.rm = TRUE)
      }
      utility_index <- utility_index + 1L
      utility_rows[[utility_index]] <- data.table(
        source = source_name,
        stratum = stratum_name,
        stratum_label = stratum_labels[[stratum_name]],
        model_id = model_id,
        model_label = model_labels[[model_id]],
        n = nrow(work),
        mae = mean(abs(work$pkd - work$prediction)),
        baseline_mae = mean(abs(work$pkd - work$prediction__development_median)),
        delta_mae = mean(work$delta_mae),
        delta_mae_document_ci95_low = delta_document[["lower"]],
        delta_mae_document_ci95_high = delta_document[["upper"]],
        delta_mae_scaffold_ci95_low = delta_scaffold[["lower"]],
        delta_mae_scaffold_ci95_high = delta_scaffold[["upper"]],
        delta_mae_protein_cluster_ci95_low = delta_protein[["lower"]],
        delta_mae_protein_cluster_ci95_high = delta_protein[["upper"]],
        delta_mae_three_way_ci95_low = delta_multiway[["lower"]],
        delta_mae_three_way_ci95_high = delta_multiway[["upper"]],
        delta_mae_dependence_envelope_ci95_low = mean(work$delta_mae) - 1.96 * delta_envelope_se,
        delta_mae_dependence_envelope_ci95_high = mean(work$delta_mae) + 1.96 * delta_envelope_se,
        r2 = 1 - sum((work$pkd - work$prediction)^2) / sum((work$pkd - mean(work$pkd))^2),
        mean_bias = mean(work$signed_error),
        mean_bias_three_way_ci95_low = bias_multiway[["lower"]],
        mean_bias_three_way_ci95_high = bias_multiway[["upper"]],
        calibration_slope = calibration_multiway[["estimate"]],
        calibration_slope_three_way_ci95_low = calibration_multiway[["lower"]],
        calibration_slope_three_way_ci95_high = calibration_multiway[["upper"]],
        coverage_90 = coverage_multiway[["estimate"]],
        coverage_90_document_ci95_low = max(0, coverage[["lower"]]),
        coverage_90_document_ci95_high = min(1, coverage[["upper"]]),
        coverage_90_three_way_ci95_low = max(0, coverage_multiway[["lower"]]),
        coverage_90_three_way_ci95_high = min(1, coverage_multiway[["upper"]]),
        coverage_90_dependence_envelope_ci95_low = ifelse(
          is.na(coverage_envelope_se), NA_real_,
          max(0, coverage_multiway[["estimate"]] - 1.96 * coverage_envelope_se)
        ),
        coverage_90_dependence_envelope_ci95_high = ifelse(
          is.na(coverage_envelope_se), NA_real_,
          min(1, coverage_multiway[["estimate"]] + 1.96 * coverage_envelope_se)
        )
      )
    }
  }
}
utility <- rbindlist(utility_rows, use.names = TRUE, fill = TRUE)
fwrite(utility, file.path(derived_dir, "V1_7_UTILITY_CALIBRATION_INFERENCE.tsv"), sep = "\t")

largest_group_sensitivity <- function(frame, model_id, group_column, group_type) {
  prediction_column <- paste0("prediction__", model_id)
  group_counts <- frame[, .N, by = group_column][order(-N, get(group_column))]
  dominant_group <- group_counts[[group_column]][[1L]]
  reduced <- frame[get(group_column) != dominant_group]
  metrics_for <- function(data) {
    prediction <- data[[prediction_column]]
    c(
      n = nrow(data),
      spearman = safe_cor(data$pkd, prediction),
      delta_mae = mean(abs(data$pkd - prediction) - abs(data$pkd - data$prediction__development_median)),
      coverage = if (model_id == "phase15_esm2_aft_frozen") {
        mean(abs(data$pkd - prediction) <= data$phase15_interval_half_width_pkd)
      } else NA_real_
    )
  }
  full_metrics <- metrics_for(frame)
  reduced_metrics <- metrics_for(reduced)
  data.table(
    group_type = group_type,
    removed_group_id = as.character(dominant_group),
    removed_n = group_counts$N[[1L]],
    removed_fraction = group_counts$N[[1L]] / nrow(frame),
    effective_groups = effective_clusters(frame[[group_column]]),
    full_n = as.integer(full_metrics[["n"]]),
    reduced_n = as.integer(reduced_metrics[["n"]]),
    full_spearman = full_metrics[["spearman"]],
    reduced_spearman = reduced_metrics[["spearman"]],
    spearman_shift = reduced_metrics[["spearman"]] - full_metrics[["spearman"]],
    full_delta_mae = full_metrics[["delta_mae"]],
    reduced_delta_mae = reduced_metrics[["delta_mae"]],
    delta_mae_shift = reduced_metrics[["delta_mae"]] - full_metrics[["delta_mae"]],
    full_coverage_90 = full_metrics[["coverage"]],
    reduced_coverage_90 = reduced_metrics[["coverage"]]
  )
}

influence <- rbindlist(lapply(sort(unique(analysis$source)), function(source_name) {
  source_frame <- analysis[source == source_name]
  rbindlist(lapply(strata, function(stratum_name) {
    frame <- subset_stratum(source_frame, stratum_name)
    if (nrow(frame) < 20L) return(NULL)
    rbindlist(lapply(model_ids, function(model_id) {
      result <- rbindlist(list(
        largest_group_sensitivity(frame, model_id, "document_component_id", "document_component"),
        largest_group_sensitivity(frame, model_id, "protein_identity_cluster_id", "protein_identity_cluster"),
        largest_group_sensitivity(frame, model_id, "frozen_scaffold_id", "ligand_scaffold")
      ))
      result[, `:=`(
        source = source_name,
        stratum = stratum_name,
        stratum_label = stratum_labels[[stratum_name]],
        model_id = model_id,
        model_label = model_labels[[model_id]]
      )]
      result
    }))
  }))
}), use.names = TRUE, fill = TRUE)
setcolorder(influence, c(
  "source", "stratum", "stratum_label", "model_id", "model_label",
  setdiff(names(influence), c("source", "stratum", "stratum_label", "model_id", "model_label"))
))
fwrite(influence, file.path(derived_dir, "V1_7_DOMINANT_CLUSTER_INFLUENCE.tsv"), sep = "\t")

source_profile <- analysis[, {
  document_counts <- table(document_component_id)
  .(
    n = .N,
    unique_ligands = uniqueN(canonical_smiles),
    unique_scaffolds = uniqueN(frozen_scaffold_id),
    unique_proteins = uniqueN(protein_sha256),
    protein_identity_clusters = uniqueN(protein_identity_cluster_id),
    document_components = uniqueN(document_component_id),
    assay_clusters = uniqueN(assay_cluster_id, na.rm = TRUE),
    repeated_ligand_row_fraction = mean(duplicated(canonical_smiles) | duplicated(canonical_smiles, fromLast = TRUE)),
    repeated_protein_row_fraction = mean(duplicated(protein_sha256) | duplicated(protein_sha256, fromLast = TRUE)),
    repeated_document_row_fraction = mean(document_component_size > 1L),
    largest_document_component_n = max(document_counts),
    largest_document_component_fraction = max(document_counts) / .N,
    effective_document_components = effective_clusters(document_component_id),
    outcome_mean = mean(pkd),
    outcome_sd = sd(pkd),
    outcome_min = min(pkd),
    outcome_max = max(pkd),
    outcomes_below_2 = sum(pkd < 2),
    outcomes_above_12 = sum(pkd > 12),
    multi_measurement_fraction = mean(provenance_candidate_count > 1, na.rm = TRUE),
    nonzero_label_range_fraction = mean(provenance_label_range > 0, na.rm = TRUE)
  )
}, by = source]
fwrite(source_profile, file.path(derived_dir, "V1_7_SOURCE_DATA_PROFILE.tsv"), sep = "\t")

source_names <- sort(unique(analysis$source))
source_a <- analysis[source == source_names[[1L]]]
source_b <- analysis[source == source_names[[2L]]]
overlap <- data.table(
  entity = c("standardized_ligand", "ligand_scaffold", "exact_protein", "protein_identity_cluster", "exact_pair"),
  source_a = source_names[[1L]],
  source_b = source_names[[2L]],
  overlap_n = c(
    length(intersect(unique(source_a$canonical_smiles), unique(source_b$canonical_smiles))),
    length(intersect(unique(source_a$frozen_scaffold_id), unique(source_b$frozen_scaffold_id))),
    length(intersect(unique(source_a$protein_sha256), unique(source_b$protein_sha256))),
    length(intersect(unique(source_a$protein_identity_cluster_id), unique(source_b$protein_identity_cluster_id))),
    length(intersect(unique(source_a$exact_pair_key), unique(source_b$exact_pair_key)))
  ),
  interpretation = c(
    "Exact standardized ligand overlap",
    "Hash overlap; valid only because the same Bemis-Murcko hashing policy was used",
    "Exact sequence-hash overlap",
    "Cluster labels were generated in source-specific combined runs and are not cross-source identifiers",
    "Exact standardized ligand-protein pair overlap"
  )
)
fwrite(overlap, file.path(derived_dir, "V1_7_CROSS_SOURCE_OVERLAP.tsv"), sep = "\t")

weighted_cor <- function(x, y, weights) {
  keep <- is.finite(x) & is.finite(y) & is.finite(weights) & weights > 0
  x <- x[keep]
  y <- y[keep]
  weights <- weights[keep] / sum(weights[keep])
  x_centered <- x - sum(weights * x)
  y_centered <- y - sum(weights * y)
  sum(weights * x_centered * y_centered) /
    sqrt(sum(weights * x_centered^2) * sum(weights * y_centered^2))
}

analysis[, standardization_cell := paste(external_cold_category, affinity_regime, sep = "|")]
cell_support <- analysis[, .N, by = .(source, standardization_cell)]
common_cells <- Reduce(intersect, lapply(source_names, function(source_name) {
  cell_support[source == source_name & N > 0, standardization_cell]
}))
common_cell_support <- copy(cell_support[standardization_cell %in% common_cells])
common_cell_support[, source_fraction := N / sum(N), by = source]
target_distribution <- common_cell_support[,
  .(target_fraction = mean(source_fraction)), by = standardization_cell
]
target_distribution[, target_fraction := target_fraction / sum(target_fraction)]

standardized_rows <- rbindlist(lapply(source_names, function(source_name) {
  source_frame <- analysis[source == source_name & standardization_cell %in% common_cells]
  source_distribution <- source_frame[, .N, by = standardization_cell][, source_fraction := N / sum(N)]
  weights <- merge(source_distribution, target_distribution, by = "standardization_cell")
  weights[, cell_weight := target_fraction / source_fraction]
  source_frame <- merge(source_frame, weights[, .(standardization_cell, cell_weight)], by = "standardization_cell")
  rbindlist(lapply(model_ids, function(model_id) {
    prediction <- source_frame[[paste0("prediction__", model_id)]]
    truth_rank <- rank(source_frame$pkd, ties.method = "average")
    prediction_rank <- rank(prediction, ties.method = "average")
    data.table(
      source = source_name,
      model_id = model_id,
      model_label = model_labels[[model_id]],
      analysis = c("raw_common_support", "affinity_and_cold_poststratified"),
      n = nrow(source_frame),
      common_support_cells = length(common_cells),
      mae = c(
        mean(abs(source_frame$pkd - prediction)),
        weighted.mean(abs(source_frame$pkd - prediction), source_frame$cell_weight)
      ),
      delta_mae_vs_development_median = c(
        mean(abs(source_frame$pkd - prediction) - abs(source_frame$pkd - source_frame$prediction__development_median)),
        weighted.mean(
          abs(source_frame$pkd - prediction) - abs(source_frame$pkd - source_frame$prediction__development_median),
          source_frame$cell_weight
        )
      ),
      rank_correlation = c(
        safe_cor(source_frame$pkd, prediction),
        weighted_cor(truth_rank, prediction_rank, source_frame$cell_weight)
      ),
      limitation = "Poststratifies observed affinity regime and cold category only; target-family confounding remains because exact proteins do not overlap"
    )
  }))
}), use.names = TRUE, fill = TRUE)
fwrite(
  standardized_rows,
  file.path(derived_dir, "V1_7_SOURCE_POSTSTRATIFICATION_SENSITIVITY.tsv"),
  sep = "\t"
)

development <- fread(paths$development)
development <- development[censoring == "exact" & is.finite(pkd)]

regularized_lookup <- function(values, keys, minimum_n = 2L) {
  summary <- data.table(group_key = keys, value = values)[
    , .(median = median(value), n = .N), by = group_key
  ]
  summary[n < minimum_n, median := NA_real_]
  setNames(summary$median, summary$group_key)
}

entity_baseline_for_source <- function(source_name, cluster_path) {
  cluster_map <- fread(cluster_path)[, .(protein_sha256, protein_identity_cluster_id)]
  development_linked <- merge(
    development,
    cluster_map,
    by = "protein_sha256",
    all.x = TRUE,
    sort = FALSE
  )
  global_median <- median(development_linked$pkd)
  molecule_medians <- regularized_lookup(development_linked$pkd, development_linked$canonical_smiles)
  protein_medians <- regularized_lookup(development_linked$pkd, development_linked$protein_sha256)
  cluster_medians <- regularized_lookup(
    development_linked$pkd, development_linked$protein_identity_cluster_id
  )
  frame <- copy(analysis[source == source_name])
  lookup <- function(map, keys) {
    values <- unname(map[keys])
    available <- !is.na(values)
    values[!available] <- global_median
    list(values = values, available = available)
  }
  molecule <- lookup(molecule_medians, frame$canonical_smiles)
  protein <- lookup(protein_medians, frame$protein_sha256)
  cluster <- lookup(cluster_medians, frame$protein_identity_cluster_id)
  frame[, molecule_prior := molecule$values]
  frame[, protein_prior := protein$values]
  frame[, cluster_prior := cluster$values]
  frame[, molecule_prior_available := molecule$available]
  frame[, protein_prior_available := protein$available]
  frame[, cluster_prior_available := cluster$available]
  frame[, combined_entity_prior := fifelse(
    molecule_prior_available & cluster_prior_available,
    (molecule_prior + cluster_prior) / 2,
    fifelse(
      molecule_prior_available,
      molecule_prior,
      fifelse(cluster_prior_available, cluster_prior, global_median)
    )
  )]
  rbindlist(lapply(strata, function(stratum_name) {
    subset <- subset_stratum(frame, stratum_name)
    if (nrow(subset) < 20L) return(NULL)
    data.table(
      source = source_name,
      stratum = stratum_name,
      n = nrow(subset),
      global_median_mae = mean(abs(subset$pkd - global_median)),
      exact_molecule_prior_mae = mean(abs(subset$pkd - subset$molecule_prior)),
      exact_protein_prior_mae = mean(abs(subset$pkd - subset$protein_prior)),
      protein_cluster_prior_mae = mean(abs(subset$pkd - subset$cluster_prior)),
      combined_entity_prior_mae = mean(abs(subset$pkd - subset$combined_entity_prior)),
      exact_molecule_prior_coverage = mean(subset$molecule_prior_available),
      exact_protein_prior_coverage = mean(subset$protein_prior_available),
      protein_cluster_prior_coverage = mean(subset$cluster_prior_available),
      note = "Exploratory development-only priors; unavailable entities fall back to the frozen global median"
    )
  }))
}

entity_baselines <- rbindlist(list(
  entity_baseline_for_source("ChEMBL37", paths$chembl_clusters),
  entity_baseline_for_source("BioLiP/MOAD", paths$biolip_clusters)
), use.names = TRUE, fill = TRUE)
fwrite(entity_baselines, file.path(derived_dir, "V1_7_ENTITY_PRIOR_BASELINES.tsv"), sep = "\t")

v15_metrics <- fread(paths$v15_metrics)
v15_permutations <- fread(paths$v15_permutations)
null_summary <- v15_permutations[, .(
  null_replicates = .N,
  null_spearman_mean = mean(spearman),
  null_spearman_q025 = as.numeric(quantile(spearman, 0.025)),
  null_spearman_q975 = as.numeric(quantile(spearman, 0.975)),
  null_mae_mean = mean(mae),
  null_mae_q025 = as.numeric(quantile(mae, 0.025)),
  null_mae_q975 = as.numeric(quantile(mae, 0.975))
), by = .(stratum, control)]
null_summary <- merge(
  null_summary,
  v15_metrics[, .(stratum, observed_spearman = spearman, observed_mae = mae)],
  by = "stratum",
  all.x = TRUE
)
null_empirical <- v15_permutations[, .(
  spearman_empirical_p_greater = {
    observed <- v15_metrics[stratum == .BY$stratum, spearman][[1L]]
    (1 + sum(spearman >= observed)) / (.N + 1)
  },
  mae_empirical_p_lower = {
    observed <- v15_metrics[stratum == .BY$stratum, mae][[1L]]
    (1 + sum(mae <= observed)) / (.N + 1)
  }
), by = .(stratum, control)]
null_summary <- merge(null_summary, null_empirical, by = c("stratum", "control"))
null_summary[, aligned_feature_necessity_supported := spearman_empirical_p_greater < 0.05]
null_summary[, interpretation := fifelse(
  aligned_feature_necessity_supported,
  "Observed rank exceeds this identity-permutation null in this stratum",
  "Observed rank does not exceed this identity-permutation null in this stratum"
)]
fwrite(null_summary, file.path(derived_dir, "V1_7_PERMUTATION_NULL_HETEROGENEITY.tsv"), sep = "\t")

uncertainty_frame <- copy(analysis)
uncertainty_frame[, prediction := prediction__phase15_esm2_aft_frozen]
uncertainty_frame[, absolute_error := abs(pkd - prediction)]
uncertainty_frame[, interval_covered := as.numeric(absolute_error <= phase15_interval_half_width_pkd)]
uncertainty_frame[, prediction_decile := cut(
  prediction,
  breaks = unique(quantile(prediction, seq(0, 1, 0.1), na.rm = TRUE)),
  include.lowest = TRUE
), by = source]
uncertainty_frame[, document_size_quartile := cut(
  document_component_size,
  breaks = unique(quantile(document_component_size, seq(0, 1, 0.25), na.rm = TRUE)),
  include.lowest = TRUE
), by = source]

summarize_conditional_coverage <- function(frame, group_column, subgroup_type_value) {
  result <- frame[, .(
    n = .N,
    coverage_90 = mean(interval_covered),
    empirical_q90_half_width = as.numeric(quantile(absolute_error, 0.90)),
    fixed_half_width = mean(phase15_interval_half_width_pkd)
  ), by = .(source, subgroup = as.character(get(group_column)))]
  result[, subgroup_type := subgroup_type_value]
  setcolorder(result, c("source", "subgroup_type", "subgroup", setdiff(
    names(result), c("source", "subgroup_type", "subgroup")
  )))
  result
}

conditional_coverage <- rbindlist(list(
  summarize_conditional_coverage(uncertainty_frame, "external_cold_category", "cold_category"),
  summarize_conditional_coverage(uncertainty_frame, "affinity_regime", "affinity_regime"),
  summarize_conditional_coverage(uncertainty_frame, "prediction_decile", "prediction_decile"),
  summarize_conditional_coverage(
    uncertainty_frame, "document_size_quartile", "document_size_quartile"
  )
), use.names = TRUE, fill = TRUE)
conditional_coverage[, half_width_inflation := empirical_q90_half_width / fixed_half_width]
fwrite(
  conditional_coverage,
  file.path(derived_dir, "V1_7_CONDITIONAL_COVERAGE_ATLAS.tsv"),
  sep = "\t"
)

conformal_quantile <- function(errors, level = 0.90) {
  probability <- min(1, ceiling((length(errors) + 1) * level) / length(errors))
  as.numeric(quantile(errors, probability, type = 1, names = FALSE))
}

label_budget_curve <- function(frame, budgets, replicates = 100L, seed = 1L) {
  set.seed(seed)
  clusters <- split(seq_len(nrow(frame)), frame$document_component_id)
  cluster_names <- names(clusters)
  cluster_sizes <- lengths(clusters)
  rbindlist(lapply(budgets, function(budget) {
    rbindlist(lapply(seq_len(replicates), function(replicate) {
      order <- sample(seq_along(cluster_names))
      cumulative <- cumsum(cluster_sizes[order])
      selected_count <- which(cumulative >= budget)[1L]
      calibration_clusters <- cluster_names[order[seq_len(selected_count)]]
      calibration_rows <- unlist(clusters[calibration_clusters], use.names = FALSE)
      test_rows <- setdiff(seq_len(nrow(frame)), calibration_rows)
      if (length(test_rows) < 20L) return(NULL)
      half_width <- conformal_quantile(frame$absolute_error[calibration_rows], 0.90)
      data.table(
        requested_label_budget = budget,
        calibration_n = length(calibration_rows),
        calibration_document_components = length(calibration_clusters),
        test_n = length(test_rows),
        half_width = half_width,
        test_coverage = mean(frame$absolute_error[test_rows] <= half_width),
        replicate = replicate
      )
    }))
  }))
}

budget_replicates <- rbindlist(lapply(seq_along(source_names), function(index) {
  source_name <- source_names[[index]]
  frame <- uncertainty_frame[source == source_name]
  budgets <- c(50L, 100L, 250L, 500L, 1000L)
  if (nrow(frame) >= 5000L) budgets <- c(budgets, 2500L, 5000L)
  result <- label_budget_curve(
    frame,
    budgets = budgets[budgets < 0.6 * nrow(frame)],
    replicates = 100L,
    seed = 20260802L + index
  )
  result[, source := source_name]
  result
}), use.names = TRUE, fill = TRUE)
budget_summary <- budget_replicates[, .(
  repeats = .N,
  median_calibration_n = median(calibration_n),
  median_calibration_document_components = median(calibration_document_components),
  median_half_width = median(half_width),
  half_width_q025 = as.numeric(quantile(half_width, 0.025)),
  half_width_q975 = as.numeric(quantile(half_width, 0.975)),
  median_test_coverage = median(test_coverage),
  coverage_q025 = as.numeric(quantile(test_coverage, 0.025)),
  coverage_q975 = as.numeric(quantile(test_coverage, 0.975))
), by = .(source, requested_label_budget)]
budget_summary[, analysis_role := "Retrospective document-disjoint diagnostic repair curve; not confirmatory"]
fwrite(
  budget_replicates,
  file.path(derived_dir, "V1_7_LOCAL_RECALIBRATION_BUDGET_REPLICATES.tsv"),
  sep = "\t"
)
fwrite(
  budget_summary,
  file.path(derived_dir, "V1_7_LOCAL_RECALIBRATION_BUDGET_SUMMARY.tsv"),
  sep = "\t"
)

biolip_risk <- merge(
  biolip_meta,
  analysis[source == "BioLiP/MOAD", .(
    pair_id,
    document_component_size,
    scaffold_size,
    protein_cluster_size
  )],
  by = "pair_id",
  all.x = TRUE
)
biolip_risk[, extreme_outcome := pkd_outcome < 2 | pkd_outcome > 12]
biolip_risk[, large_label_range := fifelse(is.na(pkd_range), FALSE, pkd_range > 0.30)]
biolip_risk[, multiple_candidates := fifelse(
  is.na(contributing_candidate_count), FALSE, contributing_candidate_count > 1
)]
document_threshold <- as.numeric(quantile(biolip_risk$document_component_size, 0.95, na.rm = TRUE))
biolip_risk[, high_document_concentration := document_component_size >= document_threshold]
biolip_risk[, risk_score :=
  4L * as.integer(extreme_outcome) +
  3L * as.integer(large_label_range) +
  2L * as.integer(multiple_candidates) +
  1L * as.integer(high_document_concentration)
]
setorder(biolip_risk, -risk_score, pair_id)
risk_queue <- biolip_risk[seq_len(min(100L, .N)), .(
  review_order = seq_len(.N),
  pair_id,
  source_class,
  candidate_ids,
  pdb_ids,
  pubmed_ids,
  canonical_smiles,
  protein_sha256,
  pkd_reported = pkd_outcome,
  pkd_min,
  pkd_max,
  pkd_range,
  exact_kd_measurement_count,
  contributing_candidate_count,
  document_component_size,
  risk_score,
  risk_reasons = paste(
    ifelse(extreme_outcome, "extreme_outcome", ""),
    ifelse(large_label_range, "label_range_gt_0.30", ""),
    ifelse(multiple_candidates, "multiple_candidates", ""),
    ifelse(high_document_concentration, "large_document_component", ""),
    sep = ";"
  ),
  reviewer_identity = "",
  source_opened_at = "",
  protein_identity_confirmed = "",
  ligand_identity_confirmed = "",
  exact_kd_confirmed = "",
  value_and_unit_confirmed = "",
  qualifier_exclusion_confirmed = "",
  primary_citation_confirmed = "",
  final_decision = "",
  reviewer_notes = ""
)]
fwrite(
  risk_queue,
  file.path(private_verification_dir, "V1_7_RISK_STRATIFIED_SOURCE_REVIEW_QUEUE_PRIVATE.tsv"),
  sep = "\t"
)

model_panel <- data.table(
  model_id = c("development_median", model_ids),
  model_label = c("Frozen development median", unname(model_labels[model_ids])),
  algorithm_family = c("Constant reference", unname(model_families[model_ids])),
  representation = c(
    "None",
    "Morgan projection and molecular descriptors",
    "Hashed protein trigrams, composition, and length",
    "Concatenated ligand and protein features",
    "Additive features plus projected interaction terms",
    "Morgan projection, ESM-2 PCA, and auxiliary protein features"
  ),
  training_endpoint = c(
    "Median of exact development pKd",
    rep("Exact development pKd", 4),
    "Qualifier-aware exact and censored development labels"
  ),
  feature_dimension_or_state = c(
    "0",
    "72",
    "86",
    "158",
    "287",
    "Frozen serialized pipeline; 64 ESM-2 principal components plus ligand and protein auxiliaries"
  ),
  loss_or_objective = c(
    "None",
    rep("absolute_error", 4),
    "XGBoost survival:aft with normal AFT loss; exact-label Huber post-calibration"
  ),
  fixed_hyperparameters = c(
    "Development exact-label median = 6.174574 pKd",
    rep(
      "learning_rate=0.05; max_iter=220; max_leaf_nodes=31; min_samples_leaf=30; l2_regularization=1.0; early_stopping=FALSE; random_state=42",
      4
    ),
    "tree_method=hist; learning_rate=0.03; max_depth=4; min_child_weight=5; subsample=0.80; colsample_bytree=0.80; lambda=1.0; alpha=0.05; AFT scale=1.0; best_iteration=1199"
  ),
  selection_or_tuning_scope = c(
    "Frozen development statistic",
    rep("No hyperparameter search; specification fixed after v1.5 unblinding and before v1.6 execution", 4),
    "Development-only inner-role selection; model and calibration artifacts frozen before external outcomes"
  ),
  external_outcome_tuning = FALSE,
  analysis_role = c(
    "Simple reference",
    rep("Post-unblinding retrospective specification", 4),
    "Frozen one-time prediction specification"
  )
)
fwrite(model_panel, file.path(table_dir, "Table_S1_MODEL_SPECIFICATIONS.tsv"), sep = "\t")

table_1 <- copy(source_profile)
table_1[, exact_pair_cross_source_overlap := overlap[
  entity == "exact_pair", overlap_n
][[1L]]]
table_1[, analytic_role := fifelse(
  source == "BioLiP/MOAD",
  "Restricted one-time blind source plus retrospective specification analysis",
  "Document- and pair-disjoint external source plus retrospective specification analysis"
)]
table_1[, rights_status := fifelse(
  source == "BioLiP/MOAD",
  "Row-level source restricted; aggregate release only",
  "Public database identifiers and row-level prediction release permitted subject to source terms"
)]
fwrite(table_1, file.path(table_dir, "Table_1_EXTERNAL_SOURCE_DESIGN.tsv"), sep = "\t")
fwrite(utility, file.path(table_dir, "Table_S2_COMPLETE_METRICS_AND_INTERVALS.tsv"), sep = "\t")
fwrite(influence, file.path(table_dir, "Table_S3_DEPENDENCE_AND_INFLUENCE.tsv"), sep = "\t")
fwrite(null_summary, file.path(table_dir, "Table_S4_COMPLETE_NULL_CONTROLS.tsv"), sep = "\t")

lineage <- data.table(
  stage = c(
    "Development labels", "External source 1", "External source 2",
    "v1.5 frozen evaluation", "v1.6 retrospective model panel",
    "v1.7 dependence and decomposition audit", "Primary-source human review"
  ),
  source = c(
    "BindingDB frozen 2024 development snapshot",
    "ChEMBL37",
    "Restricted BioLiP/MOAD exact-Kd reconstruction",
    "Frozen ESM2-AFT predictions",
    "Four HGB specifications plus frozen ESM2-AFT",
    "Frozen v1.6 pair predictions plus source metadata",
    "Frozen random 100-pair queue plus risk-stratified 100-pair queue"
  ),
  record_count_or_unit = c(
    "2,524 rows; 2,364 exact pKd labels",
    "25,436 exact-Kd pairs",
    "1,382 exact-Kd pairs",
    "1,382 restricted-source frozen predictions",
    "26,818 external pairs across five nonreference specifications",
    "26,818 external pairs with hierarchical source metadata",
    "100 random plus 100 risk-stratified review assignments"
  ),
  outcome_access_timing = c(
    "Available for model fitting",
    "Blocked until external prediction freeze",
    "Blocked until v1.5 prediction freeze",
    "One-time blind",
    "Outcomes already open; retrospective",
    "Outcomes already open; retrospective",
    "Pending independent signed adjudication"
  ),
  exclusion_or_aggregation_rule = c(
    "Qualifier-aware; exact labels for HGB and censored labels retained for AFT",
    "Kd transformed to pKd; record aggregation and exact-pair/document separation frozen upstream",
    "Explicit exact Kd only; Ki, IC50, Ka, inequalities, and ambiguous values excluded",
    "No recalibration or post-prediction row deletion",
    "No external-outcome model tuning or post-hoc model omission",
    "No reclassification of v1.5 blind status",
    "Prediction-blind identity, endpoint, value, unit, qualifier, and primary-citation review"
  ),
  public_release_scope = c(
    "Source citation, transformation code, and aggregate development summaries",
    "Row-level ChEMBL37 prediction ledger plus aggregate analyses",
    "Aggregate analyses only; combined row-level reconstruction excluded",
    "Aggregate metrics and immutable hashes; restricted rows excluded",
    "Code, model specifications, aggregate metrics, and public-source rows",
    "Complete aggregate result and figure-source tables",
    "Private until reviewer identities and source-text rights are resolved"
  ),
  confirmatory_status = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
)
fwrite(lineage, file.path(table_dir, "Table_S6_DATA_LINEAGE.tsv"), sep = "\t")

status <- list(
  schema_version = "science_advances_v1_7_scientific_upgrade_v1",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  immutable_prior_result = "v1.5 restricted source shift transport failed",
  retrospective_status = TRUE,
  model_specifications = 5,
  nonreference_algorithm_families = 2,
  external_sources = uniqueN(analysis$source),
  external_pairs = nrow(analysis),
  document_cluster_bootstrap_replicates = 1000,
  local_recalibration_repeats_per_budget = 100,
  independent_random_sample_human_adjudications_signed = 0,
  risk_stratified_human_adjudications_signed = 0,
  prospective_new_source_completed = FALSE,
  interpretation = paste(
    "Source-wide rank association is evaluated jointly with within-document and within-target",
    "rank decomposition, crossed dependence, cluster influence, point utility, and uncertainty transport."
  )
)
write_json(
  status,
  file.path(derived_dir, "V1_7_ANALYSIS_STATUS.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("v1.7 scientific upgrade tables written to: ", package_root)
