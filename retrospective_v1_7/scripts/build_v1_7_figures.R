#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(svglite)
  library(ragg)
})

options(stringsAsFactors = FALSE)
set.seed(20260802)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
script_dir <- dirname(script_path)
package_root <- dirname(script_dir)
derived_dir <- file.path(package_root, "derived")
source_dir <- file.path(package_root, "source_data")
main_dir <- file.path(package_root, "figures", "main")
supp_dir <- file.path(package_root, "figures", "supplementary")
figure_source_dir <- file.path(package_root, "figure_source_data")
qc_dir <- file.path(package_root, "qc")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  profile = file.path(derived_dir, "V1_7_SOURCE_DATA_PROFILE.tsv"),
  overlap = file.path(derived_dir, "V1_7_CROSS_SOURCE_OVERLAP.tsv"),
  rank = file.path(derived_dir, "V1_7_RANK_SIGNAL_DECOMPOSITION.tsv"),
  rank_bootstrap = file.path(derived_dir, "V1_7_OVERALL_DOCUMENT_BOOTSTRAP_RANK.tsv"),
  utility = file.path(derived_dir, "V1_7_UTILITY_CALIBRATION_INFERENCE.tsv"),
  influence = file.path(derived_dir, "V1_7_DOMINANT_CLUSTER_INFLUENCE.tsv"),
  poststratification = file.path(derived_dir, "V1_7_SOURCE_POSTSTRATIFICATION_SENSITIVITY.tsv"),
  entity_priors = file.path(derived_dir, "V1_7_ENTITY_PRIOR_BASELINES.tsv"),
  coverage = file.path(derived_dir, "V1_7_CONDITIONAL_COVERAGE_ATLAS.tsv"),
  budget_summary = file.path(derived_dir, "V1_7_LOCAL_RECALIBRATION_BUDGET_SUMMARY.tsv"),
  budget_replicates = file.path(derived_dir, "V1_7_LOCAL_RECALIBRATION_BUDGET_REPLICATES.tsv"),
  null_summary = file.path(derived_dir, "V1_7_PERMUTATION_NULL_HETEROGENEITY.tsv"),
  predictions = file.path(source_dir, "v1_6", "V1_6_PAIR_LEVEL_MODEL_PREDICTIONS.tsv"),
  agreement = file.path(source_dir, "v1_6", "V1_6_MODEL_AGREEMENT.tsv"),
  v15_null = file.path(source_dir, "v1_5", "V1_5_PERMUTATION_NULL_METRICS.tsv"),
  v15_metrics = file.path(source_dir, "v1_5", "V1_5_STRATUM_METRICS.tsv"),
  assay_forest = file.path(source_dir, "v1_4", "V1_4_MAIN_PAIRED_FOREST.tsv"),
  descriptor_support = file.path(source_dir, "v1_4", "V1_4_MAIN_DESCRIPTOR_SUPPORT.tsv"),
  geometry_contrast = file.path(source_dir, "v0_9", "PHASE26_OBSERVED_VS_PERMUTED_V0_9.tsv"),
  geometry_controls = file.path(source_dir, "v0_9", "PHASE26_NULL_CONTROL_PAIRED_EFFECTS_V0_9.tsv")
)
required_paths <- paths[names(paths) != "predictions"]
missing_paths <- names(required_paths)[!file.exists(unlist(required_paths))]
if (length(missing_paths) > 0L) {
  stop("Missing figure inputs: ", paste(missing_paths, collapse = ", "))
}
pair_level_available <- file.exists(paths$predictions)

profile <- fread(paths$profile)
overlap <- fread(paths$overlap)
rank_data <- fread(paths$rank)
rank_bootstrap <- fread(paths$rank_bootstrap)
utility <- fread(paths$utility)
influence <- fread(paths$influence)
poststratification <- fread(paths$poststratification)
entity_priors <- fread(paths$entity_priors)
coverage <- fread(paths$coverage)
budget_summary <- fread(paths$budget_summary)
budget_replicates <- fread(paths$budget_replicates)
null_summary <- fread(paths$null_summary)
predictions <- if (pair_level_available) fread(paths$predictions) else NULL
agreement <- fread(paths$agreement)
v15_null <- fread(paths$v15_null)
v15_metrics <- fread(paths$v15_metrics)
assay_forest <- fread(paths$assay_forest)
descriptor_support <- fread(paths$descriptor_support)
geometry_contrast <- fread(paths$geometry_contrast)
geometry_controls <- fread(paths$geometry_controls)

model_order <- c(
  "hgb_molecule",
  "hgb_protein",
  "hgb_fusion_additive",
  "hgb_fusion_interaction",
  "phase15_esm2_aft_frozen"
)
model_labels <- c(
  hgb_molecule = "HGB ligand",
  hgb_protein = "HGB protein",
  hgb_fusion_additive = "HGB additive",
  hgb_fusion_interaction = "HGB interaction",
  phase15_esm2_aft_frozen = "Frozen ESM2-AFT"
)
model_colors <- c(
  hgb_molecule = "#1F4E79",
  hgb_protein = "#3B82A0",
  hgb_fusion_additive = "#7A5195",
  hgb_fusion_interaction = "#A05195",
  phase15_esm2_aft_frozen = "#3F3F46"
)
source_order <- c("ChEMBL37", "BioLiP/MOAD")
source_colors <- c("ChEMBL37" = "#D55E00", "BioLiP/MOAD" = "#009E73")
stratum_order <- c(
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
stratum_colors <- c(
  overall = "#4B5563",
  drug_cold = "#0072B2",
  protein_identity_cold = "#CC79A7",
  strict_scaffold_protein_identity_double_cold = "#6C5B9B"
)

label_models <- function(values) factor(
  unname(model_labels[values]), levels = rev(unname(model_labels[model_order]))
)
label_sources <- function(values) factor(values, levels = source_order)
label_strata <- function(values) factor(
  unname(stratum_labels[values]), levels = unname(stratum_labels[stratum_order])
)

theme_pub <- function(base_size = 7.4) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#374151"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#374151"),
      axis.text = element_text(size = 6.7, colour = "#252525"),
      axis.title = element_text(size = 7.1, colour = "#111827"),
      legend.title = element_text(size = 6.8, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key.height = unit(3.5, "mm"),
      strip.background = element_rect(fill = "#F3F4F6", colour = NA),
      strip.text = element_text(size = 6.9, face = "bold", colour = "#111827", margin = margin(1.5, 1.5, 1.5, 1.5)),
      plot.title = element_text(size = 8.0, face = "bold", colour = "#111827", margin = margin(0, 0, 1, 0)),
      plot.subtitle = element_text(size = 6.7, colour = "#4B5563", margin = margin(0, 0, 2, 0)),
      plot.tag = element_text(size = 9.5, face = "bold", colour = "#111827"),
      plot.tag.position = c(0.01, 0.99),
      panel.grid.major.x = element_line(linewidth = 0.25, colour = "#E5E7EB"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 3, 2, 3)
    )
}

theme_set(theme_pub())

save_figure <- function(plot, stem, directory, width_mm = 180, height_mm = 165, tiff = TRUE) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svg_path <- file.path(directory, paste0(stem, ".svg"))
  pdf_path <- file.path(directory, paste0(stem, ".pdf"))
  png_path <- file.path(directory, paste0(stem, ".png"))
  svglite(svg_path, width = width_in, height = height_in, bg = "white")
  print(plot)
  dev.off()
  cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Helvetica", bg = "white")
  print(plot)
  dev.off()
  agg_png(png_path, width = width_in, height = height_in, units = "in", res = 450, background = "white")
  print(plot)
  dev.off()
  if (tiff) {
    tiff_path <- file.path(directory, paste0(stem, ".tiff"))
    agg_tiff(
      tiff_path, width = width_in, height = height_in, units = "in", res = 600,
      background = "white", compression = "lzw"
    )
    print(plot)
    dev.off()
  }
}

write_source <- function(data, name) {
  fwrite(as.data.table(data), file.path(figure_source_dir, paste0(name, ".tsv")), sep = "\t")
}

# Figure 1: scientific question and hierarchical evidence base.
workflow <- data.table(
  x = c(0.52, 1.77, 3.03, 4.28),
  stage = sprintf("%02d", 1:4),
  kicker = c("MODEL DATA", "EXTERNAL DATA", "COMPARISON UNIT", "VALIDATION CONTROL"),
  title = c(
    "Training source",
    "External sources",
    "Hierarchical\nestimands",
    "Temporal\nvalidation"
  ),
  detail = c(
    "BindingDB\n2,364 exact\ntraining labels",
    "ChEMBL37 +\nBioLiP/MOAD\n26,818 external pairs",
    "Pair -> ligand ->\ndocument/scaffold ->\nexact-protein context",
    "Later source\noutcome-blind;\nlocked specifications"
  ),
  fill = c("#E8F1F8", "#FCE9E2", "#E7F5EF", "#F3EAF4"),
  accent = c("#4C78A8", "#D26A3A", "#4F9D86", "#8A6BA8")
)
p1a <- ggplot(workflow, aes(x, 1)) +
  geom_tile(aes(fill = fill), width = 1.00, height = 0.82, colour = "#94A3B8", linewidth = 0.42) +
  geom_rect(
    aes(xmin = x - 0.50, xmax = x + 0.50, ymin = 1.31, ymax = 1.40, fill = accent),
    inherit.aes = FALSE, colour = NA
  ) +
  geom_text(
    aes(x = x - 0.42, y = 1.355, label = stage),
    inherit.aes = FALSE, hjust = 0, size = 2.18, fontface = "bold", colour = "white"
  ) +
  geom_text(
    aes(label = kicker), y = 1.235, size = 1.78, fontface = "bold", colour = "#52606D"
  ) +
  geom_text(aes(label = title), y = 1.075, size = 2.55, lineheight = 0.92, fontface = "bold", colour = "#111827") +
  geom_text(aes(label = detail), y = 0.805, size = 2.12, lineheight = 0.95, colour = "#374151") +
  geom_segment(
    data = data.table(x = c(1.07, 2.32, 3.57), xend = c(1.22, 2.47, 3.72)),
    aes(x = x, xend = xend, y = 1.0, yend = 1.0),
    inherit.aes = FALSE, linewidth = 0.45, colour = "#4B5563",
    arrow = arrow(length = unit(1.8, "mm"), type = "closed")
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.24, 4.56), ylim = c(0.46, 1.56), clip = "off") +
  labs(
    title = "The estimand changes with the comparison unit",
    subtitle = "Data source -> comparison unit -> outcome-blind validation"
  ) +
  theme_void(base_family = "Helvetica", base_size = 7.4) +
  theme(
    plot.title = element_text(size = 8, face = "bold", colour = "#111827"),
    plot.subtitle = element_text(size = 6.7, colour = "#4B5563"),
    plot.margin = margin(3, 5, 2, 5)
  )

profile_long <- melt(
  profile,
  id.vars = "source",
  measure.vars = c(
    "n", "unique_ligands", "unique_proteins", "document_components",
    "effective_document_components"
  ),
  variable.name = "measure",
  value.name = "value"
)
measure_labels <- c(
  n = "Pairs",
  unique_ligands = "Ligands",
  unique_proteins = "Proteins",
  document_components = "Document components",
  effective_document_components = "Effective documents"
)
profile_long[, measure_label := factor(
  unname(measure_labels[measure]), levels = rev(unname(measure_labels))
)]
profile_long[, source := label_sources(source)]
p1b <- ggplot(profile_long, aes(value, measure_label, colour = source)) +
  geom_segment(
    aes(x = 1, xend = value, yend = measure_label), linewidth = 0.55, alpha = 0.55,
    position = position_dodge(width = 0.38)
  ) +
  geom_point(size = 2.3, position = position_dodge(width = 0.38)) +
  geom_text(
    aes(label = ifelse(
      measure == "effective_document_components",
      sprintf("%.1f", value), comma(round(value))
    )),
    hjust = -0.12, size = 2.28, show.legend = FALSE,
    position = position_dodge(width = 0.38)
  ) +
  scale_x_log10(
    labels = function(values) comma(values, accuracy = 1),
    expand = expansion(mult = c(0.02, 0.35))
  ) +
  scale_colour_manual(values = source_colors, drop = FALSE) +
  labs(
    title = "Source support",
    subtitle = "Nominal vs concentration-equivalent counts",
    x = "Count (log scale)", y = NULL, colour = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

if (pair_level_available) {
  predictions[, source := label_sources(source)]
  affinity_density <- rbindlist(lapply(source_order, function(source_name) {
    observed <- predictions[source == source_name, pkd]
    estimate <- density(observed, from = 0, to = 14, n = 512, adjust = 1.05)
    data.table(source = source_name, pkd = estimate$x, density = estimate$y)
  }))
} else {
  density_path <- file.path(figure_source_dir, "Figure_1C_affinity_density.tsv")
  if (!file.exists(density_path)) stop("Missing public aggregate input: ", density_path)
  affinity_density <- fread(density_path)
}
affinity_density[, source := label_sources(as.character(source))]
affinity_guides <- data.table(
  pkd = c(2, 12),
  density = max(affinity_density$density) * 0.88,
  label = c("pKd 2", "pKd 12")
)
p1c <- ggplot(affinity_density, aes(pkd, density, fill = source, colour = source)) +
  geom_area(alpha = 0.20, position = "identity") +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = c(2, 12), linetype = "dotted", linewidth = 0.35, colour = "#6B7280") +
  geom_text(
    data = affinity_guides,
    aes(pkd, density, label = label),
    inherit.aes = FALSE, vjust = -0.35, size = 2.15, colour = "#4B5563"
  ) +
  scale_fill_manual(values = source_colors, drop = FALSE) +
  scale_colour_manual(values = source_colors, drop = FALSE) +
  labs(
    title = "Affinity shifts by source",
    subtitle = "Dotted guides: observed pKd 2 and 12",
    x = "Observed pKd", y = "Density", fill = "Source", colour = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "none")

repeat_data <- profile[, .(
  source,
  Ligand = repeated_ligand_row_fraction,
  Protein = repeated_protein_row_fraction,
  Document = repeated_document_row_fraction
)]
repeat_data <- melt(repeat_data, id.vars = "source", variable.name = "unit", value.name = "fraction")
repeat_data[, source := label_sources(source)]
p1d <- ggplot(repeat_data, aes(unit, fraction, fill = source)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.64, colour = "white", linewidth = 0.2) +
  geom_text(
    aes(label = percent(fraction, accuracy = 1)),
    position = position_dodge(width = 0.74), vjust = -0.35, size = 2.25
  ) +
  scale_fill_manual(values = source_colors, drop = FALSE) +
  scale_y_continuous(labels = percent, limits = c(0, 1.08), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Repeated entities are common",
    x = NULL, y = "Rows in repeated groups", fill = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

figure1 <- p1a / (p1b | p1c | p1d) +
  plot_layout(heights = c(0.76, 1.24), widths = c(1.22, 1, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(figure1, "Figure_1_v1_7_hierarchical_evidence_base", main_dir, 180, 148)
write_source(workflow, "Figure_1A_workflow")
write_source(profile_long, "Figure_1B_source_support")
write_source(affinity_density, "Figure_1C_affinity_density")
write_source(repeat_data, "Figure_1D_repeated_units")

# Figure 2: specification-consistent source-wide rank and inconsistent point accuracy.
rank_bootstrap[, model_label_plot := label_models(model_id)]
rank_bootstrap[, source := label_sources(source)]
p2a <- ggplot(rank_bootstrap, aes(estimate, model_label_plot, colour = source)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#6B7280") +
  geom_errorbarh(
    aes(xmin = simultaneous_ci95_low, xmax = simultaneous_ci95_high),
    height = 0, linewidth = 0.75, position = position_dodge(width = 0.45)
  ) +
  geom_point(size = 2.25, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = source_colors, drop = FALSE) +
  labs(
    title = "Source-wide rank survives joint inference",
    subtitle = "1,000 component bootstraps; max-t 95% CI",
    x = "Spearman rho", y = NULL, colour = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

rank_heat <- copy(rank_data)
rank_heat[, model_label_plot := factor(
  unname(model_labels[model_id]), levels = rev(unname(model_labels[model_order]))
)]
rank_heat[, stratum_plot := factor(
  unname(stratum_labels[stratum]), levels = unname(stratum_labels[stratum_order])
)]
rank_heat[, source := label_sources(source)]
p2b <- ggplot(rank_heat, aes(stratum_plot, model_label_plot, fill = spearman)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", spearman)), size = 2.28, colour = "#111827") +
  facet_wrap(~source, ncol = 1) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
    limits = c(-0.1, 0.65), oob = squish
  ) +
  labs(
    title = "Rank varies by source and shift",
    x = NULL, y = NULL, fill = "Spearman\nrho"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

utility_plot <- utility[stratum %in% c("overall", "protein_identity_cold", "strict_scaffold_protein_identity_double_cold")]
utility_plot[, model_label_plot := label_models(model_id)]
utility_plot[, source := label_sources(source)]
utility_plot[, stratum_plot := factor(
  unname(stratum_labels[stratum]),
  levels = unname(stratum_labels[c("overall", "protein_identity_cold", "strict_scaffold_protein_identity_double_cold")])
)]
p2c <- ggplot(utility_plot, aes(delta_mae, model_label_plot, colour = stratum_plot)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#111827") +
  geom_errorbarh(
    aes(
      xmin = delta_mae_dependence_envelope_ci95_low,
      xmax = delta_mae_dependence_envelope_ci95_high
    ),
    height = 0, linewidth = 0.55,
    position = position_dodge(width = 0.56)
  ) +
  geom_point(size = 1.9, position = position_dodge(width = 0.56)) +
  facet_wrap(~source) +
  scale_colour_manual(values = unname(stratum_colors[c(
    "overall", "protein_identity_cold", "strict_scaffold_protein_identity_double_cold"
  )])) +
  labs(
    title = "Point accuracy is not consistent",
    subtitle = "Negative favors model; dependence-envelope 95% CI",
    x = expression(Delta*"MAE versus frozen median (pKd)"), y = NULL, colour = "Stratum"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

family_nodes <- data.table(
  x = c(0.90, 2.70),
  y = c(5.45, 5.45),
  label = c("HGB retrospective family\nn = 4 related variants", "Frozen comparator\nn = 1"),
  fill = c("#E7EEF6", "#ECECF0"),
  text_colour = c("#1F4E79", "#3F3F46")
)
spec_nodes <- data.table(
  x = c(rep(1.45, 4), 2.70),
  y = c(4.35, 3.30, 2.25, 1.20, 3.30),
  label = unname(model_labels[model_order]),
  fill = unname(model_colors[model_order])
)
branch_data <- rbind(
  data.table(x = 0.90, xend = 0.90, y = 5.00, yend = 1.20),
  data.table(x = 0.90, xend = 1.45, y = c(4.35, 3.30, 2.25, 1.20), yend = c(4.35, 3.30, 2.25, 1.20)),
  data.table(x = 2.70, xend = 2.70, y = 5.00, yend = 3.30)
)
p2d <- ggplot() +
  geom_segment(
    data = branch_data,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, linewidth = 0.55, colour = "#9CA3AF"
  ) +
  geom_label(
    data = family_nodes,
    aes(x, y, label = label, fill = fill, colour = text_colour),
    inherit.aes = FALSE, fontface = "bold", size = 2.25,
    label.size = 0.25, label.padding = unit(0.18, "lines")
  ) +
  geom_label(
    data = spec_nodes,
    aes(x, y, label = label, fill = fill),
    inherit.aes = FALSE, colour = "white", fontface = "bold", size = 2.25,
    label.size = 0.20, label.padding = unit(0.16, "lines")
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0.42, 3.18), ylim = c(0.62, 6.02), clip = "off") +
  labs(
    title = "Model specifications form two families",
    subtitle = "Shared HGB lineage versus one frozen comparator"
  ) +
  theme_void(base_family = "Helvetica", base_size = 7.4) +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.7, colour = "#4B5563"),
    plot.margin = margin(2, 3, 2, 3)
  )

figure2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(widths = c(1.05, 0.95), heights = c(1.08, 0.92)) +
  plot_annotation(tag_levels = "A")
save_figure(figure2, "Figure_2_v1_7_rank_and_point_inference", main_dir, 180, 178)
write_source(rank_bootstrap, "Figure_2A_overall_rank_bootstrap")
write_source(rank_heat, "Figure_2B_rank_heatmap")
write_source(utility_plot, "Figure_2C_point_accuracy_forest")

# Figure 3: hero decomposition from source-wide to within-context rank.
overall_rank <- rank_data[stratum == "overall"]
overall_rank[, model_label_plot := label_models(model_id)]
overall_rank[, source := label_sources(source)]
attenuation <- melt(
  overall_rank,
  id.vars = c("source", "model_id", "model_label_plot"),
  measure.vars = c("spearman", "document_residual_rank_rho"),
  variable.name = "estimand",
  value.name = "rho"
)
attenuation[, estimand_label := factor(
  c(spearman = "Source wide", document_residual_rank_rho = "Within document")[estimand],
  levels = c("Source wide", "Within document")
)]
p3a <- ggplot() +
  geom_segment(
    data = dcast(attenuation, source + model_id + model_label_plot ~ estimand, value.var = "rho"),
    aes(
      x = spearman, xend = document_residual_rank_rho,
      y = model_label_plot, yend = model_label_plot, colour = source
    ), linewidth = 0.8, alpha = 0.65,
    arrow = arrow(length = unit(1.5, "mm"), type = "closed")
  ) +
  geom_point(
    data = attenuation,
    aes(rho, model_label_plot, shape = estimand_label, fill = source, colour = source),
    size = 2.35, stroke = 0.55
  ) +
  facet_wrap(~source) +
  scale_colour_manual(values = source_colors) +
  scale_y_discrete(limits = rev(unname(model_labels[model_order])), drop = FALSE) +
  scale_fill_manual(values = source_colors) +
  scale_shape_manual(values = c("Source wide" = 21, "Within document" = 24)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#6B7280") +
  labs(
    title = "Source-wide rank attenuates within documents",
    x = "Spearman or residual-rank correlation", y = NULL,
    shape = NULL, fill = "Source", colour = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

between_data <- rbindlist(list(
  overall_rank[, .(
    source, model_id, model_label_plot,
    group = "Document component",
    fraction = document_between_covariance_fraction
  )],
  overall_rank[, .(
    source, model_id, model_label_plot,
    group = "Exact protein",
    fraction = exact_protein_between_covariance_fraction
  )]
))
between_data[, group := factor(group, levels = c("Document component", "Exact protein"))]
p3b <- ggplot(between_data, aes(fraction, model_label_plot, colour = group, shape = group)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.4, colour = "#9CA3AF") +
  geom_point(size = 2.15, position = position_dodge(width = 0.45)) +
  facet_wrap(~source) +
  scale_colour_manual(values = c("Document component" = "#7A5195", "Exact protein" = "#3B82A0")) +
  scale_shape_manual(values = c("Document component" = 16, "Exact protein" = 17)) +
  scale_x_continuous(labels = percent, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    title = "Rank covariance is mostly between contexts",
    subtitle = "Descriptive covariance fraction; not causal attribution",
    x = "Between-group fraction of rank covariance", y = NULL, colour = NULL, shape = NULL
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

within_protein <- overall_rank[model_id != "hgb_protein"]
not_estimable <- unique(overall_rank[model_id == "hgb_protein", .(source)])
not_estimable[, `:=`(
  model_id = "hgb_protein",
  model_label_plot = factor("HGB protein", levels = rev(unname(model_labels[model_order]))),
  exact_protein_macro_spearman = 0,
  label = "not estimable"
)]
p3c <- ggplot(within_protein, aes(exact_protein_macro_spearman, model_label_plot, colour = source)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#6B7280") +
  geom_segment(aes(x = 0, xend = exact_protein_macro_spearman, yend = model_label_plot), linewidth = 0.65, alpha = 0.6) +
  geom_point(size = 2.35) +
  geom_text(
    data = not_estimable,
    aes(x = 0.02, y = model_label_plot, label = label),
    inherit.aes = FALSE, colour = "#6B7280", hjust = 0, size = 2.05
  ) +
  facet_wrap(~source) +
  scale_colour_manual(values = source_colors) +
  scale_y_discrete(limits = rev(unname(model_labels[model_order])), drop = FALSE) +
  scale_x_continuous(limits = c(0, 0.35), breaks = c(0, 0.1, 0.2, 0.3)) +
  labs(
    title = "Within-protein ligand ordering is weak",
    subtitle = "Macro Spearman among exact proteins with at least three pairs",
    x = "Macro within-protein Spearman rho", y = NULL, colour = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

double_rank <- rank_data[stratum == "strict_scaffold_protein_identity_double_cold"]
double_rank[, model_label_plot := label_models(model_id)]
double_rank[, source := label_sources(source)]
double_long <- melt(
  double_rank,
  id.vars = c("source", "model_id", "model_label_plot"),
  measure.vars = c("spearman", "document_residual_rank_rho"),
  variable.name = "estimand", value.name = "rho"
)
double_long[, estimand_label := factor(
  c(spearman = "Source-stratum", document_residual_rank_rho = "Within document")[estimand],
  levels = c("Source-stratum", "Within document")
)]
p3d <- ggplot(double_long, aes(rho, model_label_plot, colour = estimand_label, shape = estimand_label)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  geom_point(size = 2.1, position = position_dodge(width = 0.46)) +
  facet_wrap(~source) +
  scale_colour_manual(values = c("Source-stratum" = "#4B5563", "Within document" = "#7A5195")) +
  scale_shape_manual(values = c("Source-stratum" = 16, "Within document" = 17)) +
  labs(
    title = "Double-cold rank collapses within documents",
    x = "Rank correlation", y = NULL, colour = NULL, shape = NULL
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

figure3 <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(widths = c(1.05, 0.95), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(figure3, "Figure_3_v1_7_rank_signal_decomposition", main_dir, 180, 174)
write_source(attenuation, "Figure_3A_rank_attenuation")
write_source(between_data, "Figure_3B_between_context_covariance")
write_source(within_protein, "Figure_3C_within_protein_rank")
write_source(not_estimable, "Figure_3C_not_estimable")
write_source(double_long, "Figure_3D_double_cold_rank")

# Figure 4: source composition and point-accuracy sensitivity.
influence_overall <- influence[stratum == "overall" & group_type == "document_component"]
influence_overall[, model_label_plot := label_models(model_id)]
influence_overall[, source := label_sources(source)]
influence_long <- melt(
  influence_overall,
  id.vars = c("source", "model_id", "model_label_plot", "removed_fraction"),
  measure.vars = c("full_delta_mae", "reduced_delta_mae"),
  variable.name = "analysis", value.name = "delta_mae"
)
influence_long[, analysis_label := factor(
  c(full_delta_mae = "Full source", reduced_delta_mae = "Largest component removed")[analysis],
  levels = c("Full source", "Largest component removed")
)]
p4a <- ggplot() +
  geom_segment(
    data = influence_overall,
    aes(
      x = full_delta_mae, xend = reduced_delta_mae,
      y = model_label_plot, yend = model_label_plot, colour = source
    ), linewidth = 0.8, alpha = 0.65,
    arrow = arrow(length = unit(1.5, "mm"), type = "closed")
  ) +
  geom_point(
    data = influence_long,
    aes(delta_mae, model_label_plot, shape = analysis_label, fill = source, colour = source),
    size = 2.25, stroke = 0.5
  ) +
  facet_wrap(~source) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  scale_colour_manual(values = source_colors) +
  scale_fill_manual(values = source_colors) +
  scale_shape_manual(values = c("Full source" = 21, "Largest component removed" = 24)) +
  labs(
    title = "Dominant document shifts MAE",
    subtitle = "Largest component: ChEMBL37 22.2%; BioLiP/MOAD 3.0%",
    x = expression(Delta*"MAE versus frozen median (pKd)"), y = NULL,
    shape = NULL, fill = "Source", colour = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank()) +
  guides(colour = "none", fill = "none", shape = guide_legend(title = NULL))

post_plot <- poststratification
post_plot[, model_label_plot := label_models(model_id)]
post_plot[, source := label_sources(source)]
post_plot[, analysis_label := factor(
  c(raw_common_support = "Raw common support", affinity_and_cold_poststratified = "Poststratified")[analysis],
  levels = c("Raw common support", "Poststratified")
)]
p4b <- ggplot(post_plot, aes(delta_mae_vs_development_median, model_label_plot, colour = analysis_label, shape = analysis_label)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  geom_point(size = 2.05, position = position_dodge(width = 0.48)) +
  facet_wrap(~source) +
  scale_colour_manual(values = c("Raw common support" = "#6B7280", "Poststratified" = "#0072B2")) +
  scale_shape_manual(values = c("Raw common support" = 16, "Poststratified" = 17)) +
  labs(
    title = "Poststratification shrinks differences",
    subtitle = "Affinity and cold-category support balanced",
    x = expression(Delta*"MAE versus frozen median (pKd)"), y = NULL, colour = NULL, shape = NULL
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

prior_plot <- rbindlist(list(
  entity_priors[, .(
    source, stratum, baseline = "Exact molecule prior",
    coverage = exact_molecule_prior_coverage, mae = exact_molecule_prior_mae
  )],
  entity_priors[, .(
    source, stratum, baseline = "Protein-cluster prior",
    coverage = protein_cluster_prior_coverage, mae = protein_cluster_prior_mae
  )]
))
prior_plot <- prior_plot[stratum %in% c("overall", "drug_cold", "protein_identity_cold")]
prior_plot[, source := label_sources(source)]
prior_plot[, stratum_plot := label_strata(stratum)]
p4c <- ggplot(prior_plot, aes(coverage, mae, colour = stratum_plot, shape = baseline)) +
  geom_point(size = 2.15, alpha = 0.9) +
  facet_wrap(~source, scales = "free_y") +
  scale_x_continuous(labels = percent, limits = c(-0.03, 1.03)) +
  scale_colour_manual(values = unname(stratum_colors[stratum_order])) +
  scale_shape_manual(values = c("Exact molecule prior" = 16, "Protein-cluster prior" = 17)) +
  labs(
    title = "Entity priors are sparse",
    x = "Development-prior coverage", y = "MAE (pKd)", colour = "Stratum", shape = "Baseline"
  ) +
  theme_pub() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

if (pair_level_available) {
  pred_long <- melt(
    predictions,
    id.vars = c("source", "pair_id", "pkd"),
    measure.vars = paste0("prediction__", model_order),
    variable.name = "prediction_column", value.name = "prediction"
  )
  pred_long[, model_id := sub("^prediction__", "", prediction_column)]
  pred_long[, affinity_regime := cut(
    pkd, breaks = c(-Inf, 5, 6, 7, 8, Inf),
    labels = c("<5", "5-<6", "6-<7", "7-<8", ">=8"), right = FALSE
  )]
  error_regime <- pred_long[, .(mae = mean(abs(pkd - prediction)), n = .N), by = .(source, model_id, affinity_regime)]
} else {
  error_path <- file.path(figure_source_dir, "Figure_4D_affinity_regime_error.tsv")
  if (!file.exists(error_path)) stop("Missing public aggregate input: ", error_path)
  error_regime <- fread(error_path)
}
error_regime[, model_label_plot := factor(
  unname(model_labels[model_id]), levels = rev(unname(model_labels[model_order]))
)]
error_regime[, source := label_sources(as.character(source))]
error_regime[, affinity_regime := factor(
  as.character(affinity_regime), levels = c("<5", "5-<6", "6-<7", "7-<8", ">=8")
)]
regime_counts <- dcast(
  error_regime[, .(n = first(n)), by = .(affinity_regime, source)],
  affinity_regime ~ source,
  value.var = "n"
)
compact_count <- function(values) {
  ifelse(values >= 1000, sprintf("%.1fk", values / 1000), comma(values))
}
regime_count_labels <- setNames(
  paste0(
    as.character(regime_counts$affinity_regime),
    "\n",
    compact_count(regime_counts[["ChEMBL37"]]),
    " / ",
    compact_count(regime_counts[["BioLiP/MOAD"]])
  ),
  as.character(regime_counts$affinity_regime)
)
p4d <- ggplot(error_regime, aes(affinity_regime, model_label_plot, fill = mae)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", mae)), size = 2.12, colour = "#111827") +
  facet_wrap(~source, ncol = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#B2182B") +
  labs(
    title = "Errors concentrate in affinity tails",
    x = "Observed pKd regime (n = ChEMBL37 / BioLiP/MOAD)", y = NULL, fill = "MAE"
  ) +
  scale_x_discrete(labels = regime_count_labels) +
  theme_pub() +
  theme(
    panel.grid = element_blank(), legend.position = "right",
    axis.text.x = element_text(angle = 25, hjust = 1, size = 5.6)
  )

figure4 <- (p4a | p4b) / (p4c | p4d) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(figure4, "Figure_4_v1_7_composition_and_point_accuracy", main_dir, 180, 178)
write_source(influence_long, "Figure_4A_dominant_component")
write_source(post_plot, "Figure_4B_poststratification")
write_source(prior_plot, "Figure_4C_entity_priors")
write_source(error_regime, "Figure_4D_affinity_regime_error")

# Figure 5: predictive uncertainty and local recalibration.
coverage_overall <- utility[stratum == "overall" & model_id == "phase15_esm2_aft_frozen"]
coverage_overall[, source := label_sources(source)]
p5a <- ggplot(coverage_overall, aes(coverage_90, source, colour = source)) +
  geom_vline(xintercept = 0.9, linetype = "dashed", linewidth = 0.55, colour = "#111827") +
  geom_errorbarh(
    aes(
      xmin = coverage_90_dependence_envelope_ci95_low,
      xmax = coverage_90_dependence_envelope_ci95_high
    ), height = 0, linewidth = 0.8
  ) +
  geom_point(size = 2.6) +
  geom_text(aes(label = percent(coverage_90, accuracy = 0.1)), hjust = -0.35, size = 2.35) +
  scale_colour_manual(values = source_colors, guide = "none") +
  scale_x_continuous(labels = percent, limits = c(0.5, 0.94), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9)) +
  labs(
    title = "The frozen nominal 90% interval under-covers",
    subtitle = "Widest dependence-envelope 95% interval; dashed line = 90% target",
    x = "External coverage", y = NULL
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

coverage_affinity <- coverage[subgroup_type == "affinity_regime"]
coverage_affinity[, source := label_sources(source)]
coverage_affinity[, subgroup := factor(subgroup, levels = c("<5", "5-<6", "6-<7", "7-<8", ">=8"))]
p5b <- ggplot(coverage_affinity, aes(subgroup, coverage_90, fill = source)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", linewidth = 0.5, colour = "#111827") +
  geom_col(position = position_dodge(width = 0.74), width = 0.65, colour = "white", linewidth = 0.2) +
  geom_text(
    aes(label = percent(coverage_90, accuracy = 1)),
    position = position_dodge(width = 0.74), vjust = -0.35, size = 2.2
  ) +
  scale_fill_manual(values = source_colors, drop = FALSE) +
  scale_y_continuous(labels = percent, limits = c(0, 1.04), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Miscoverage is concentrated in both affinity tails",
    x = "Observed pKd regime", y = "Coverage", fill = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.x = element_blank())

budget_summary[, source := label_sources(source)]
p5c <- ggplot(budget_summary, aes(requested_label_budget, median_test_coverage, colour = source, fill = source)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", linewidth = 0.5, colour = "#111827") +
  geom_ribbon(aes(ymin = coverage_q025, ymax = coverage_q975), alpha = 0.13, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.9) +
  scale_x_log10(breaks = c(50, 100, 250, 500, 1000, 2500, 5000), labels = comma) +
  scale_y_continuous(labels = percent, limits = c(0.65, 1.0)) +
  scale_colour_manual(values = source_colors) +
  scale_fill_manual(values = source_colors) +
  labs(
    title = "Local labels can restore median coverage",
    subtitle = "Median across 100 document-disjoint splits; ribbons show 95% ranges",
    x = "Requested local calibration labels", y = "Test coverage", colour = "Source", fill = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "top")

p5d <- ggplot(budget_summary, aes(requested_label_budget, median_half_width, colour = source, fill = source)) +
  geom_hline(yintercept = 1.455871, linetype = "dashed", linewidth = 0.5, colour = "#111827") +
  geom_ribbon(aes(ymin = half_width_q025, ymax = half_width_q975), alpha = 0.13, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.9) +
  scale_x_log10(breaks = c(50, 100, 250, 500, 1000, 2500, 5000), labels = comma) +
  scale_colour_manual(values = source_colors) +
  scale_fill_manual(values = source_colors) +
  labs(
    title = "Coverage repair costs interval sharpness",
    subtitle = "Median half-width; dashed line = frozen 1.456-pKd half-width",
    x = "Requested local calibration labels", y = "Recalibrated half-width (pKd)", colour = "Source", fill = "Source"
  ) +
  theme_pub() +
  theme(legend.position = "none")

figure5 <- (p5a | p5b) / (p5c | p5d) +
  plot_layout(widths = c(0.9, 1.1), heights = c(0.92, 1.08)) +
  plot_annotation(tag_levels = "A")
save_figure(figure5, "Figure_5_v1_7_uncertainty_transport_and_repair", main_dir, 180, 166)
write_source(coverage_overall, "Figure_5A_overall_coverage")
write_source(coverage_affinity, "Figure_5B_affinity_coverage")
write_source(budget_summary, "Figure_5CD_recalibration_budget")

# Figure 6: feature-necessity heterogeneity and bounded conclusions.
null_plot <- null_summary[stratum %in% c(
  "all_exact_kd", "drug_cold", "protein_identity_cold", "strict_double_cold"
)]
null_stratum_labels <- c(
  all_exact_kd = "Overall",
  drug_cold = "Ligand cold",
  protein_identity_cold = "Target cold",
  strict_double_cold = "Double cold"
)
null_control_labels <- c(
  protein_permutation = "Protein",
  ligand_permutation = "Ligand",
  joint_permutation = "Joint"
)
null_plot[, stratum_plot := factor(
  unname(null_stratum_labels[stratum]), levels = unname(null_stratum_labels)
)]
null_plot[, control_plot := factor(
  unname(null_control_labels[control]), levels = unname(null_control_labels)
)]
null_plot[, signed_strength := pmin(2.1, -log10(spearman_empirical_p_greater))]
p6a <- ggplot(null_plot, aes(control_plot, stratum_plot, fill = signed_strength)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(
    aes(label = ifelse(spearman_empirical_p_greater < 0.05, sprintf("P=%.3f", spearman_empirical_p_greater), sprintf("P=%.2f", spearman_empirical_p_greater))),
    size = 2.25, colour = "#111827"
  ) +
  scale_fill_gradient(low = "#F3F4F6", high = "#0072B2", limits = c(0, 2.1)) +
  labs(
    title = "Feature necessity varies by shift",
    subtitle = "One-sided finite-null P values for observed Spearman rho",
    x = NULL, y = NULL, fill = "-log10(P)"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

p6b <- ggplot(null_plot, aes(null_spearman_mean, stratum_plot, colour = control_plot)) +
  geom_errorbarh(
    aes(xmin = null_spearman_q025, xmax = null_spearman_q975),
    height = 0, linewidth = 0.7, position = position_dodge(width = 0.48)
  ) +
  geom_point(size = 1.9, position = position_dodge(width = 0.48)) +
  geom_point(
    aes(x = observed_spearman), inherit.aes = TRUE,
    shape = 23, fill = "white", colour = "#111827", size = 2.1,
    position = position_dodge(width = 0.48)
  ) +
  facet_wrap(~control_plot) +
  scale_colour_manual(values = c(
    "Protein" = "#56B4E9",
    "Ligand" = "#E69F00",
    "Joint" = "#6B7280"
  ), guide = "none") +
  labs(
    title = "Observed rank exceeds only some nulls",
    subtitle = "Lines: null 95%; circles: null mean; diamonds: observed",
    x = "Spearman rho", y = NULL
  ) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

assay_specificity <- assay_forest[reference_label == "Joint permutation", .(
  stream = "Assay context",
  comparison = category_label,
  estimate = mae_difference_candidate_minus_reference,
  low = ci_95_low,
  high = ci_95_high,
  n = paired_documents
)]
geometry_specificity <- geometry_contrast[, .(
  stream = "Geometry",
  comparison = challenge_label,
  estimate = observed_minus_permuted_delta_mae,
  low = ci95_lower,
  high = ci95_upper,
  n
)]
specificity <- rbindlist(list(assay_specificity, geometry_specificity), fill = TRUE)
specificity[, comparison_plot := factor(
  paste(stream, comparison, sep = ": "),
  levels = rev(paste(stream, comparison, sep = ": "))
)]
p6c <- ggplot(specificity, aes(estimate, comparison_plot, colour = stream)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#111827") +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = 0.75) +
  geom_point(size = 2.3) +
  scale_colour_manual(values = c("Assay context" = "#009E73", "Geometry" = "#CC79A7")) +
  labs(
    title = "Context and geometry lack specificity",
    subtitle = "Observed minus aligned-feature permutation; 95% CI",
    x = expression(Delta*"MAE contrast (pKd)"), y = NULL, colour = "Stream"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

claim_map <- data.table(
  claim = factor(
    c(
      "Source-wide rank association",
      "Within-document rank",
      "Within-target ligand ordering",
      "Point accuracy vs constant",
      "Nominal 90% coverage",
      "Aligned-feature necessity",
      "Assay / geometry specificity"
    ),
    levels = rev(c(
      "Source-wide rank association",
      "Within-document rank",
      "Within-target ligand ordering",
      "Point accuracy vs constant",
      "Nominal 90% coverage",
      "Aligned-feature necessity",
      "Assay / geometry specificity"
    ))
  ),
  status = factor(
    c("Supported", "Attenuated", "Weak", "Source dependent", "Rejected", "Heterogeneous", "Rejected"),
    levels = c("Supported", "Attenuated", "Weak", "Source dependent", "Heterogeneous", "Rejected")
  ),
  interpretation = c(
    "Supported",
    "Attenuated",
    "Weak",
    "Source dependent",
    "Under-covered",
    "Heterogeneous",
    "Not specific"
  ),
  evidence = c(
    "10/10 joint CIs > 0",
    "rho 0.005-0.338",
    "macro rho 0.074-0.300",
    "composition sensitive",
    "64-66% overall",
    "P varies by stratum",
    "null-specific CIs cross 0"
  )
)
status_colors <- c(
  Supported = "#009E73",
  Attenuated = "#56B4E9",
  Weak = "#E69F00",
  `Source dependent` = "#D55E00",
  Heterogeneous = "#CC79A7",
  Rejected = "#6B7280"
)
p6d <- ggplot(claim_map, aes(1, claim, fill = status)) +
  geom_tile(width = 0.22, height = 0.70, colour = "white") +
  geom_text(
    aes(x = 1.20, label = paste0(interpretation, "; ", evidence)),
    hjust = 0, size = 2.15, colour = "#111827"
  ) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  coord_cartesian(xlim = c(0.82, 3.00), clip = "off") +
  labs(
    title = "Claim strength follows the estimand",
    subtitle = "Next gate: temporal, family, and target matched",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_pub() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
    panel.grid = element_blank(), legend.position = "none"
  )

figure6 <- (p6a | p6b) / (p6c | p6d) +
  plot_layout(widths = c(0.95, 1.05), heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "A")
save_figure(figure6, "Figure_6_v1_7_feature_necessity_and_claim_boundary", main_dir, 180, 174)
write_source(null_plot, "Figure_6AB_permutation_heterogeneity")
write_source(specificity, "Figure_6C_specificity_controls")
write_source(claim_map, "Figure_6D_claim_map")

# Supplementary Figure 1: source concentration and overlap.
overlap_plot <- copy(overlap)
overlap_plot[, entity_label := factor(
  c(
    standardized_ligand = "Exact ligand",
    ligand_scaffold = "Ligand scaffold",
    exact_protein = "Exact protein",
    protein_identity_cluster = "Protein cluster label",
    exact_pair = "Exact pair"
  )[entity],
  levels = rev(c("Exact ligand", "Ligand scaffold", "Exact protein", "Protein cluster label", "Exact pair"))
)]
ps1a <- ggplot(overlap_plot, aes(overlap_n, entity_label)) +
  geom_segment(aes(x = 0, xend = overlap_n, yend = entity_label), linewidth = 0.7, colour = "#9CA3AF") +
  geom_point(size = 2.6, colour = "#0072B2") +
  geom_text(aes(label = overlap_n), hjust = -0.35, size = 2.35) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Cross-source overlap", x = "Shared identifiers", y = NULL) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

concentration <- profile[, .(
  source,
  largest_fraction = largest_document_component_fraction,
  effective = effective_document_components,
  nominal = document_components
)]
concentration[, source := label_sources(source)]
ps1b <- ggplot(concentration, aes(nominal, effective, colour = source, size = largest_fraction)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9CA3AF") +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = source), nudge_y = c(-12, 12), size = 2.35, show.legend = FALSE) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_colour_manual(values = source_colors, guide = "none") +
  scale_size_continuous(labels = percent, range = c(3, 7)) +
  labs(
    title = "Nominal and effective document units diverge",
    x = "Nominal connected components", y = "Effective components", size = "Largest\ncomponent"
  ) +
  theme_pub()

ps1c <- p1c + labs(title = "Complete external outcome distributions") + theme(legend.position = "top")
ps1d <- p1d + labs(title = "Repeated-unit fractions")
supp1 <- (ps1a | ps1b) / (ps1c | ps1d) + plot_annotation(tag_levels = "A")
save_figure(supp1, "Supplementary_Figure_1_v1_7_source_structure", supp_dir, 180, 145, tiff = FALSE)

# Supplementary Figure 2: complete rank intervals under each dependence scheme.
rank_intervals <- rbindlist(list(
  rank_data[, .(source, stratum, model_id, estimate = spearman, low = document_cluster_ci95_low, high = document_cluster_ci95_high, scheme = "Document")],
  rank_data[, .(source, stratum, model_id, estimate = spearman, low = scaffold_cluster_ci95_low, high = scaffold_cluster_ci95_high, scheme = "Scaffold")],
  rank_data[, .(source, stratum, model_id, estimate = spearman, low = protein_cluster_ci95_low, high = protein_cluster_ci95_high, scheme = "Protein cluster")],
  rank_data[, .(source, stratum, model_id, estimate = spearman, low = three_way_cluster_ci95_low, high = three_way_cluster_ci95_high, scheme = "Three way")],
  rank_data[, .(source, stratum, model_id, estimate = spearman, low = dependence_envelope_ci95_low, high = dependence_envelope_ci95_high, scheme = "SE envelope")]
))
rank_intervals[, model_label_plot := label_models(model_id)]
rank_intervals[, source := label_sources(source)]
rank_intervals[, stratum_plot := label_strata(stratum)]
rank_intervals[, scheme := factor(scheme, levels = c("Document", "Scaffold", "Protein cluster", "Three way", "SE envelope"))]
ps2 <- ggplot(rank_intervals, aes(estimate, model_label_plot, colour = scheme)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = 0.45, position = position_dodge(width = 0.66)) +
  geom_point(size = 1.25, position = position_dodge(width = 0.66)) +
  facet_grid(stratum_plot ~ source, scales = "free_x") +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#111827")) +
  labs(
    title = "Complete rank-correlation dependence sensitivity",
    x = "Spearman rho with 95% interval", y = NULL, colour = "Dependence scheme"
  ) +
  theme_pub() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())
save_figure(ps2, "Supplementary_Figure_2_v1_7_rank_dependence_atlas", supp_dir, 180, 205, tiff = FALSE)
write_source(rank_intervals, "Supplementary_Figure_2_rank_intervals")

# Supplementary Figure 3: model agreement.
agreement_plot <- agreement[stratum == "overall"]
agreement_plot[, model_a_label := factor(unname(model_labels[model_a]), levels = unname(model_labels[model_order]))]
agreement_plot[, model_b_label := factor(unname(model_labels[model_b]), levels = rev(unname(model_labels[model_order])))]
agreement_plot[, source := label_sources(source)]
ps3a <- ggplot(agreement_plot, aes(model_a_label, model_b_label, fill = prediction_spearman)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", prediction_spearman)), size = 2.3) +
  facet_wrap(~source) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-0.2, 1), oob = squish) +
  labs(title = "Prediction-rank agreement", x = NULL, y = NULL, fill = "rho") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
ps3b <- ggplot(agreement_plot, aes(model_a_label, model_b_label, fill = absolute_error_spearman)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", absolute_error_spearman)), size = 2.3) +
  facet_wrap(~source) +
  scale_fill_gradient(low = "white", high = "#D55E00", limits = c(0, 1), oob = squish) +
  labs(title = "Absolute-error agreement", x = NULL, y = NULL, fill = "rho") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
supp3 <- ps3a / ps3b + plot_annotation(tag_levels = "A")
save_figure(supp3, "Supplementary_Figure_3_v1_7_model_agreement", supp_dir, 180, 175, tiff = FALSE)

# Supplementary Figure 4: complete point-accuracy forest.
utility_complete <- copy(utility)
utility_complete[, model_label_plot := label_models(model_id)]
utility_complete[, source := label_sources(source)]
utility_complete[, stratum_plot := label_strata(stratum)]
ps4 <- ggplot(utility_complete, aes(delta_mae, model_label_plot, colour = stratum_plot)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#111827") +
  geom_errorbarh(
    aes(xmin = delta_mae_dependence_envelope_ci95_low, xmax = delta_mae_dependence_envelope_ci95_high),
    height = 0, linewidth = 0.55, position = position_dodge(width = 0.64)
  ) +
  geom_point(size = 1.6, position = position_dodge(width = 0.64)) +
  facet_grid(stratum_plot ~ source, scales = "free_x") +
  scale_colour_manual(values = unname(stratum_colors[stratum_order]), guide = "none") +
  labs(
    title = "Complete incremental point-accuracy atlas",
    subtitle = "Negative values favor the model; conservative dependence-envelope 95% intervals",
    x = expression(Delta*"MAE versus frozen median (pKd)"), y = NULL
  ) +
  theme_pub() + theme(panel.grid.major.y = element_blank())
save_figure(ps4, "Supplementary_Figure_4_v1_7_complete_point_accuracy", supp_dir, 180, 205, tiff = FALSE)

# Supplementary Figure 5: dominant document, protein-cluster, and scaffold influence.
influence_plot <- copy(influence)
influence_plot[, model_label_plot := label_models(model_id)]
influence_plot[, source := label_sources(source)]
influence_plot[, stratum_plot := label_strata(stratum)]
influence_plot[, group_label := factor(
  c(document_component = "Document", protein_identity_cluster = "Protein cluster", ligand_scaffold = "Scaffold")[group_type],
  levels = c("Document", "Protein cluster", "Scaffold")
)]
ps5 <- ggplot(influence_plot, aes(spearman_shift, model_label_plot, colour = group_label)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  geom_point(size = 1.55, position = position_dodge(width = 0.56)) +
  facet_grid(stratum_plot ~ source, scales = "free_x") +
  scale_colour_manual(values = c("Document" = "#D55E00", "Protein cluster" = "#0072B2", "Scaffold" = "#009E73")) +
  labs(
    title = "Largest-cluster deletion changes rank estimates",
    x = "Change in Spearman rho after deletion", y = NULL, colour = "Removed unit"
  ) +
  theme_pub() + theme(legend.position = "top", panel.grid.major.y = element_blank())
save_figure(ps5, "Supplementary_Figure_5_v1_7_cluster_influence", supp_dir, 180, 205, tiff = FALSE)

# Supplementary Figure 6: full conditional coverage atlas.
coverage_plot <- copy(coverage)
coverage_plot[, source := label_sources(source)]
coverage_plot[, subgroup_type_label := factor(
  c(
    cold_category = "Cold category",
    affinity_regime = "Affinity regime",
    prediction_decile = "Prediction decile",
    document_size_quartile = "Document-size quartile"
  )[subgroup_type],
  levels = c("Cold category", "Affinity regime", "Prediction decile", "Document-size quartile")
)]
ps6 <- ggplot(coverage_plot, aes(subgroup, coverage_90, colour = source, group = source)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", linewidth = 0.45, colour = "#111827") +
  geom_point(size = 1.7) +
  facet_grid(subgroup_type_label ~ source, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = source_colors, guide = "none") +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title = "Conditional coverage atlas", x = "Subgroup", y = "Coverage") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank())
save_figure(ps6, "Supplementary_Figure_6_v1_7_conditional_coverage", supp_dir, 180, 205, tiff = FALSE)

# Supplementary Figure 7: recalibration replicate distributions.
budget_replicates[, source := label_sources(source)]
budget_replicates[, budget_label := factor(
  comma(requested_label_budget), levels = comma(sort(unique(requested_label_budget)))
)]
ps7a <- ggplot(budget_replicates, aes(budget_label, test_coverage, fill = source)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", linewidth = 0.45) +
  geom_boxplot(outlier.size = 0.35, linewidth = 0.4, width = 0.68, position = position_dodge(width = 0.74)) +
  scale_fill_manual(values = source_colors) +
  scale_y_continuous(labels = percent, limits = c(0.45, 1)) +
  labs(title = "Coverage across document-disjoint recalibration repeats", x = "Requested labels", y = "Test coverage", fill = "Source") +
  theme_pub() + theme(legend.position = "top")
ps7b <- ggplot(budget_replicates, aes(budget_label, half_width, fill = source)) +
  geom_hline(yintercept = 1.455871, linetype = "dashed", linewidth = 0.45) +
  geom_boxplot(outlier.size = 0.35, linewidth = 0.4, width = 0.68, position = position_dodge(width = 0.74)) +
  scale_fill_manual(values = source_colors) +
  labs(title = "Half-width across recalibration repeats", x = "Requested labels", y = "Half-width (pKd)", fill = "Source") +
  theme_pub() + theme(legend.position = "none")
supp7 <- ps7a / ps7b + plot_annotation(tag_levels = "A")
save_figure(supp7, "Supplementary_Figure_7_v1_7_recalibration_repeats", supp_dir, 180, 165, tiff = FALSE)

# Supplementary Figure 8: complete identity-permutation distributions.
v15_null[, stratum_plot := factor(
  c(
    all_exact_kd = "Overall",
    drug_cold = "Ligand cold",
    protein_identity_cold = "Target cold",
    protein_identity_cold_including_double = "Target incl. double",
    seen_molecule_and_protein_cluster_new_pair = "New relation",
    strict_double_cold = "Double cold"
  )[stratum]
)]
v15_null[, control_plot := factor(
  unname(null_control_labels[control]), levels = unname(null_control_labels)
)]
observed_null <- v15_metrics[, .(
  stratum,
  stratum_plot = factor(c(
    all_exact_kd = "Overall", drug_cold = "Ligand cold", protein_identity_cold = "Target cold",
    protein_identity_cold_including_double = "Target incl. double",
    seen_molecule_and_protein_cluster_new_pair = "New relation", strict_double_cold = "Double cold"
  )[stratum]),
  observed_spearman = spearman
)]
ps8 <- ggplot(v15_null, aes(spearman, fill = control_plot, colour = control_plot)) +
  geom_density(alpha = 0.16, linewidth = 0.55) +
  geom_vline(
    data = observed_null,
    aes(xintercept = observed_spearman), inherit.aes = FALSE,
    linewidth = 0.55, colour = "#111827"
  ) +
  facet_wrap(~stratum_plot, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Protein" = "#56B4E9", "Ligand" = "#E69F00", "Joint" = "#6B7280")) +
  scale_colour_manual(values = c("Protein" = "#56B4E9", "Ligand" = "#E69F00", "Joint" = "#6B7280")) +
  labs(
    title = "Complete identity-permutation null distributions",
    subtitle = "Black vertical line: observed Spearman rho; 100 permutations per control and stratum",
    x = "Spearman rho", y = "Density", fill = "Control", colour = "Control"
  ) +
  theme_pub() + theme(legend.position = "top")
save_figure(ps8, "Supplementary_Figure_8_v1_7_permutation_nulls", supp_dir, 180, 150, tiff = FALSE)

# Supplementary Figure 9: full assay-context and geometry controls.
assay_all <- assay_forest[, .(
  stream = "Assay context",
  comparison = paste(category_label, "vs", reference_label),
  estimate = mae_difference_candidate_minus_reference,
  low = ci_95_low,
  high = ci_95_high
)]
geometry_all <- geometry_controls[, .(
  stream = "Geometry controls",
  comparison = paste(challenge_label, control_label, sep = ": "),
  estimate = delta_mae,
  low = ci95_lower,
  high = ci95_upper
)]
ps9a <- ggplot(assay_all, aes(estimate, reorder(comparison, estimate), colour = stream)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = 0.65) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c("Assay context" = "#009E73"), guide = "none") +
  labs(title = "Assay-context contrasts", x = expression(Delta*"MAE (pKd)"), y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank())
ps9b <- ggplot(geometry_all, aes(estimate, reorder(comparison, estimate), colour = stream)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = 0.65) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c("Geometry controls" = "#CC79A7"), guide = "none") +
  labs(title = "Geometry and matched controls", x = expression(Delta*"MAE (pKd)"), y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_blank())
supp9 <- ps9a | ps9b + plot_annotation(tag_levels = "A")
save_figure(supp9, "Supplementary_Figure_9_v1_7_context_geometry_controls", supp_dir, 180, 170, tiff = FALSE)

# Supplementary Figure 10: data-quality and human-review status.
quality <- profile[, .(
  source,
  `Outside pKd 2-12` = outcomes_below_2 + outcomes_above_12,
  `Multiple measurements / candidates` = round(multi_measurement_fraction * n),
  `Nonzero aggregated range` = round(nonzero_label_range_fraction * n)
)]
quality <- melt(quality, id.vars = "source", variable.name = "flag", value.name = "rows")
quality[, source := label_sources(source)]
ps10a <- ggplot(quality, aes(flag, rows, fill = source)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.64) +
  geom_text(aes(label = rows), position = position_dodge(width = 0.74), vjust = -0.35, size = 2.25) +
  scale_fill_manual(values = source_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Provenance-risk flags", x = NULL, y = "Rows", fill = "Source") +
  theme_pub() + theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.major.x = element_blank())

review_status <- data.table(
  queue = factor(c("Frozen random sample", "Risk-stratified sample"), levels = c("Frozen random sample", "Risk-stratified sample")),
  machine_reconstructed = c(100, 100),
  independently_signed = c(0, 0)
)
review_long <- melt(review_status, id.vars = "queue", variable.name = "stage", value.name = "pairs")
review_long[, stage_label := factor(
  c(machine_reconstructed = "Machine prepared", independently_signed = "Independent signed decision")[stage],
  levels = c("Machine prepared", "Independent signed decision")
)]
ps10b <- ggplot(review_long, aes(queue, pairs, fill = stage_label)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = pairs), position = position_dodge(width = 0.72), vjust = -0.35, size = 2.3) +
  scale_fill_manual(values = c("Machine prepared" = "#56B4E9", "Independent signed decision" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 110), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Machine reconstruction does not close source adjudication",
    x = NULL, y = "Pairs", fill = "Review stage"
  ) +
  theme_pub() + theme(legend.position = "top", panel.grid.major.x = element_blank())

descriptor_top <- descriptor_support[order(-row_fraction)][1:min(12, .N)]
descriptor_top[, descriptor_label := factor(descriptor_label, levels = rev(descriptor_label))]
ps10c <- ggplot(descriptor_top, aes(row_fraction, descriptor_label, colour = family)) +
  geom_segment(aes(x = 0, xend = row_fraction, yend = descriptor_label), linewidth = 0.65, alpha = 0.6) +
  geom_point(size = 2) +
  scale_x_continuous(labels = percent) +
  labs(title = "Assay-description concepts are unevenly supported", x = "Rows positive", y = NULL, colour = "Descriptor family") +
  theme_pub() + theme(legend.position = "top", panel.grid.major.y = element_blank())

supp10 <- (ps10a | ps10b) / ps10c +
  plot_layout(heights = c(0.9, 1.1)) + plot_annotation(tag_levels = "A")
save_figure(supp10, "Supplementary_Figure_10_v1_7_data_quality_and_review", supp_dir, 180, 165, tiff = FALSE)

# Supplementary Figure 11: source poststratification and entity-prior details.
ps11a <- ggplot(post_plot, aes(rank_correlation, model_label_plot, colour = analysis_label, shape = analysis_label)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#111827") +
  geom_point(size = 2, position = position_dodge(width = 0.48)) +
  facet_wrap(~source, scales = "free_x") +
  scale_colour_manual(values = c("Raw common support" = "#6B7280", "Poststratified" = "#0072B2")) +
  scale_shape_manual(values = c("Raw common support" = 16, "Poststratified" = 17)) +
  labs(title = "Rank sensitivity to observed-source poststratification", x = "Weighted rank correlation", y = NULL, colour = NULL, shape = NULL) +
  theme_pub() + theme(legend.position = "top", panel.grid.major.y = element_blank())

entity_long <- melt(
  entity_priors,
  id.vars = c("source", "stratum", "n"),
  measure.vars = c(
    "global_median_mae", "exact_molecule_prior_mae", "exact_protein_prior_mae",
    "protein_cluster_prior_mae", "combined_entity_prior_mae"
  ),
  variable.name = "baseline", value.name = "mae"
)
entity_long[, baseline_label := factor(c(
  global_median_mae = "Global median",
  exact_molecule_prior_mae = "Exact molecule prior",
  exact_protein_prior_mae = "Exact protein prior",
  protein_cluster_prior_mae = "Protein-cluster prior",
  combined_entity_prior_mae = "Combined prior"
)[baseline])]
entity_long[, source := label_sources(source)]
entity_long[, stratum_plot := label_strata(stratum)]
ps11b <- ggplot(entity_long, aes(mae, baseline_label, colour = stratum_plot)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.55)) +
  facet_wrap(~source, scales = "free_x") +
  scale_colour_manual(values = unname(stratum_colors[stratum_order])) +
  labs(title = "Development-only entity-prior baseline MAE", x = "MAE (pKd)", y = NULL, colour = "Stratum") +
  theme_pub() + theme(legend.position = "top", panel.grid.major.y = element_blank())
supp11 <- ps11a / ps11b + plot_annotation(tag_levels = "A")
save_figure(supp11, "Supplementary_Figure_11_v1_7_standardization_and_priors", supp_dir, 180, 165, tiff = FALSE)

session_path <- file.path(qc_dir, "V1_7_R_SESSION_INFO.txt")
writeLines(capture.output(sessionInfo()), session_path)
build_status <- list(
  schema_version = "science_advances_v1_7_figure_build_v1",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  main_figures = 6,
  supplementary_figures = 11,
  main_width_mm = 180,
  main_height_range_mm = c(148, 178),
  minimum_theme_text_pt = 6.5,
  formats_main = c("svg", "pdf", "png_450dpi", "tiff_600dpi"),
  formats_supplementary = c("svg", "pdf", "png_450dpi")
)
jsonlite::write_json(
  build_status,
  file.path(qc_dir, "V1_7_FIGURE_BUILD_STATUS.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("v1.7 main and supplementary figures written to: ", package_root)
