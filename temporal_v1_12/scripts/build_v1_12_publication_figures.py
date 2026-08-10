#!/usr/bin/env python3
"""Build the pre-specified v1.12 main and supplementary figures."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


PRIMARY_MODELS = [
    "global_top5_morgan_similarity_v1_8",
    "hgb_fusion_additive_pre2018_v1_8",
    "ligand_mpnn_frozen_esm2_v1_8",
    "smiles_protein_cross_attention_v1_8",
]
MODEL_LABELS = {
    "global_top5_morgan_similarity_v1_8": "Morgan local",
    "hgb_fusion_additive_pre2018_v1_8": "HGB fusion",
    "ligand_mpnn_frozen_esm2_v1_8": "Ligand MPNN + ESM2",
    "smiles_protein_cross_attention_v1_8": "Cross-attention",
}
MODEL_COLORS = {
    "global_top5_morgan_similarity_v1_8": "#1F4E79",
    "hgb_fusion_additive_pre2018_v1_8": "#7A5195",
    "ligand_mpnn_frozen_esm2_v1_8": "#3B82A0",
    "smiles_protein_cross_attention_v1_8": "#A05195",
}
SOURCE_LABELS = {
    "strict_dual_timestamp_patent": "Strict-date patent",
    "database_entry_patent": "Entry-time patent",
    "database_entry_residual_patent": "Residual patent",
    "database_entry_chembl_import": "ChEMBL import",
    "database_entry_pubchem_assay": "PubChem assay",
    "database_entry_curated_article": "Curated article",
    "database_entry_other_documented": "Other documented",
}
SENSITIVITY_LABELS = {
    "exclude_exact_10000_nm": "Exclude 10,000 nM",
    "leave_one_origin_out__database_entry_chembl_import": "Exclude ChEMBL import",
    "leave_one_origin_out__database_entry_patent": "Exclude entry-time patent",
    "leave_one_origin_out__database_entry_residual_patent": "Exclude residual patent",
    "leave_one_origin_out__strict_dual_timestamp_patent": "Exclude strict-date patent",
    "nested_document_cap__10": "Document cap 10",
    "nested_document_cap__15": "Document cap 15",
    "nested_document_cap__20": "Document cap 20",
    "publication_to_entry_lag__0_90_days": "Entry lag 0-90 d",
    "publication_to_entry_lag__91_365_days": "Entry lag 91-365 d",
    "publication_to_entry_lag__over_365_days": "Entry lag >365 d",
    "publication_to_entry_lag__unknown_or_negative": "Entry lag unknown/negative",
    "replicated_pairs": "Replicated pairs",
    "single_measurement_pairs": "Single-measurement pairs",
    "source_origin__database_entry_chembl_import": "ChEMBL import only",
    "source_origin__database_entry_patent": "Entry-time patent only",
    "source_origin__database_entry_residual_patent": "Residual patent only",
    "source_origin__strict_dual_timestamp_patent": "Strict-date patent only",
    "strict_publication_and_entry_post_2024": "Strict post-2024",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected JSON object: {path}")
    return value


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.2,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.2,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.3,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def clean_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.12,
        1.06,
        label,
        transform=axis.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
        ha="left",
    )


def padded_limits(values: np.ndarray, include: tuple[float, ...]) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    finite = np.r_[finite, np.asarray(include, dtype=float)]
    low = float(np.min(finite))
    high = float(np.max(finite))
    span = max(high - low, 0.2)
    return low - 0.08 * span, high + 0.08 * span


def save_figure(figure: plt.Figure, stem: Path) -> list[Path]:
    outputs: list[Path] = []
    for suffix, dpi in ((".pdf", 600), (".svg", 600), (".png", 600), (".tiff", 600)):
        path = stem.with_suffix(suffix)
        figure.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.03)
        outputs.append(path)
    return outputs


def gate_panel(axis: plt.Axes, manifest: dict[str, Any]) -> pd.DataFrame:
    gate = manifest["structural_adequacy_gate"]
    required = gate["required"]
    observed = gate["observed"]
    labels = {
        "pairs": "Pairs",
        "document_components": "Document components",
        "effective_document_components": "Effective components",
        "exact_proteins": "Exact proteins",
        "eligible_within_target_groups": "Eligible targets",
    }
    order = list(labels)
    ratios = np.asarray([float(observed[key]) / float(required[key]) for key in order])
    y = np.arange(len(order))[::-1]
    colors = np.where(ratios >= 1, "#3F7D6B", "#9B4B4B")
    axis.axvline(1.0, color="#374151", lw=0.8, ls="--")
    axis.hlines(y, 0, ratios, color="#C7C7C7", lw=1.2)
    axis.scatter(ratios, y, c=colors, s=32, edgecolor="white", linewidth=0.5, zorder=3)
    def format_count(key: str, value: Any) -> str:
        if key == "effective_document_components":
            return f"{float(value):,.1f}"
        return f"{int(value):,}"

    for ratio, position, key in zip(ratios, y, order):
        axis.text(
            ratio,
            position + 0.15,
            f"{format_count(key, observed[key])}/{format_count(key, required[key])}\n{ratio:.2f}x",
            ha="center",
            va="bottom",
            fontsize=5.2,
            color="#333333",
        )
    axis.set_yticks(y, [labels[key] for key in order])
    axis.set_xlabel("Observed / minimum (ratio)")
    axis.set_xlim(0, max(1.25, float(np.max(ratios)) * 1.25))
    axis.set_title("Structure-adapted gate", loc="left", fontweight="bold")
    clean_axis(axis)
    return pd.DataFrame(
        {
            "gate": order,
            "required": [required[key] for key in order],
            "observed": [observed[key] for key in order],
            "observed_to_required": ratios,
            "passed": [bool(gate["checks"][key]) for key in order],
        }
    )


def paired_rank_panel(axis: plt.Axes, primary: pd.DataFrame) -> pd.DataFrame:
    endpoints = ["source_wide_spearman", "within_target_centered_spearman"]
    subset = primary.loc[primary["endpoint"].isin(endpoints)].copy()
    y_base = np.arange(len(PRIMARY_MODELS))[::-1]
    offsets = {endpoints[0]: 0.12, endpoints[1]: -0.12}
    styles = {
        endpoints[0]: ("o", "Source-wide"),
        endpoints[1]: ("s", "Within target"),
    }
    axis.axvline(0, color="#555555", lw=0.7, ls="--")
    for endpoint in endpoints:
        rows = subset.set_index(["model_id", "endpoint"])
        estimates = np.asarray(
            [rows.loc[(model, endpoint), "estimate"] for model in PRIMARY_MODELS],
            dtype=float,
        )
        lows = np.asarray(
            [rows.loc[(model, endpoint), "simultaneous_95_ci_low"] for model in PRIMARY_MODELS],
            dtype=float,
        )
        highs = np.asarray(
            [rows.loc[(model, endpoint), "simultaneous_95_ci_high"] for model in PRIMARY_MODELS],
            dtype=float,
        )
        positions = y_base + offsets[endpoint]
        for index, model in enumerate(PRIMARY_MODELS):
            axis.errorbar(
                estimates[index],
                positions[index],
                xerr=np.asarray(
                    [[estimates[index] - lows[index]], [highs[index] - estimates[index]]]
                ),
                fmt=styles[endpoint][0],
                color=MODEL_COLORS[model],
                mfc=MODEL_COLORS[model] if endpoint == endpoints[0] else "white",
                mec=MODEL_COLORS[model],
                ms=4.2,
                mew=0.8,
                lw=0.9,
                capsize=1.8,
            )
    all_values = subset[
        ["estimate", "simultaneous_95_ci_low", "simultaneous_95_ci_high"]
    ].to_numpy(float).ravel()
    axis.set_xlim(*padded_limits(all_values, (0.0,)))
    axis.set_yticks(y_base, [MODEL_LABELS[model] for model in PRIMARY_MODELS])
    axis.set_xlabel("Spearman rho (simultaneous 95% CI)")
    axis.set_title("Rank association by scale", loc="left", fontweight="bold")
    axis.legend(
        handles=[
            Line2D(
                [], [], marker="o", color="#374151", markerfacecolor="#374151",
                markeredgecolor="#374151", linestyle="None", markersize=4.3,
                label="Source-wide",
            ),
            Line2D(
                [], [], marker="s", color="#374151", markerfacecolor="white",
                markeredgecolor="#374151", linestyle="None", markersize=4.3,
                label="Within target",
            ),
        ],
        loc="lower left",
        bbox_to_anchor=(0.01, 0.96),
        frameon=False,
        handletextpad=0.35,
        borderpad=0.2,
        ncol=2,
    )
    clean_axis(axis)
    return subset


def attenuation_panel(axis: plt.Axes, primary: pd.DataFrame) -> pd.DataFrame:
    subset = primary.loc[primary["endpoint"].eq("target_rank_attenuation")].copy()
    rows = subset.set_index("model_id")
    y = np.arange(len(PRIMARY_MODELS))[::-1]
    axis.axvline(0, color="#555555", lw=0.7, ls="--")
    axis.axvline(0.10, color="#7A5195", lw=0.9, ls=":")
    values = []
    for position, model in zip(y, PRIMARY_MODELS):
        row = rows.loc[model]
        estimate = float(row["estimate"])
        low = float(row["simultaneous_95_ci_low"])
        high = float(row["simultaneous_95_ci_high"])
        values.extend([estimate, low, high])
        axis.errorbar(
            estimate,
            position,
            xerr=np.asarray([[estimate - low], [high - estimate]]),
            fmt="o",
            color=MODEL_COLORS[model],
            ms=4.5,
            lw=1.0,
            capsize=2,
        )
    low_limit, high_limit = padded_limits(np.asarray(values), (0.0, 0.10))
    axis.set_xlim(low_limit, high_limit)
    axis.set_yticks(y, [MODEL_LABELS[model] for model in PRIMARY_MODELS])
    axis.set_xlabel("Source-wide minus within-target rho")
    axis.set_title("Pre-specified attenuation", loc="left", fontweight="bold")
    clean_axis(axis)
    return subset


def source_heatmap_panel(axis: plt.Axes, sensitivities: pd.DataFrame) -> pd.DataFrame:
    subset = sensitivities.loc[
        sensitivities["category"].isin(
            {"mandatory_source_origin", "mandatory_leave_one_origin_out"}
        )
    ].copy()
    order = []
    for source in SOURCE_LABELS:
        for prefix in ("source_origin__", "leave_one_origin_out__"):
            key = prefix + source
            if subset["sensitivity"].eq(key).any():
                order.append(key)
    matrix = (
        subset.pivot(index="sensitivity", columns="model_id", values="target_rank_attenuation")
        .reindex(index=order, columns=PRIMARY_MODELS)
    )
    values = matrix.to_numpy(float)
    finite = np.abs(values[np.isfinite(values)])
    bound = max(0.15, float(np.max(finite)) if len(finite) else 0.15)
    image = axis.imshow(
        values,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound),
        interpolation="none",
    )
    labels = []
    for key in order:
        if key.startswith("source_origin__"):
            labels.append(SOURCE_LABELS.get(key.split("__", 1)[1], key) + " only")
        else:
            labels.append("Exclude " + SOURCE_LABELS.get(key.split("__", 1)[1], key))
    axis.set_yticks(np.arange(len(order)), labels)
    axis.set_xticks(
        np.arange(len(PRIMARY_MODELS)),
        [MODEL_LABELS[model] for model in PRIMARY_MODELS],
        rotation=35,
        ha="right",
    )
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = values[row_index, column_index]
            if np.isfinite(value):
                axis.text(
                    column_index,
                    row_index,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=5.4,
                    color="white" if abs(value) > 0.55 * bound else "#222222",
                )
    colorbar = axis.figure.colorbar(image, ax=axis, fraction=0.035, pad=0.025)
    colorbar.set_label("Rank attenuation", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=5.8)
    axis.set_title("Source-origin tests", loc="left", fontweight="bold")
    return subset


def build_main(
    root: Path,
    manifest: dict[str, Any],
    primary: pd.DataFrame,
    sensitivities: pd.DataFrame,
    source_data_dir: Path,
) -> list[Path]:
    figure = plt.figure(figsize=(183 / 25.4, 154 / 25.4), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=[0.92, 1.25], width_ratios=[0.9, 1.1])
    axes = [
        figure.add_subplot(grid[0, 0]),
        figure.add_subplot(grid[0, 1]),
        figure.add_subplot(grid[1, 0]),
        figure.add_subplot(grid[1, 1]),
    ]
    gate_data = gate_panel(axes[0], manifest)
    rank_data = paired_rank_panel(axes[1], primary)
    attenuation_data = attenuation_panel(axes[2], primary)
    source_data = source_heatmap_panel(axes[3], sensitivities)
    for axis, label in zip(axes, "ABCD"):
        panel_label(axis, label)
    gate_data.to_csv(source_data_dir / "Figure_T1A_evidence_gate.tsv", sep="\t", index=False)
    rank_data.to_csv(source_data_dir / "Figure_T1B_rank_scale.tsv", sep="\t", index=False)
    attenuation_data.to_csv(source_data_dir / "Figure_T1C_attenuation.tsv", sep="\t", index=False)
    source_data.to_csv(source_data_dir / "Figure_T1D_source_stress.tsv", sep="\t", index=False)
    outputs = save_figure(
        figure, root / "figures" / "main" / "Figure_T1_v1_12_temporal_architecture_test"
    )
    outputs.extend(
        save_figure(
            figure,
            root / "figures" / "main" / "Figure_6_outcome_blind_temporal_validation",
        )
    )
    plt.close(figure)
    return outputs


def build_supplement(
    root: Path,
    manifest: dict[str, Any],
    blind: pd.DataFrame,
    target_support: pd.DataFrame,
    secondary: pd.DataFrame,
    covariance: pd.DataFrame,
    sensitivities: pd.DataFrame,
    source_data_dir: Path,
) -> list[Path]:
    figure = plt.figure(figsize=(183 / 25.4, 205 / 25.4), constrained_layout=True)
    grid = figure.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.15])
    axes = [figure.add_subplot(grid[row, column]) for row in range(3) for column in range(2)]

    source_counts = (
        blind["source_origin_stratum"].value_counts().rename_axis("source_origin").reset_index(name="pairs")
    )
    source_counts["label"] = source_counts["source_origin"].map(SOURCE_LABELS).fillna(source_counts["source_origin"])
    source_counts = source_counts.sort_values("pairs")
    axes[0].barh(source_counts["label"], source_counts["pairs"], color="#56B4E9")
    axes[0].set_xlabel("Frozen pairs")
    axes[0].set_title("Source composition", loc="left", fontweight="bold")
    clean_axis(axes[0])

    eligible = target_support["eligible_within_target"].astype(str).str.lower().eq("true")
    axes[1].scatter(
        target_support.loc[~eligible, "unique_ligands"],
        target_support.loc[~eligible, "unique_scaffolds"],
        s=8,
        color="#BDBDBD",
        alpha=0.65,
        linewidth=0,
        rasterized=True,
    )
    axes[1].scatter(
        target_support.loc[eligible, "unique_ligands"],
        target_support.loc[eligible, "unique_scaffolds"],
        s=11,
        color="#009E73",
        alpha=0.8,
        linewidth=0,
        rasterized=True,
    )
    axes[1].axvline(5, color="#555555", lw=0.7, ls="--")
    axes[1].axhline(3, color="#555555", lw=0.7, ls="--")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Unique ligands per exact protein")
    axes[1].set_ylabel("Unique scaffolds")
    axes[1].set_title("Within-target support", loc="left", fontweight="bold")
    clean_axis(axes[1])

    covariance_primary = covariance.loc[
        covariance["model_id"].isin(PRIMARY_MODELS)
        & covariance["grouping"].eq("exact_protein")
    ].set_index("model_id").reindex(PRIMARY_MODELS)
    y = np.arange(len(PRIMARY_MODELS))[::-1]
    axes[2].barh(
        y,
        covariance_primary["between_group_rank_covariance"],
        color="#0072B2",
        label="Between protein",
    )
    axes[2].barh(
        y,
        covariance_primary["within_group_rank_covariance"],
        left=covariance_primary["between_group_rank_covariance"],
        color="#E69F00",
        label="Within protein",
    )
    axes[2].set_yticks(y, [MODEL_LABELS[model] for model in PRIMARY_MODELS])
    axes[2].set_xlabel("Rank covariance")
    axes[2].set_title("Exact-protein covariance decomposition", loc="left", fontweight="bold")
    clean_axis(axes[2])

    secondary_primary = secondary.set_index("model_id").reindex(PRIMARY_MODELS)
    metric_columns = [
        "mae",
        "calibration_slope_observed_on_predicted",
        "dispersion_ratio_predicted_to_observed",
        "interval_90_coverage",
        "best_ligand_selection_regret_pkd",
    ]
    metric_labels = ["MAE", "Slope", "Dispersion", "90% coverage", "Best-ligand regret"]
    normalized = secondary_primary[metric_columns].astype(float)
    z = (normalized - normalized.mean(axis=0)) / normalized.std(axis=0, ddof=0).replace(0, np.nan)
    image = axes[3].imshow(z.to_numpy(float), aspect="auto", cmap="cividis", interpolation="none")
    axes[3].set_yticks(np.arange(len(PRIMARY_MODELS)), [MODEL_LABELS[model] for model in PRIMARY_MODELS])
    axes[3].set_xticks(np.arange(len(metric_labels)), metric_labels, rotation=35, ha="right")
    axes[3].set_title("Calibration and selection diagnostics", loc="left", fontweight="bold")
    figure.colorbar(image, ax=axes[3], fraction=0.035, pad=0.025, label="Column z-score")

    sensitivity_matrix = (
        sensitivities.pivot(index="sensitivity", columns="model_id", values="target_rank_attenuation")
        .reindex(columns=PRIMARY_MODELS)
    )
    values = sensitivity_matrix.to_numpy(float)
    finite = np.abs(values[np.isfinite(values)])
    bound = max(0.15, float(np.max(finite)) if len(finite) else 0.15)
    image = axes[4].imshow(
        values,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-bound, vcenter=0, vmax=bound),
        interpolation="none",
    )
    sensitivity_labels = [
        SENSITIVITY_LABELS.get(value, value.replace("__", ": "))
        for value in sensitivity_matrix.index
    ]
    axes[4].set_yticks(np.arange(len(sensitivity_matrix.index)), sensitivity_labels)
    axes[4].set_xticks(np.arange(len(PRIMARY_MODELS)), [MODEL_LABELS[model] for model in PRIMARY_MODELS], rotation=35, ha="right")
    axes[4].set_title("All locked attenuation sensitivities", loc="left", fontweight="bold")
    figure.colorbar(image, ax=axes[4], fraction=0.035, pad=0.025, label="Rank attenuation")

    source_components = blind.groupby("document_component_id").size().sort_values(ascending=False)
    cumulative = source_components.cumsum() / source_components.sum()
    axes[5].plot(np.arange(1, len(cumulative) + 1), cumulative, color="#0072B2", lw=1.3)
    axes[5].axhline(0.5, color="#555555", lw=0.7, ls="--")
    axes[5].axhline(0.9, color="#555555", lw=0.7, ls=":")
    axes[5].set_xscale("log")
    axes[5].set_ylim(0, 1.02)
    axes[5].set_xlabel("Document components, ranked by size")
    axes[5].set_ylabel("Cumulative pair fraction")
    axes[5].set_title("Dependence concentration", loc="left", fontweight="bold")
    clean_axis(axes[5])

    for axis, label in zip(axes, "ABCDEF"):
        panel_label(axis, label)
    source_counts.to_csv(source_data_dir / "Figure_TS1A_source_composition.tsv", sep="\t", index=False)
    target_support.to_csv(source_data_dir / "Figure_TS1B_target_support.tsv", sep="\t", index=False)
    covariance_primary.reset_index().to_csv(source_data_dir / "Figure_TS1C_covariance.tsv", sep="\t", index=False)
    secondary_primary.reset_index().to_csv(source_data_dir / "Figure_TS1D_calibration_utility.tsv", sep="\t", index=False)
    sensitivity_matrix.reset_index().to_csv(source_data_dir / "Figure_TS1E_all_sensitivities.tsv", sep="\t", index=False)
    source_components.rename("pairs").rename_axis("document_component_id").reset_index().to_csv(
        source_data_dir / "Figure_TS1F_component_concentration.tsv", sep="\t", index=False
    )
    outputs = save_figure(
        figure,
        root / "figures" / "supplementary" / "Figure_TS1_v1_12_temporal_architecture_diagnostics",
    )
    outputs.extend(
        save_figure(
            figure,
            root / "figures" / "supplementary" / "Figure_S12_temporal_architecture_diagnostics",
        )
    )
    plt.close(figure)
    return outputs


def annotated_matrix(
    axis: plt.Axes,
    matrix: pd.DataFrame,
    title: str,
    colorbar_label: str,
) -> None:
    values = matrix.to_numpy(float)
    finite = np.abs(values[np.isfinite(values)])
    bound = max(0.15, float(np.max(finite)) if len(finite) else 0.15)
    image = axis.imshow(
        values,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound),
        interpolation="none",
    )
    axis.set_yticks(
        np.arange(len(matrix.index)),
        [SOURCE_LABELS.get(value, value) for value in matrix.index],
    )
    axis.set_xticks(
        np.arange(len(PRIMARY_MODELS)),
        [MODEL_LABELS[model] for model in PRIMARY_MODELS],
        rotation=35,
        ha="right",
    )
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = values[row_index, column_index]
            if np.isfinite(value):
                axis.text(
                    column_index,
                    row_index,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=5.5,
                    color="white" if abs(value) > 0.55 * bound else "#222222",
                )
    colorbar = axis.figure.colorbar(image, ax=axis, fraction=0.035, pad=0.025)
    colorbar.set_label(colorbar_label, fontsize=6.4)
    colorbar.ax.tick_params(labelsize=5.8)
    axis.set_title(title, loc="left", fontweight="bold")


def build_resampling_supplement(
    root: Path,
    primary: pd.DataFrame,
    sensitivities: pd.DataFrame,
    resampling: pd.DataFrame,
    source_intervals: pd.DataFrame,
    source_data_dir: Path,
) -> list[Path]:
    figure = plt.figure(figsize=(183 / 25.4, 183 / 25.4), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=[1.25, 1.0])
    axes = [
        figure.add_subplot(grid[0, :]),
        figure.add_subplot(grid[1, 0]),
        figure.add_subplot(grid[1, 1]),
    ]

    document = primary.loc[
        primary["endpoint"].eq("target_rank_attenuation"),
        [
            "model_id",
            "estimate",
            "simultaneous_95_ci_low",
            "simultaneous_95_ci_high",
        ],
    ].rename(
        columns={
            "simultaneous_95_ci_low": "ci_low",
            "simultaneous_95_ci_high": "ci_high",
        }
    )
    document["resampling_unit"] = "document_component"
    document["interval_type"] = "simultaneous"
    alternatives = resampling.loc[
        resampling["endpoint"].eq("target_rank_attenuation"),
        [
            "model_id",
            "resampling_unit",
            "estimate",
            "ordinary_95_ci_low",
            "ordinary_95_ci_high",
        ],
    ].rename(
        columns={
            "ordinary_95_ci_low": "ci_low",
            "ordinary_95_ci_high": "ci_high",
        }
    )
    alternatives["interval_type"] = "ordinary"
    forest = pd.concat([document, alternatives], ignore_index=True)
    unit_order = ["document_component", "exact_protein", "scaffold", "source_origin"]
    unit_labels = {
        "document_component": "Document",
        "exact_protein": "Exact protein",
        "scaffold": "Scaffold",
        "source_origin": "Source origin",
    }
    forest["unit_order"] = forest["resampling_unit"].map(
        {value: index for index, value in enumerate(unit_order)}
    )
    forest["model_order"] = forest["model_id"].map(
        {value: index for index, value in enumerate(PRIMARY_MODELS)}
    )
    forest = forest.sort_values(["unit_order", "model_order"]).reset_index(drop=True)
    y = np.arange(len(forest))[::-1]
    axes[0].axvspan(0.10, 0.50, color="#E6F4EA", alpha=0.65, lw=0)
    axes[0].axvline(0, color="#555555", lw=0.7, ls="--")
    axes[0].axvline(0.10, color="#009E73", lw=0.9, ls=":")
    for position, row in zip(y, forest.itertuples(index=False)):
        axes[0].errorbar(
            float(row.estimate),
            position,
            xerr=np.asarray(
                [[float(row.estimate) - float(row.ci_low)], [float(row.ci_high) - float(row.estimate)]]
            ),
            fmt="o",
            color=MODEL_COLORS[str(row.model_id)],
            ms=3.8,
            lw=0.8,
            capsize=1.6,
        )
    axes[0].set_yticks(
        y,
        [
            f"{unit_labels[row.resampling_unit]} | {MODEL_LABELS[row.model_id]}"
            for row in forest.itertuples(index=False)
        ],
    )
    axes[0].set_xlabel("Source-wide minus within-target rho (95% CI)")
    axes[0].set_title(
        "Model-specific attenuation across resampling units",
        loc="left",
        fontweight="bold",
    )
    clean_axis(axes[0])

    strict = sensitivities.loc[
        sensitivities["sensitivity"].eq("strict_publication_and_entry_post_2024")
    ].set_index("model_id").reindex(PRIMARY_MODELS)
    strict_y = np.arange(len(PRIMARY_MODELS))[::-1]
    axes[1].axvline(0, color="#555555", lw=0.7, ls="--")
    for position, model in zip(strict_y, PRIMARY_MODELS):
        source = float(strict.loc[model, "source_wide_spearman"])
        centered = float(strict.loc[model, "within_target_centered_spearman"])
        axes[1].plot([source, centered], [position, position], color="#AFAFAF", lw=1.1)
        axes[1].scatter(source, position, marker="o", s=25, color=MODEL_COLORS[model], zorder=3)
        axes[1].scatter(
            centered,
            position,
            marker="s",
            s=23,
            facecolor="white",
            edgecolor=MODEL_COLORS[model],
            linewidth=0.9,
            zorder=3,
        )
    axes[1].set_yticks(strict_y, [MODEL_LABELS[model] for model in PRIMARY_MODELS])
    axes[1].set_xlabel("Spearman rho")
    axes[1].set_title(
        "Strict publication + entry time subset (n=1,915)",
        loc="left",
        fontweight="bold",
    )
    axes[1].legend(
        handles=[
            mpl.lines.Line2D([], [], marker="o", color="#333333", lw=0, label="Source-wide", markersize=4),
            mpl.lines.Line2D([], [], marker="s", mfc="white", mec="#333333", color="#333333", lw=0, label="Within target", markersize=4),
        ],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
    )
    clean_axis(axes[1])

    source_only = source_intervals.loc[
        source_intervals["category"].eq("mandatory_source_origin")
    ].copy()
    source_only["source_origin"] = source_only["sensitivity"].str.replace(
        "source_origin__", "", regex=False
    )
    attenuation_matrix = (
        source_only.loc[source_only["endpoint"].eq("target_rank_attenuation")]
        .pivot(index="source_origin", columns="model_id", values="estimate")
        .reindex(columns=PRIMARY_MODELS)
    )
    annotated_matrix(
        axes[2],
        attenuation_matrix,
        "Source-origin heterogeneity",
        "Rank attenuation",
    )

    for axis, label in zip(axes, "ABC"):
        panel_label(axis, label)
    forest.drop(columns=["unit_order", "model_order"]).to_csv(
        source_data_dir / "Figure_TS2A_resampling_forest.tsv", sep="\t", index=False
    )
    strict.reset_index().to_csv(
        source_data_dir / "Figure_TS2B_strict_temporal_subset.tsv", sep="\t", index=False
    )
    attenuation_matrix.reset_index().to_csv(
        source_data_dir / "Figure_TS2C_source_origin_attenuation.tsv", sep="\t", index=False
    )
    outputs = save_figure(
        figure,
        root
        / "figures"
        / "supplementary"
        / "Figure_TS2_v1_12_resampling_and_source_heterogeneity",
    )
    outputs.extend(
        save_figure(
            figure,
            root
            / "figures"
            / "supplementary"
            / "Figure_S13_resampling_and_source_heterogeneity",
        )
    )
    plt.close(figure)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser()
    default_root = Path(__file__).resolve().parent.parent
    parser.add_argument("--package-root", type=Path, default=default_root)
    args = parser.parse_args()
    root = args.package_root.resolve()
    complete_path = root / "results_locked" / "V1_12_ANALYSIS_COMPLETE.json"
    if not complete_path.is_file():
        raise RuntimeError("Locked analysis is not complete")
    manifest = load_json(root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json")
    primary = pd.read_csv(root / "results_locked" / "V1_12_PRIMARY_ENDPOINTS.tsv", sep="\t")
    sensitivities = pd.read_csv(
        root / "results_locked" / "V1_12_LOCKED_SENSITIVITY_POINT_ESTIMATES.tsv", sep="\t"
    )
    secondary = pd.read_csv(root / "results_locked" / "V1_12_SECONDARY_METRICS.tsv", sep="\t")
    covariance = pd.read_csv(
        root / "results_locked" / "V1_12_RANK_COVARIANCE_DECOMPOSITION.tsv", sep="\t"
    )
    resampling = pd.read_csv(
        root
        / "results_locked"
        / "extended_sensitivity"
        / "V1_12_RESAMPLING_UNIT_SENSITIVITY.tsv",
        sep="\t",
    )
    source_intervals = pd.read_csv(
        root
        / "results_locked"
        / "extended_sensitivity"
        / "V1_12_SOURCE_ORIGIN_BOOTSTRAP_INTERVALS.tsv",
        sep="\t",
    )
    blind = pd.read_csv(root / "source_frozen" / "V1_12_BLIND_STRUCTURES.tsv", sep="\t", dtype=str)
    target_support = pd.read_csv(
        root / "source_frozen" / "V1_12_TARGET_SUPPORT.tsv", sep="\t"
    )
    if set(primary["model_id"]) != set(PRIMARY_MODELS):
        raise RuntimeError("Primary figure input omits or adds a frozen model")
    observed_source_sensitivities = set(
        sensitivities.loc[
            sensitivities["category"].isin(
                {"mandatory_source_origin", "mandatory_leave_one_origin_out"}
            ),
            "sensitivity",
        ]
    )
    expected = set()
    for source in blind["source_origin_stratum"].dropna().unique():
        expected.add(f"source_origin__{source}")
        expected.add(f"leave_one_origin_out__{source}")
    if not expected.issubset(observed_source_sensitivities):
        raise RuntimeError("Mandatory source-origin or leave-one-origin-out rows are missing")
    configure_style()
    for directory in (
        root / "figures" / "main",
        root / "figures" / "supplementary",
        root / "figure_source_data",
        root / "qc",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    source_data_dir = root / "figure_source_data"
    outputs = []
    outputs.extend(build_main(root, manifest, primary, sensitivities, source_data_dir))
    outputs.extend(
        build_supplement(
            root,
            manifest,
            blind,
            target_support,
            secondary,
            covariance,
            sensitivities,
            source_data_dir,
        )
    )
    outputs.extend(
        build_resampling_supplement(
            root,
            primary,
            sensitivities,
            resampling,
            source_intervals,
            source_data_dir,
        )
    )
    receipt = {
        "schema_version": "science_advances_v1_12_figure_build_receipt_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "FIGURES_BUILT_FROM_COMPLETE_LOCKED_RESULTS",
        "display_policy": "All four frozen families, all observed source origins, and every leave-one-origin-out row are displayed without directional filtering.",
        "figure_code_sha256": sha256_file(Path(__file__).resolve()),
        "analysis_complete_sha256": sha256_file(complete_path),
        "outputs": {
            str(path.relative_to(root)): {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in outputs
        },
    }
    receipt_path = root / "qc" / "V1_12_FIGURE_BUILD_RECEIPT.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "files": len(outputs)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
