#!/usr/bin/env python3
"""Run the precommitted v1.12 analysis after the Stage 2 firewall opens."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import stat
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import rankdata


ALL_MODELS = [
    "development_exact_median_v1_8",
    "global_top5_morgan_similarity_v1_8",
    "hgb_fusion_additive_pre2018_v1_8",
    "ligand_mpnn_frozen_esm2_v1_8",
    "smiles_protein_cross_attention_v1_8",
]
PRIMARY_MODELS = ALL_MODELS[1:]
PRIMARY_ENDPOINTS = [
    "source_wide_spearman",
    "within_target_centered_spearman",
    "within_target_pairwise_concordance",
    "target_rank_attenuation",
]
ENDPOINT_NULLS = {
    "source_wide_spearman": 0.0,
    "within_target_centered_spearman": 0.0,
    "within_target_pairwise_concordance": 0.5,
    "target_rank_attenuation": 0.0,
}
SOURCE_STRATA = [
    "strict_dual_timestamp_patent",
    "database_entry_patent",
    "database_entry_residual_patent",
    "database_entry_chembl_import",
    "database_entry_pubchem_assay",
    "database_entry_curated_article",
    "database_entry_other_documented",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected a JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def finite_or_none(value: float) -> float | None:
    return float(value) if np.isfinite(value) else None


def parse_bool(series: pd.Series) -> np.ndarray:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes"})
        .to_numpy(bool)
    )


def weighted_rank(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Return frequency-weighted midranks; zero-weight observations are harmless."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    sorted_weights = weights[order]
    starts = np.r_[0, np.flatnonzero(sorted_values[1:] != sorted_values[:-1]) + 1]
    group_weights = np.add.reduceat(sorted_weights, starts)
    before = np.cumsum(group_weights) - group_weights
    group_ranks = before + 0.5 * group_weights
    lengths = np.diff(np.r_[starts, len(values)])
    ranked_sorted = np.repeat(group_ranks, lengths)
    ranked = np.empty(len(values), dtype=float)
    ranked[order] = ranked_sorted
    return ranked


def weighted_corr(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
    if valid.sum() < 3:
        return math.nan
    x = x[valid]
    y = y[valid]
    weights = weights[valid]
    total = weights.sum()
    x_centered = x - np.dot(weights, x) / total
    y_centered = y - np.dot(weights, y) / total
    covariance = np.dot(weights, x_centered * y_centered)
    denominator = math.sqrt(
        float(np.dot(weights, x_centered * x_centered))
        * float(np.dot(weights, y_centered * y_centered))
    )
    return float(covariance / denominator) if denominator > 0 else math.nan


def weighted_spearman(
    x: np.ndarray, y: np.ndarray, weights: np.ndarray
) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
    if valid.sum() < 3:
        return math.nan
    x = x[valid]
    y = y[valid]
    weights = weights[valid]
    return weighted_corr(
        weighted_rank(x, weights), weighted_rank(y, weights), weights
    )


def weighted_group_center(
    values: np.ndarray, group_codes: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    groups = int(group_codes.max()) + 1 if len(group_codes) else 0
    totals = np.bincount(group_codes, weights=weights, minlength=groups)
    sums = np.bincount(
        group_codes, weights=weights * values, minlength=groups
    )
    means = np.divide(
        sums,
        totals,
        out=np.full(groups, np.nan, dtype=float),
        where=totals > 0,
    )
    return values - means[group_codes]


def centered_spearman(
    observed: np.ndarray,
    predicted: np.ndarray,
    group_codes: np.ndarray,
    eligible_rows: np.ndarray,
    weights: np.ndarray,
) -> float:
    mask = eligible_rows & (weights > 0)
    if mask.sum() < 3:
        return math.nan
    local_groups, _ = pd.factorize(group_codes[mask], sort=True)
    local_weights = weights[mask]
    observed_centered = weighted_group_center(
        observed[mask], local_groups, local_weights
    )
    predicted_centered = weighted_group_center(
        predicted[mask], local_groups, local_weights
    )
    return weighted_spearman(observed_centered, predicted_centered, local_weights)


@dataclass(frozen=True)
class ConcordanceSpec:
    target_code: int
    component_codes: np.ndarray
    component_row_counts: np.ndarray
    signed_pair_matrix: np.ndarray


@dataclass
class ModelData:
    model_id: str
    pair_ids: np.ndarray
    observed: np.ndarray
    predicted: np.ndarray
    interval_low: np.ndarray
    interval_high: np.ndarray
    protein_labels: np.ndarray
    target_codes: np.ndarray
    scaffold_labels: np.ndarray
    component_codes: np.ndarray
    source_labels: np.ndarray
    eligible_target_codes: np.ndarray
    eligible_rows: np.ndarray
    concordance_specs: list[ConcordanceSpec]


def aggregate_signed_pairs(
    observed: np.ndarray,
    predicted: np.ndarray,
    local_component_codes: np.ndarray,
    component_count: int,
    chunk_size: int = 384,
) -> np.ndarray:
    matrix = np.zeros((component_count, component_count), dtype=float)
    n_rows = len(observed)
    column_components = local_component_codes
    for start in range(0, n_rows, chunk_size):
        stop = min(start + chunk_size, n_rows)
        signed = np.sign(observed[start:stop, None] - observed[None, :])
        signed *= np.sign(predicted[start:stop, None] - predicted[None, :])
        row_components = np.repeat(
            local_component_codes[start:stop], n_rows
        )
        tiled_columns = np.tile(column_components, stop - start)
        np.add.at(matrix, (row_components, tiled_columns), signed.ravel())
    return matrix


def build_concordance_specs(data: ModelData) -> list[ConcordanceSpec]:
    specs: list[ConcordanceSpec] = []
    for target_code in data.eligible_target_codes:
        indices = np.flatnonzero(data.target_codes == target_code)
        components, local_codes = np.unique(
            data.component_codes[indices], return_inverse=True
        )
        component_rows = np.bincount(
            local_codes, minlength=len(components)
        ).astype(float)
        signed_matrix = aggregate_signed_pairs(
            data.observed[indices],
            data.predicted[indices],
            local_codes,
            len(components),
        )
        specs.append(
            ConcordanceSpec(
                target_code=int(target_code),
                component_codes=components.astype(int),
                component_row_counts=component_rows,
                signed_pair_matrix=signed_matrix,
            )
        )
    return specs


def concordance_from_component_multiplicities(
    specs: list[ConcordanceSpec], component_multiplicities: np.ndarray
) -> np.ndarray:
    multiplicities = np.asarray(component_multiplicities, dtype=float)
    if multiplicities.ndim == 1:
        multiplicities = multiplicities[None, :]
    sums = np.zeros(multiplicities.shape[0], dtype=float)
    counts = np.zeros(multiplicities.shape[0], dtype=int)
    for spec in specs:
        local = multiplicities[:, spec.component_codes]
        signed_sum = 0.5 * np.sum(
            (local @ spec.signed_pair_matrix) * local, axis=1
        )
        expanded_rows = local @ spec.component_row_counts
        squared_row_weights = (local * local) @ spec.component_row_counts
        denominator = 0.5 * (expanded_rows * expanded_rows - squared_row_weights)
        valid = denominator > 0
        values = np.full(len(denominator), np.nan, dtype=float)
        values[valid] = 0.5 + 0.5 * signed_sum[valid] / denominator[valid]
        sums[valid] += values[valid]
        counts[valid] += 1
    return np.divide(
        sums,
        counts,
        out=np.full(len(sums), np.nan, dtype=float),
        where=counts > 0,
    )


def prepare_model_data(
    frame: pd.DataFrame,
    model_id: str,
    global_component_map: dict[str, int],
) -> ModelData:
    subset = frame.loc[frame["model_id"].eq(model_id)].copy()
    abstained = parse_bool(subset["abstained"])
    subset["predicted_pkd"] = pd.to_numeric(
        subset["predicted_pkd"], errors="coerce"
    )
    valid = (~abstained) & np.isfinite(subset["predicted_pkd"].to_numpy(float))
    subset = subset.loc[valid].reset_index(drop=True)
    if len(subset) < 3:
        raise RuntimeError(f"Too few nonabstained predictions for {model_id}")
    protein_codes, protein_labels = pd.factorize(
        subset["protein_sha256"], sort=True
    )
    support = (
        subset.assign(target_code=protein_codes)
        .groupby("target_code", sort=True)
        .agg(
            ligands=("canonical_smiles", "nunique"),
            scaffolds=("scaffold_sha256", "nunique"),
        )
    )
    eligible_codes = support.index[
        support["ligands"].ge(5) & support["scaffolds"].ge(3)
    ].to_numpy(int)
    data = ModelData(
        model_id=model_id,
        pair_ids=subset["blind_pair_id"].astype(str).to_numpy(),
        observed=pd.to_numeric(subset["observed_pkd"], errors="raise").to_numpy(float),
        predicted=subset["predicted_pkd"].to_numpy(float),
        interval_low=pd.to_numeric(subset["interval_low"], errors="coerce").to_numpy(float),
        interval_high=pd.to_numeric(subset["interval_high"], errors="coerce").to_numpy(float),
        protein_labels=subset["protein_sha256"].astype(str).to_numpy(),
        target_codes=protein_codes.astype(int),
        scaffold_labels=subset["scaffold_sha256"].astype(str).to_numpy(),
        component_codes=subset["document_component_id"]
        .map(global_component_map)
        .to_numpy(int),
        source_labels=subset["source_origin_stratum"].astype(str).to_numpy(),
        eligible_target_codes=eligible_codes,
        eligible_rows=np.isin(protein_codes, eligible_codes),
        concordance_specs=[],
    )
    data.concordance_specs = build_concordance_specs(data)
    return data


def point_core_metrics(data: ModelData, component_count: int) -> dict[str, float]:
    row_weights = np.ones(len(data.observed), dtype=float)
    source = weighted_spearman(data.observed, data.predicted, row_weights)
    centered = centered_spearman(
        data.observed,
        data.predicted,
        data.target_codes,
        data.eligible_rows,
        row_weights,
    )
    concordance = concordance_from_component_multiplicities(
        data.concordance_specs, np.ones(component_count, dtype=float)
    )[0]
    return {
        "source_wide_spearman": source,
        "within_target_centered_spearman": centered,
        "within_target_pairwise_concordance": concordance,
        "target_rank_attenuation": source - centered,
    }


def bootstrap_core_metrics(
    data: ModelData, component_draws: np.ndarray
) -> dict[str, np.ndarray]:
    n_replicates = component_draws.shape[0]
    source = np.full(n_replicates, np.nan, dtype=float)
    centered = np.full(n_replicates, np.nan, dtype=float)
    for replicate in range(n_replicates):
        weights = component_draws[replicate, data.component_codes].astype(float)
        source[replicate] = weighted_spearman(
            data.observed, data.predicted, weights
        )
        centered[replicate] = centered_spearman(
            data.observed,
            data.predicted,
            data.target_codes,
            data.eligible_rows,
            weights,
        )
    concordance = concordance_from_component_multiplicities(
        data.concordance_specs, component_draws
    )
    return {
        "source_wide_spearman": source,
        "within_target_centered_spearman": centered,
        "within_target_pairwise_concordance": concordance,
        "target_rank_attenuation": source - centered,
    }


def holm_adjust(values: Iterable[float]) -> list[float]:
    p_values = np.asarray(list(values), dtype=float)
    adjusted = np.full(len(p_values), np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(p_values))
    if len(finite) == 0:
        return adjusted.tolist()
    ordered = finite[np.argsort(p_values[finite])]
    running = 0.0
    m = len(ordered)
    for position, index in enumerate(ordered):
        candidate = min(1.0, (m - position) * p_values[index])
        running = max(running, candidate)
        adjusted[index] = running
    return adjusted.tolist()


def simultaneous_primary_inference(
    points: dict[str, dict[str, float]],
    bootstraps: dict[str, dict[str, np.ndarray]],
    minimum_valid: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    records: list[dict[str, Any]] = []
    endpoint_meta: dict[str, Any] = {}
    p_lookup: dict[tuple[str, str], float] = {}
    for endpoint in PRIMARY_ENDPOINTS:
        matrix = np.column_stack(
            [bootstraps[model][endpoint] for model in PRIMARY_MODELS]
        )
        point_vector = np.asarray(
            [points[model][endpoint] for model in PRIMARY_MODELS], dtype=float
        )
        standard_errors = np.nanstd(matrix, axis=0, ddof=1)
        complete = np.all(np.isfinite(matrix), axis=1)
        complete &= np.all(np.isfinite(standard_errors) & (standard_errors > 0))
        valid_count = int(complete.sum())
        critical = math.nan
        max_statistics = np.asarray([], dtype=float)
        if valid_count >= minimum_valid:
            studentized = (
                matrix[complete] - point_vector[None, :]
            ) / standard_errors[None, :]
            max_statistics = np.max(np.abs(studentized), axis=1)
            critical = float(np.quantile(max_statistics, 0.95, method="higher"))
        endpoint_meta[endpoint] = {
            "valid_joint_replicates": valid_count,
            "minimum_required": minimum_valid,
            "simultaneous_critical_value": finite_or_none(critical),
        }
        for model_index, model_id in enumerate(PRIMARY_MODELS):
            estimate = point_vector[model_index]
            standard_error = standard_errors[model_index]
            lower = estimate - critical * standard_error
            upper = estimate + critical * standard_error
            observed_statistic = abs(
                (estimate - ENDPOINT_NULLS[endpoint]) / standard_error
            )
            p_value = math.nan
            if len(max_statistics):
                p_value = float(
                    (1 + np.sum(max_statistics >= observed_statistic))
                    / (len(max_statistics) + 1)
                )
            p_lookup[(model_id, endpoint)] = p_value
            records.append(
                {
                    "model_id": model_id,
                    "endpoint": endpoint,
                    "null_value": ENDPOINT_NULLS[endpoint],
                    "estimate": estimate,
                    "bootstrap_se": standard_error,
                    "simultaneous_95_ci_low": lower,
                    "simultaneous_95_ci_high": upper,
                    "max_t_adjusted_p": p_value,
                    "valid_joint_bootstrap_replicates": valid_count,
                }
            )
    for model_id in PRIMARY_MODELS:
        adjusted = holm_adjust(
            p_lookup[(model_id, endpoint)] for endpoint in PRIMARY_ENDPOINTS
        )
        for endpoint, value in zip(PRIMARY_ENDPOINTS, adjusted):
            for record in records:
                if record["model_id"] == model_id and record["endpoint"] == endpoint:
                    record["holm_across_endpoints_p"] = value
                    break
    return pd.DataFrame(records), endpoint_meta


def decision_from_primary(primary: pd.DataFrame) -> dict[str, Any]:
    by_key = primary.set_index(["model_id", "endpoint"])
    family_rows: list[dict[str, Any]] = []
    for model_id in PRIMARY_MODELS:
        source = by_key.loc[(model_id, "source_wide_spearman")]
        centered = by_key.loc[(model_id, "within_target_centered_spearman")]
        concordance = by_key.loc[
            (model_id, "within_target_pairwise_concordance")
        ]
        attenuation = by_key.loc[(model_id, "target_rank_attenuation")]
        source_support = bool(source["simultaneous_95_ci_low"] > 0)
        centered_support = bool(centered["simultaneous_95_ci_low"] > 0)
        concordance_support = bool(concordance["simultaneous_95_ci_low"] > 0.5)
        material_attenuation = bool(
            attenuation["estimate"] >= 0.10
            and attenuation["simultaneous_95_ci_low"] > 0
        )
        centered_lower = bool(centered["estimate"] < source["estimate"])
        family_rows.append(
            {
                "model_id": model_id,
                "source_wide_support": source_support,
                "within_target_centered_support": centered_support,
                "within_target_concordance_support": concordance_support,
                "material_attenuation": material_attenuation,
                "within_target_centered_lower_than_source_wide": centered_lower,
                "family_structural_support": source_support
                and material_attenuation,
            }
        )
    source_count = sum(row["source_wide_support"] for row in family_rows)
    centered_count = sum(
        row["within_target_centered_support"] for row in family_rows
    )
    attenuation_count = sum(row["material_attenuation"] for row in family_rows)
    structural_count = sum(row["family_structural_support"] for row in family_rows)
    all_centered_lower = all(
        row["within_target_centered_lower_than_source_wide"]
        for row in family_rows
    )
    if structural_count >= 3 and all_centered_lower:
        branch = "composition_dominated_replication"
    elif source_count <= 1:
        branch = "no_temporal_transport"
    elif (
        source_count >= 3
        and centered_count >= 3
        and attenuation_count < 3
        and not all_centered_lower
    ):
        branch = "transported_within_target_ordering"
    else:
        branch = "heterogeneous_or_inconclusive"
    return {
        "schema_version": "science_advances_v1_12_locked_decision_v1",
        "decision_branch": branch,
        "family_results": family_rows,
        "counts": {
            "source_wide_support": source_count,
            "within_target_centered_support": centered_count,
            "material_attenuation": attenuation_count,
            "family_structural_support": structural_count,
        },
        "within_target_centered_lower_in_all_four": all_centered_lower,
        "claim_boundary": (
            "Retrospective outcome-unexposed, structure-adapted temporal aggregate validation assembled after "
            "structure-only feasibility review; not calendar-prospective, source-independent, "
            "mechanistic, therapeutic, clinical, or wet-lab validation."
        ),
    }


def grouped_centered_spearman(
    observed: np.ndarray, predicted: np.ndarray, labels: np.ndarray
) -> float:
    codes, _ = pd.factorize(labels, sort=True)
    counts = np.bincount(codes)
    eligible = counts[codes] >= 2
    weights = np.ones(len(observed), dtype=float)
    return centered_spearman(
        observed, predicted, codes, eligible, weights
    )


def target_utility(data: ModelData) -> tuple[float, float, int]:
    enrichments: list[float] = []
    regrets: list[float] = []
    for target_code in data.eligible_target_codes:
        indices = np.flatnonzero(data.target_codes == target_code)
        n_rows = len(indices)
        if n_rows < 2:
            continue
        top_n = max(1, int(math.ceil(0.10 * n_rows)))
        observed_order = indices[np.argsort(data.observed[indices], kind="stable")]
        predicted_order = indices[np.argsort(data.predicted[indices], kind="stable")]
        observed_top = set(observed_order[-top_n:])
        predicted_top = set(predicted_order[-top_n:])
        overlap = len(observed_top & predicted_top)
        enrichments.append(float(overlap * n_rows / (top_n * top_n)))
        selected = predicted_order[-1]
        regrets.append(
            float(np.max(data.observed[indices]) - data.observed[selected])
        )
    if not enrichments:
        return math.nan, math.nan, 0
    return float(np.mean(enrichments)), float(np.mean(regrets)), len(enrichments)


def secondary_metrics(
    all_frame: pd.DataFrame,
    model_data: dict[str, ModelData],
    component_count: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    covariance_rows: list[dict[str, Any]] = []
    for model_id in ALL_MODELS:
        data = model_data[model_id]
        metrics = point_core_metrics(data, component_count)
        residual = data.observed - data.predicted
        prediction_variance = float(np.var(data.predicted, ddof=1))
        slope = (
            float(np.cov(data.predicted, data.observed, ddof=1)[0, 1])
            / prediction_variance
            if prediction_variance > 0
            else math.nan
        )
        intercept = (
            float(np.mean(data.observed) - slope * np.mean(data.predicted))
            if np.isfinite(slope)
            else math.nan
        )
        interval_valid = np.isfinite(data.interval_low) & np.isfinite(data.interval_high)
        coverage = (
            float(
                np.mean(
                    (data.observed[interval_valid] >= data.interval_low[interval_valid])
                    & (data.observed[interval_valid] <= data.interval_high[interval_valid])
                )
            )
            if interval_valid.any()
            else math.nan
        )
        interval_width = (
            float(np.mean(data.interval_high[interval_valid] - data.interval_low[interval_valid]))
            if interval_valid.any()
            else math.nan
        )
        enrichment, regret, utility_targets = target_utility(data)
        model_rows = all_frame.loc[all_frame["model_id"].eq(model_id)]
        abstention_count = int(parse_bool(model_rows["abstained"]).sum())
        rows.append(
            {
                "model_id": model_id,
                "nonabstained_pairs": len(data.observed),
                "abstentions": abstention_count,
                "abstention_rate": abstention_count / len(model_rows),
                "eligible_within_target_groups": len(data.eligible_target_codes),
                **metrics,
                "within_document_component_centered_spearman": grouped_centered_spearman(
                    data.observed,
                    data.predicted,
                    np.asarray(
                        [str(value) for value in data.component_codes], dtype=object
                    ),
                ),
                "within_source_stratum_centered_spearman": grouped_centered_spearman(
                    data.observed, data.predicted, data.source_labels
                ),
                "mae": float(np.mean(np.abs(residual))),
                "median_absolute_error": float(np.median(np.abs(residual))),
                "rmse": float(np.sqrt(np.mean(residual * residual))),
                "calibration_slope_observed_on_predicted": slope,
                "calibration_intercept": intercept,
                "dispersion_ratio_predicted_to_observed": float(
                    np.std(data.predicted, ddof=1) / np.std(data.observed, ddof=1)
                ),
                "interval_90_coverage": coverage,
                "mean_interval_width": interval_width,
                "within_target_top10_enrichment": enrichment,
                "best_ligand_selection_regret_pkd": regret,
                "utility_target_groups": utility_targets,
            }
        )
        observed_ranks = rankdata(data.observed, method="average")
        predicted_ranks = rankdata(data.predicted, method="average")
        for group_name, labels in (
            ("exact_protein", data.protein_labels),
            ("document_component", data.component_codes.astype(str)),
            ("source_stratum", data.source_labels),
        ):
            codes, _ = pd.factorize(labels, sort=True)
            counts = np.bincount(codes).astype(float)
            mean_observed = np.bincount(
                codes, weights=observed_ranks
            ) / counts
            mean_predicted = np.bincount(
                codes, weights=predicted_ranks
            ) / counts
            overall_observed = float(np.mean(observed_ranks))
            overall_predicted = float(np.mean(predicted_ranks))
            between = float(
                np.sum(
                    counts
                    * (mean_observed - overall_observed)
                    * (mean_predicted - overall_predicted)
                )
                / len(codes)
            )
            within = float(
                np.mean(
                    (observed_ranks - mean_observed[codes])
                    * (predicted_ranks - mean_predicted[codes])
                )
            )
            total = float(
                np.mean(
                    (observed_ranks - overall_observed)
                    * (predicted_ranks - overall_predicted)
                )
            )
            covariance_rows.append(
                {
                    "model_id": model_id,
                    "grouping": group_name,
                    "groups": len(counts),
                    "total_rank_covariance": total,
                    "between_group_rank_covariance": between,
                    "within_group_rank_covariance": within,
                    "decomposition_error": total - between - within,
                    "between_fraction_of_total": between / total
                    if total != 0
                    else math.nan,
                }
            )
    comparison_rows: list[dict[str, Any]] = []
    references = [ALL_MODELS[0], ALL_MODELS[1]]
    indexed = {
        model_id: pd.Series(data.predicted, index=data.pair_ids)
        for model_id, data in model_data.items()
    }
    observed_indexed = pd.Series(
        model_data[ALL_MODELS[0]].observed,
        index=model_data[ALL_MODELS[0]].pair_ids,
    )
    for model_id in ALL_MODELS:
        for reference in references:
            common = indexed[model_id].index.intersection(indexed[reference].index)
            common = common.intersection(observed_indexed.index)
            observed = observed_indexed.loc[common].to_numpy(float)
            model_error = np.abs(observed - indexed[model_id].loc[common].to_numpy(float))
            reference_error = np.abs(
                observed - indexed[reference].loc[common].to_numpy(float)
            )
            comparison_rows.append(
                {
                    "model_id": model_id,
                    "reference_model_id": reference,
                    "common_pairs": len(common),
                    "mae_difference_model_minus_reference": float(
                        np.mean(model_error - reference_error)
                    ),
                    "win_rate_lower_absolute_error": float(
                        np.mean(model_error < reference_error)
                    ),
                }
            )
    return (
        pd.DataFrame(rows),
        pd.DataFrame(covariance_rows),
        pd.DataFrame(comparison_rows),
    )


def sensitivity_masks(base: pd.DataFrame) -> list[tuple[str, str, np.ndarray]]:
    publication = pd.to_datetime(base["publication_date_min"], errors="coerce")
    entry = pd.to_datetime(base["bindingdb_date_min"], errors="coerce")
    cutoff = pd.Timestamp("2024-01-02")
    measurement_count = pd.to_numeric(base["measurement_count"], errors="coerce")
    exact_10000 = parse_bool(base["is_exact_10000_nm"])
    component_rank = pd.to_numeric(
        base["component_selection_rank"], errors="coerce"
    )
    masks: list[tuple[str, str, np.ndarray]] = [
        (
            "strict_publication_and_entry_post_2024",
            "locked_temporal",
            (publication >= cutoff).to_numpy(bool) & (entry >= cutoff).to_numpy(bool),
        ),
        ("exclude_exact_10000_nm", "locked_endpoint", ~exact_10000),
        (
            "single_measurement_pairs",
            "locked_measurement",
            measurement_count.eq(1).to_numpy(bool),
        ),
        (
            "replicated_pairs",
            "locked_measurement",
            measurement_count.gt(1).to_numpy(bool),
        ),
    ]
    observed_strata = [
        stratum
        for stratum in SOURCE_STRATA
        if base["source_origin_stratum"].eq(stratum).any()
    ]
    for stratum in observed_strata:
        masks.append(
            (
                f"source_origin__{stratum}",
                "mandatory_source_origin",
                base["source_origin_stratum"].eq(stratum).to_numpy(bool),
            )
        )
        masks.append(
            (
                f"leave_one_origin_out__{stratum}",
                "mandatory_leave_one_origin_out",
                base["source_origin_stratum"].ne(stratum).to_numpy(bool),
            )
        )
    for cap in (10, 15, 20):
        masks.append(
            (
                f"nested_document_cap__{cap}",
                "locked_document_balance",
                component_rank.le(cap).fillna(False).to_numpy(bool),
            )
        )
    lag = (entry - publication).dt.days
    masks.extend(
        [
            (
                "publication_to_entry_lag__0_90_days",
                "locked_lag",
                lag.between(0, 90, inclusive="both").fillna(False).to_numpy(bool),
            ),
            (
                "publication_to_entry_lag__91_365_days",
                "locked_lag",
                lag.between(91, 365, inclusive="both").fillna(False).to_numpy(bool),
            ),
            (
                "publication_to_entry_lag__over_365_days",
                "locked_lag",
                lag.gt(365).fillna(False).to_numpy(bool),
            ),
            (
                "publication_to_entry_lag__unknown_or_negative",
                "locked_lag",
                (lag.isna() | lag.lt(0)).to_numpy(bool),
            ),
        ]
    )
    return masks


def point_sensitivity_analysis(
    base: pd.DataFrame,
    predictions: pd.DataFrame,
    global_component_map: dict[str, int],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    base_indexed = base.set_index("blind_pair_id", drop=False)
    component_count = len(global_component_map)
    for sensitivity, category, mask in sensitivity_masks(base):
        pair_ids = set(base.loc[mask, "blind_pair_id"].astype(str))
        if len(pair_ids) == 0:
            for model_id in PRIMARY_MODELS:
                rows.append(
                    {
                        "sensitivity": sensitivity,
                        "category": category,
                        "model_id": model_id,
                        "nonabstained_pairs": 0,
                        "exact_proteins": 0,
                        "eligible_within_target_groups": 0,
                        **{endpoint: math.nan for endpoint in PRIMARY_ENDPOINTS},
                    }
                )
            continue
        local_base = base_indexed.loc[sorted(pair_ids)].reset_index(drop=True)
        local_predictions = predictions.loc[
            predictions["blind_pair_id"].astype(str).isin(pair_ids)
        ]
        local_long = local_predictions.merge(
            local_base, on="blind_pair_id", how="inner", validate="many_to_one"
        )
        for model_id in PRIMARY_MODELS:
            try:
                data = prepare_model_data(
                    local_long, model_id, global_component_map
                )
                metrics = point_core_metrics(data, component_count)
                rows.append(
                    {
                        "sensitivity": sensitivity,
                        "category": category,
                        "model_id": model_id,
                        "nonabstained_pairs": len(data.observed),
                        "exact_proteins": len(np.unique(data.protein_labels)),
                        "eligible_within_target_groups": len(
                            data.eligible_target_codes
                        ),
                        **metrics,
                    }
                )
            except RuntimeError:
                rows.append(
                    {
                        "sensitivity": sensitivity,
                        "category": category,
                        "model_id": model_id,
                        "nonabstained_pairs": 0,
                        "exact_proteins": 0,
                        "eligible_within_target_groups": 0,
                        **{endpoint: math.nan for endpoint in PRIMARY_ENDPOINTS},
                    }
                )
    return pd.DataFrame(rows)


def validate_and_unblind(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    paths = {
        "protocol": root
        / "protocol"
        / "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json",
        "implementation_spec": root / "protocol" / "V1_12_ANALYSIS_IMPLEMENTATION_SPEC.json",
        "source_manifest": root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json",
        "blind": root / "source_frozen" / "V1_12_BLIND_STRUCTURES.tsv",
        "join_map": root / "source_frozen" / "V1_12_OUTCOME_JOIN_MAP.tsv",
        "predictions": root / "predictions_frozen" / "V1_12_FROZEN_PREDICTIONS.tsv",
        "prediction_receipt": root / "predictions_frozen" / "V1_12_PREDICTIONS_FROZEN_RECEIPT.json",
        "stage2_lock": root / "locks" / "V1_12_STAGE2_LOCK.json",
        "analysis_code_lock": root / "locks" / "V1_12_ANALYSIS_CODE_LOCK.json",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise RuntimeError(f"Required locked inputs are missing: {missing}")
    complete_path = root / "results_locked" / "V1_12_ANALYSIS_COMPLETE.json"
    if complete_path.exists():
        raise RuntimeError("Locked analysis already completed; rerun is prohibited")
    if (root / "locks" / "V1_12_UNBLIND_RECEIPT.json").exists():
        raise RuntimeError("An unblind receipt already exists; automatic rerun is prohibited")
    protocol = load_json(paths["protocol"])
    implementation = load_json(paths["implementation_spec"])
    manifest = load_json(paths["source_manifest"])
    prediction_receipt = load_json(paths["prediction_receipt"])
    stage2 = load_json(paths["stage2_lock"])
    code_lock = load_json(paths["analysis_code_lock"])
    manifest_quarantines = {
        str(item["cohort"]): str(item["precommitted_sha256"])
        for item in manifest.get("quarantine_sources", [])
    }
    stage2_quarantines = {
        str(item["cohort"]): str(item["precommitted_sha256"])
        for item in stage2.get("quarantined_outcomes", [])
    }
    checks = {
        "protocol_status_locked": protocol.get("status")
        == "OUTCOME_BLIND_ANALYSIS_LOCKED_AFTER_STRUCTURE_ONLY_FEASIBILITY",
        "implementation_preoutcome": implementation.get("status")
        == "LOCKED_BEFORE_OUTCOME_OR_PREDICTION_ACCESS",
        "source_gate_passed": manifest.get("status")
        == "STRUCTURE_ONLY_MEMBERSHIP_FROZEN",
        "calendar_prospective_false": manifest.get("calendar_prospective") is False,
        "stage2_all_checks_passed": stage2.get("all_checks_passed") is True,
        "stage2_authorized": stage2.get("outcome_access_authorized") is True,
        "analysis_code_hash_matches": code_lock.get("analysis_code_sha256")
        == sha256_file(Path(__file__).resolve()),
        "protocol_hash_matches": code_lock.get("stage1_protocol_sha256")
        == sha256_file(paths["protocol"]),
        "implementation_hash_matches": code_lock.get("implementation_spec_sha256")
        == sha256_file(paths["implementation_spec"]),
        "blind_hash_matches": manifest.get("blind_structure_table_sha256")
        == sha256_file(paths["blind"]),
        "prediction_hash_matches": prediction_receipt.get("prediction_file_sha256")
        == sha256_file(paths["predictions"]),
        "prediction_blind_hash_matches": prediction_receipt.get(
            "blind_structure_table_sha256"
        )
        == sha256_file(paths["blind"]),
        "join_map_hash_matches": manifest.get("outcome_join_map_sha256")
        == sha256_file(paths["join_map"]),
        "quarantine_hashes_precommitted_consistently": manifest_quarantines
        == stage2_quarantines
        and set(manifest_quarantines) == {"v1.9", "v1.10", "v1.11"},
    }
    quarantine_paths: dict[str, Path] = {}
    restricted_modes: dict[str, str] = {}
    for item in manifest.get("quarantine_sources", []):
        cohort = str(item["cohort"])
        outcome_path = Path(str(item["path"]))
        quarantine_paths[cohort] = outcome_path
        if not outcome_path.is_file():
            checks[f"{cohort}_outcome_file_exists"] = False
            continue
        mode = stat.S_IMODE(os.stat(outcome_path).st_mode)
        restricted_modes[cohort] = format(mode, "03o")
        checks[f"{cohort}_outcome_restricted_before_unblind"] = (
            mode & 0o377
        ) == 0
    if not all(checks.values()):
        failures = [key for key, value in checks.items() if not value]
        raise RuntimeError(f"Pre-unblind firewall failed: {failures}")

    blind = pd.read_csv(paths["blind"], sep="\t", dtype=str, keep_default_na=False)
    predictions = pd.read_csv(
        paths["predictions"], sep="\t", dtype=str, keep_default_na=False
    )
    if set(predictions["model_id"]) != set(ALL_MODELS):
        raise RuntimeError("Frozen prediction model membership changed")
    if len(predictions) != len(blind) * len(ALL_MODELS):
        raise RuntimeError("Frozen prediction row count changed")
    if predictions.duplicated(["blind_pair_id", "model_id"]).any():
        raise RuntimeError("Duplicate frozen pair-model prediction rows")
    join_map = pd.read_csv(
        paths["join_map"], sep="\t", dtype=str, keep_default_na=False
    )
    if join_map.duplicated(["origin_cohort", "origin_blind_pair_id"]).any():
        raise RuntimeError("Duplicate source identifiers in the frozen outcome join map")
    if set(join_map["blind_pair_id"]) != set(blind["blind_pair_id"]):
        raise RuntimeError("Join map and blind pair membership differ")

    observed_hashes: dict[str, str] = {}
    for cohort, outcome_path in quarantine_paths.items():
        os.chmod(outcome_path, 0o600)
        observed_hash = sha256_file(outcome_path)
        observed_hashes[cohort] = observed_hash
        if observed_hash != manifest_quarantines[cohort]:
            os.chmod(outcome_path, 0)
            raise RuntimeError(
                f"{cohort} quarantined outcome hash changed at authorized unblind"
            )
    unblind_receipt = {
        "schema_version": "science_advances_v1_12_unblind_receipt_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "MULTISOURCE_OUTCOME_ACCESS_AUTHORIZED_HASH_VERIFIED",
        "pre_unblind_checks": checks,
        "outcome_sha256": observed_hashes,
        "outcome_file_modes_before_authorization": restricted_modes,
        "outcome_file_mode_after_authorization": "600",
        "stage2_lock_sha256": sha256_file(paths["stage2_lock"]),
        "analysis_code_lock_sha256": sha256_file(paths["analysis_code_lock"]),
        "analysis_code_sha256": sha256_file(Path(__file__).resolve()),
        "outcome_values_parsed_before_receipt": False,
        "prediction_rerun_allowed": False,
    }
    write_json(root / "locks" / "V1_12_UNBLIND_RECEIPT.json", unblind_receipt)
    source_outcomes: list[pd.DataFrame] = []
    source_rows: dict[str, int] = {}
    for cohort, outcome_path in quarantine_paths.items():
        local = pd.read_csv(
            outcome_path, sep="\t", dtype=str, keep_default_na=False
        )
        local = local.rename(columns={"blind_pair_id": "origin_blind_pair_id"})
        local.insert(0, "origin_cohort", cohort)
        source_rows[cohort] = len(local)
        source_outcomes.append(local)
    all_outcomes = pd.concat(source_outcomes, ignore_index=True)
    required_outcomes = {
        "origin_cohort",
        "origin_blind_pair_id",
        "observed_pkd",
        "kd_nm_median",
        "measurement_count",
        "is_exact_10000_nm",
    }
    if not required_outcomes.issubset(all_outcomes.columns):
        raise RuntimeError("Quarantined outcome schema changed")
    if all_outcomes.duplicated(["origin_cohort", "origin_blind_pair_id"]).any():
        raise RuntimeError("Duplicate source outcome pair identifiers")
    outcomes = join_map.merge(
        all_outcomes,
        on=["origin_cohort", "origin_blind_pair_id"],
        how="left",
        validate="one_to_one",
    )
    outcomes = outcomes.drop(
        columns=["pair_sha256", "origin_cohort", "origin_blind_pair_id"]
    )
    if set(outcomes["blind_pair_id"]) != set(blind["blind_pair_id"]):
        raise RuntimeError("Outcome and blind pair membership differ")
    observed_values = pd.to_numeric(outcomes["observed_pkd"], errors="coerce")
    if not np.all(np.isfinite(observed_values)):
        raise RuntimeError("Nonfinite observed pKd after authorized unblind")
    write_json(
        root / "locks" / "V1_12_OUTCOME_READ_RECEIPT.json",
        {
            "schema_version": "science_advances_v1_12_outcome_read_receipt_v1",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "status": "OUTCOMES_PARSED_AFTER_AUTHORIZED_UNBLIND",
            "rows": len(outcomes),
            "source_rows": source_rows,
            "columns": list(outcomes.columns),
            "outcome_sha256": observed_hashes,
            "prediction_rerun_allowed": False,
        },
    )
    return blind, predictions, outcomes, {
        "paths": {
            **{key: str(value) for key, value in paths.items()},
            "quarantined_outcomes": {
                cohort: str(path) for cohort, path in quarantine_paths.items()
            },
        },
        "checks": checks,
        "protocol": protocol,
        "manifest": manifest,
        "unblind_receipt": unblind_receipt,
    }


def self_test() -> int:
    observed = np.asarray([1, 2, 3, 4, 5, 6], dtype=float)
    predicted = np.asarray([1, 3, 2, 4, 6, 5], dtype=float)
    weights = np.asarray([1, 2, 1, 1, 2, 1], dtype=float)
    expanded_observed = np.repeat(observed, weights.astype(int))
    expanded_predicted = np.repeat(predicted, weights.astype(int))
    expected = float(pd.Series(expanded_observed).corr(pd.Series(expanded_predicted), method="spearman"))
    actual = weighted_spearman(observed, predicted, weights)
    if not np.isclose(actual, expected, atol=1e-12):
        raise RuntimeError(f"Weighted Spearman self-test failed: {actual} != {expected}")
    group_codes = np.asarray([0, 0, 0, 1, 1, 1], dtype=int)
    centered = centered_spearman(
        observed, predicted, group_codes, np.ones(6, dtype=bool), weights
    )
    if not np.isfinite(centered):
        raise RuntimeError("Centered Spearman self-test failed")
    local_components = np.asarray([0, 0, 1, 2, 2, 3], dtype=int)
    dummy = ModelData(
        model_id="self_test",
        pair_ids=np.asarray([f"p{i}" for i in range(6)]),
        observed=observed,
        predicted=predicted,
        interval_low=np.full(6, np.nan),
        interval_high=np.full(6, np.nan),
        protein_labels=np.asarray(["t"] * 6),
        target_codes=np.zeros(6, dtype=int),
        scaffold_labels=np.asarray(["a", "b", "c", "d", "e", "f"]),
        component_codes=local_components,
        source_labels=np.asarray(["s"] * 6),
        eligible_target_codes=np.asarray([0]),
        eligible_rows=np.ones(6, dtype=bool),
        concordance_specs=[],
    )
    dummy.concordance_specs = build_concordance_specs(dummy)
    component_weights = np.asarray([1, 2, 1, 2], dtype=float)
    concordance = concordance_from_component_multiplicities(
        dummy.concordance_specs, component_weights
    )[0]
    pair_numerator = 0.0
    pair_denominator = 0.0
    row_weights = component_weights[local_components]
    for i in range(6):
        for j in range(i + 1, 6):
            score = 0.5
            sign_product = np.sign(observed[i] - observed[j]) * np.sign(
                predicted[i] - predicted[j]
            )
            if sign_product > 0:
                score = 1.0
            elif sign_product < 0:
                score = 0.0
            pair_weight = row_weights[i] * row_weights[j]
            pair_numerator += pair_weight * score
            pair_denominator += pair_weight
    expected_concordance = pair_numerator / pair_denominator
    if not np.isclose(concordance, expected_concordance, atol=1e-12):
        raise RuntimeError("Pairwise concordance self-test failed")
    print(
        json.dumps(
            {
                "status": "SELF_TEST_PASSED",
                "weighted_spearman": actual,
                "centered_spearman": centered,
                "pairwise_concordance": concordance,
                "outcome_rows_read": 0,
            },
            indent=2,
        )
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    root_default = Path(__file__).resolve().parent.parent
    parser.add_argument("--package-root", type=Path, default=root_default)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    root = args.package_root.resolve()
    blind, predictions, outcomes, audit = validate_and_unblind(root)
    results_dir = root / "results_locked"
    results_dir.mkdir(parents=True, exist_ok=True)

    base = blind.merge(
        outcomes, on="blind_pair_id", how="inner", validate="one_to_one",
        suffixes=("", "_outcome"),
    )
    if "measurement_count_outcome" in base.columns:
        blind_counts = pd.to_numeric(base["measurement_count"], errors="coerce")
        outcome_counts = pd.to_numeric(
            base["measurement_count_outcome"], errors="coerce"
        )
        if not blind_counts.equals(outcome_counts):
            raise RuntimeError("Blind and outcome measurement counts differ")
        base["measurement_count"] = outcome_counts
        base = base.drop(columns="measurement_count_outcome")
    long_frame = predictions.merge(
        base, on="blind_pair_id", how="inner", validate="many_to_one"
    )
    components = sorted(base["document_component_id"].astype(str).unique())
    component_map = {value: index for index, value in enumerate(components)}
    model_data = {
        model_id: prepare_model_data(long_frame, model_id, component_map)
        for model_id in ALL_MODELS
    }
    points = {
        model_id: point_core_metrics(model_data[model_id], len(components))
        for model_id in PRIMARY_MODELS
    }

    replicates = int(audit["protocol"]["inference"]["bootstrap_replicates"])
    seed = int(audit["protocol"]["inference"]["seed"])
    rng = np.random.default_rng(seed)
    probabilities = np.full(len(components), 1.0 / len(components), dtype=float)
    component_draws = rng.multinomial(
        len(components), probabilities, size=replicates
    ).astype(np.int16)
    bootstraps = {
        model_id: bootstrap_core_metrics(model_data[model_id], component_draws)
        for model_id in PRIMARY_MODELS
    }
    minimum_valid = int(math.ceil(0.95 * replicates))
    primary, endpoint_meta = simultaneous_primary_inference(
        points, bootstraps, minimum_valid
    )
    decision = decision_from_primary(primary)
    decision.update(
        {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "validation_design": audit["manifest"]["validation_design"],
            "calendar_prospective": False,
            "bootstrap_replicates": replicates,
            "bootstrap_seed": seed,
            "endpoint_inference": endpoint_meta,
        }
    )

    secondary, covariance, mae_comparisons = secondary_metrics(
        long_frame, model_data, len(components)
    )
    sensitivities = point_sensitivity_analysis(
        base, predictions, component_map
    )

    primary_path = results_dir / "V1_12_PRIMARY_ENDPOINTS.tsv"
    decision_path = results_dir / "V1_12_LOCKED_DECISION.json"
    secondary_path = results_dir / "V1_12_SECONDARY_METRICS.tsv"
    covariance_path = results_dir / "V1_12_RANK_COVARIANCE_DECOMPOSITION.tsv"
    mae_path = results_dir / "V1_12_MAE_REFERENCE_COMPARISONS.tsv"
    sensitivity_path = results_dir / "V1_12_LOCKED_SENSITIVITY_POINT_ESTIMATES.tsv"
    bootstrap_path = results_dir / "V1_12_PRIMARY_BOOTSTRAP_REPLICATES.npz"
    primary.to_csv(primary_path, sep="\t", index=False)
    write_json(decision_path, decision)
    secondary.to_csv(secondary_path, sep="\t", index=False)
    covariance.to_csv(covariance_path, sep="\t", index=False)
    mae_comparisons.to_csv(mae_path, sep="\t", index=False)
    sensitivities.to_csv(sensitivity_path, sep="\t", index=False)
    np.savez_compressed(
        bootstrap_path,
        model_ids=np.asarray(PRIMARY_MODELS),
        component_ids=np.asarray(components),
        component_multiplicities=component_draws,
        **{
            f"{endpoint}__{model_id}": bootstraps[model_id][endpoint]
            for endpoint in PRIMARY_ENDPOINTS
            for model_id in PRIMARY_MODELS
        },
    )

    branch_text = {
        "composition_dominated_replication": (
            "The frozen outcome-unexposed temporal aggregate reproduces positive source-wide "
            "association with materially attenuated exact-target ligand ordering."
        ),
        "transported_within_target_ordering": (
            "The frozen outcome-unexposed temporal aggregate supports transport of within-target "
            "ligand ordering without the prespecified material attenuation pattern."
        ),
        "no_temporal_transport": (
            "The frozen outcome-unexposed temporal aggregate does not support source-wide temporal transport."
        ),
        "heterogeneous_or_inconclusive": (
            "The frozen outcome-unexposed temporal aggregate is heterogeneous or inconclusive under "
            "the prespecified four-family decision rule."
        ),
    }[decision["decision_branch"]]
    markdown = [
        "# V1.12 locked outcome-unexposed temporal result",
        "",
        f"**Decision: `{decision['decision_branch']}`**",
        "",
        branch_text,
        "",
        "This is a retrospective, outcome-unexposed, structure-adapted temporal aggregate assembled after structure-only feasibility review. It is not fully protocol-prespecified, calendar-prospective, or source-independent validation.",
        "",
        "## Family decisions",
        "",
        "| Frozen model | Source-wide support | Centered support | Material attenuation | Structural support |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in decision["family_results"]:
        markdown.append(
            f"| {row['model_id']} | {row['source_wide_support']} | "
            f"{row['within_target_centered_support']} | {row['material_attenuation']} | "
            f"{row['family_structural_support']} |"
        )
    markdown.extend(
        [
            "",
            "All four nonreference families, all observed source-origin strata, and every locked sensitivity remain in the machine-readable outputs regardless of direction.",
        ]
    )
    (results_dir / "V1_12_LOCKED_DECISION.md").write_text(
        "\n".join(markdown) + "\n", encoding="utf-8"
    )

    output_files = [
        primary_path,
        decision_path,
        secondary_path,
        covariance_path,
        mae_path,
        sensitivity_path,
        bootstrap_path,
        results_dir / "V1_12_LOCKED_DECISION.md",
    ]
    complete = {
        "schema_version": "science_advances_v1_12_analysis_complete_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "LOCKED_ANALYSIS_COMPLETE",
        "decision_branch": decision["decision_branch"],
        "validation_design": audit["manifest"]["validation_design"],
        "calendar_prospective": False,
        "pairs": len(base),
        "document_components": len(components),
        "models": ALL_MODELS,
        "primary_models": PRIMARY_MODELS,
        "bootstrap_replicates": replicates,
        "bootstrap_seed": seed,
        "analysis_code_sha256": sha256_file(Path(__file__).resolve()),
        "unblind_receipt_sha256": sha256_file(
            root / "locks" / "V1_12_UNBLIND_RECEIPT.json"
        ),
        "outputs": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in output_files
        },
        "prediction_rerun_allowed": False,
    }
    write_json(results_dir / "V1_12_ANALYSIS_COMPLETE.json", complete)
    print(
        json.dumps(
            {
                "status": complete["status"],
                "decision_branch": complete["decision_branch"],
                "pairs": complete["pairs"],
                "document_components": complete["document_components"],
                "bootstrap_replicates": complete["bootstrap_replicates"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
