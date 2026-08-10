#!/usr/bin/env python3
"""Run prelocked source and resampling-unit sensitivities after v1.12 unblinding."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PRIMARY_MODELS = [
    "global_top5_morgan_similarity_v1_8",
    "hgb_fusion_additive_pre2018_v1_8",
    "ligand_mpnn_frozen_esm2_v1_8",
    "smiles_protein_cross_attention_v1_8",
]
ENDPOINTS = [
    "source_wide_spearman",
    "within_target_centered_spearman",
    "within_target_pairwise_concordance",
    "target_rank_attenuation",
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
        raise RuntimeError(f"Expected JSON object: {path}")
    return value


def load_analysis_module(path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(
        "v112_locked_analysis_for_sensitivity", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("Cannot load the locked v1.12 analysis module")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def percentile_summary(
    values: np.ndarray, estimate: float, minimum_valid: int
) -> dict[str, Any]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) < minimum_valid:
        return {
            "estimate": estimate,
            "bootstrap_se": math.nan,
            "ordinary_95_ci_low": math.nan,
            "ordinary_95_ci_high": math.nan,
            "valid_bootstrap_replicates": len(finite),
        }
    return {
        "estimate": estimate,
        "bootstrap_se": float(np.std(finite, ddof=1)),
        "ordinary_95_ci_low": float(
            np.quantile(finite, 0.025, method="lower")
        ),
        "ordinary_95_ci_high": float(
            np.quantile(finite, 0.975, method="higher")
        ),
        "valid_bootstrap_replicates": len(finite),
    }


def prepare_joined(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    blind = pd.read_csv(
        root / "source_frozen" / "V1_12_BLIND_STRUCTURES.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    join_map = pd.read_csv(
        root / "source_frozen" / "V1_12_OUTCOME_JOIN_MAP.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    manifest = load_json(root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json")
    unblind = load_json(root / "locks" / "V1_12_UNBLIND_RECEIPT.json")
    authorized_hashes = {
        str(key): str(value)
        for key, value in unblind.get("outcome_sha256", {}).items()
    }
    source_frames: list[pd.DataFrame] = []
    for item in manifest.get("quarantine_sources", []):
        cohort = str(item["cohort"])
        outcome_path = Path(str(item["path"]))
        observed_hash = sha256_file(outcome_path)
        if authorized_hashes.get(cohort) != observed_hash:
            raise RuntimeError(
                f"{cohort} outcome differs from the authorized unblind receipt"
            )
        local = pd.read_csv(
            outcome_path, sep="\t", dtype=str, keep_default_na=False
        ).rename(columns={"blind_pair_id": "origin_blind_pair_id"})
        local.insert(0, "origin_cohort", cohort)
        source_frames.append(local)
    all_outcomes = pd.concat(source_frames, ignore_index=True)
    outcomes = join_map.merge(
        all_outcomes,
        on=["origin_cohort", "origin_blind_pair_id"],
        how="left",
        validate="one_to_one",
    ).drop(columns=["pair_sha256", "origin_cohort", "origin_blind_pair_id"])
    predictions = pd.read_csv(
        root / "predictions_frozen" / "V1_12_FROZEN_PREDICTIONS.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    base = blind.merge(
        outcomes,
        on="blind_pair_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_outcome"),
    )
    if "measurement_count_outcome" in base.columns:
        base = base.drop(columns="measurement_count").rename(
            columns={"measurement_count_outcome": "measurement_count"}
        )
    long_frame = predictions.merge(
        base, on="blind_pair_id", how="inner", validate="many_to_one"
    )
    return base, long_frame


def source_bootstraps(
    analysis: Any,
    base: pd.DataFrame,
    long_frame: pd.DataFrame,
    component_map: dict[str, int],
    component_draws: np.ndarray,
    minimum_valid: int,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    masks = [
        (name, category, mask)
        for name, category, mask in analysis.sensitivity_masks(base)
        if category in {
            "mandatory_source_origin",
            "mandatory_leave_one_origin_out",
        }
    ]
    records: list[dict[str, Any]] = []
    distributions: dict[str, np.ndarray] = {}
    for name, category, mask in masks:
        pair_ids = set(base.loc[mask, "blind_pair_id"].astype(str))
        subset = long_frame.loc[
            long_frame["blind_pair_id"].astype(str).isin(pair_ids)
        ].copy()
        for model_id in PRIMARY_MODELS:
            data = analysis.prepare_model_data(subset, model_id, component_map)
            points = analysis.point_core_metrics(data, len(component_map))
            boot = analysis.bootstrap_core_metrics(data, component_draws)
            for endpoint in ENDPOINTS:
                summary = percentile_summary(
                    boot[endpoint], points[endpoint], minimum_valid
                )
                records.append(
                    {
                        "sensitivity": name,
                        "category": category,
                        "model_id": model_id,
                        "endpoint": endpoint,
                        "nonabstained_pairs": len(data.observed),
                        "exact_proteins": len(np.unique(data.protein_labels)),
                        "eligible_within_target_groups": len(
                            data.eligible_target_codes
                        ),
                        **summary,
                    }
                )
                distributions[f"{name}__{model_id}__{endpoint}"] = boot[endpoint]
    return pd.DataFrame(records), distributions


def unit_concordance_specs(
    analysis: Any, data: Any, row_unit_codes: np.ndarray
) -> list[Any]:
    specs = []
    for target_code in data.eligible_target_codes:
        indices = np.flatnonzero(data.target_codes == target_code)
        units, local_codes = np.unique(
            row_unit_codes[indices], return_inverse=True
        )
        unit_rows = np.bincount(local_codes, minlength=len(units)).astype(float)
        signed = analysis.aggregate_signed_pairs(
            data.observed[indices],
            data.predicted[indices],
            local_codes,
            len(units),
        )
        specs.append(
            analysis.ConcordanceSpec(
                target_code=int(target_code),
                component_codes=units.astype(int),
                component_row_counts=unit_rows,
                signed_pair_matrix=signed,
            )
        )
    return specs


def target_cluster_concordance(
    analysis: Any, data: Any, target_draws: np.ndarray
) -> np.ndarray:
    values = []
    codes = []
    for target_code in data.eligible_target_codes:
        indices = np.flatnonzero(data.target_codes == target_code)
        observed = data.observed[indices]
        predicted = data.predicted[indices]
        numerator = 0.0
        denominator = 0.0
        for left in range(len(indices)):
            for right in range(left + 1, len(indices)):
                product = np.sign(observed[left] - observed[right]) * np.sign(
                    predicted[left] - predicted[right]
                )
                numerator += 1.0 if product > 0 else 0.0 if product < 0 else 0.5
                denominator += 1.0
        if denominator > 0:
            values.append(numerator / denominator)
            codes.append(int(target_code))
    values_array = np.asarray(values, dtype=float)
    codes_array = np.asarray(codes, dtype=int)
    local = target_draws[:, codes_array].astype(float)
    denominator = local.sum(axis=1)
    return np.divide(
        local @ values_array,
        denominator,
        out=np.full(len(local), np.nan),
        where=denominator > 0,
    )


def alternative_unit_bootstraps(
    analysis: Any,
    base: pd.DataFrame,
    long_frame: pd.DataFrame,
    component_map: dict[str, int],
    replicates: int,
    seed: int,
    minimum_valid: int,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    unit_columns = [
        ("exact_protein", "protein_sha256"),
        ("scaffold", "scaffold_sha256"),
        ("source_origin", "source_origin_stratum"),
    ]
    base_index = base.set_index("blind_pair_id")
    records: list[dict[str, Any]] = []
    distributions: dict[str, np.ndarray] = {}
    for unit_position, (unit_name, column) in enumerate(unit_columns, start=1):
        for model_position, model_id in enumerate(PRIMARY_MODELS, start=1):
            data = analysis.prepare_model_data(
                long_frame, model_id, component_map
            )
            labels = base_index.loc[data.pair_ids, column].astype(str).to_numpy()
            row_unit_codes, unit_labels = pd.factorize(labels, sort=True)
            unit_count = len(unit_labels)
            rng = np.random.default_rng(seed + 1000 * unit_position + model_position)
            draws = rng.multinomial(
                unit_count,
                np.full(unit_count, 1.0 / unit_count),
                size=replicates,
            ).astype(np.int16)
            source = np.full(replicates, np.nan)
            centered = np.full(replicates, np.nan)
            for replicate in range(replicates):
                weights = draws[replicate, row_unit_codes].astype(float)
                source[replicate] = analysis.weighted_spearman(
                    data.observed, data.predicted, weights
                )
                centered[replicate] = analysis.centered_spearman(
                    data.observed,
                    data.predicted,
                    data.target_codes,
                    data.eligible_rows,
                    weights,
                )
            if unit_name == "exact_protein":
                concordance = target_cluster_concordance(
                    analysis, data, draws
                )
            else:
                specs = unit_concordance_specs(
                    analysis, data, row_unit_codes
                )
                concordance = analysis.concordance_from_component_multiplicities(
                    specs, draws
                )
            boot = {
                "source_wide_spearman": source,
                "within_target_centered_spearman": centered,
                "within_target_pairwise_concordance": concordance,
                "target_rank_attenuation": source - centered,
            }
            points = analysis.point_core_metrics(data, len(component_map))
            fragility = (
                "fewer_than_10_resampling_units"
                if unit_count < 10
                else "none"
            )
            for endpoint in ENDPOINTS:
                summary = percentile_summary(
                    boot[endpoint], points[endpoint], minimum_valid
                )
                records.append(
                    {
                        "resampling_unit": unit_name,
                        "model_id": model_id,
                        "endpoint": endpoint,
                        "resampling_units": unit_count,
                        "fragility_flag": fragility,
                        **summary,
                    }
                )
                distributions[f"{unit_name}__{model_id}__{endpoint}"] = boot[
                    endpoint
                ]
    return pd.DataFrame(records), distributions


def self_test() -> int:
    values = np.asarray([0.1, 0.2, np.nan, 0.3])
    summary = percentile_summary(values, 0.2, 3)
    if summary["valid_bootstrap_replicates"] != 3:
        raise RuntimeError("Extended sensitivity summary self-test failed")
    print(
        json.dumps(
            {
                "status": "SELF_TEST_PASSED",
                "finite_replicates": summary["valid_bootstrap_replicates"],
                "real_outcome_rows_read": 0,
            },
            indent=2,
        )
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    default_root = Path(__file__).resolve().parent.parent
    parser.add_argument("--package-root", type=Path, default=default_root)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    root = args.package_root.resolve()
    analysis_path = root / "scripts" / "run_v1_12_locked_analysis.py"
    sensitivity_lock_path = root / "locks" / "V1_12_SENSITIVITY_CODE_LOCK.json"
    analysis_complete_path = root / "results_locked" / "V1_12_ANALYSIS_COMPLETE.json"
    unblind_path = root / "locks" / "V1_12_UNBLIND_RECEIPT.json"
    required = [
        analysis_path,
        sensitivity_lock_path,
        analysis_complete_path,
        unblind_path,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"Extended sensitivity prerequisites missing: {missing}")
    lock = load_json(sensitivity_lock_path)
    if lock.get("sensitivity_code_sha256") != sha256_file(Path(__file__).resolve()):
        raise RuntimeError("Extended sensitivity code differs from preoutcome lock")
    if lock.get("analysis_code_sha256") != sha256_file(analysis_path):
        raise RuntimeError("Primary analysis code differs from sensitivity lock")
    complete = load_json(analysis_complete_path)
    unblind = load_json(unblind_path)
    manifest = load_json(root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json")
    authorized_hashes = {
        str(key): str(value)
        for key, value in unblind.get("outcome_sha256", {}).items()
    }
    for item in manifest.get("quarantine_sources", []):
        cohort = str(item["cohort"])
        outcome_path = Path(str(item["path"]))
        if authorized_hashes.get(cohort) != sha256_file(outcome_path):
            raise RuntimeError(
                f"{cohort} outcome hash differs from authorized unblind receipt"
            )
    output_dir = root / "results_locked" / "extended_sensitivity"
    if output_dir.exists():
        raise RuntimeError("Extended sensitivity outputs already exist; rerun prohibited")
    output_dir.mkdir(parents=True)

    analysis = load_analysis_module(analysis_path)
    base, long_frame = prepare_joined(root)
    bootstrap_path = root / "results_locked" / "V1_12_PRIMARY_BOOTSTRAP_REPLICATES.npz"
    with np.load(bootstrap_path, allow_pickle=False) as archive:
        component_ids = archive["component_ids"].astype(str)
        component_draws = archive["component_multiplicities"].astype(np.int16)
    component_map = {
        value: index for index, value in enumerate(component_ids)
    }
    observed_components = set(base["document_component_id"].astype(str))
    if observed_components != set(component_ids):
        raise RuntimeError("Saved bootstrap component membership changed")
    replicates = int(component_draws.shape[0])
    minimum_valid = int(math.ceil(0.95 * replicates))
    seed = int(complete["bootstrap_seed"])

    source_table, source_distributions = source_bootstraps(
        analysis,
        base,
        long_frame,
        component_map,
        component_draws,
        minimum_valid,
    )
    unit_table, unit_distributions = alternative_unit_bootstraps(
        analysis,
        base,
        long_frame,
        component_map,
        replicates,
        seed,
        minimum_valid,
    )
    source_path = output_dir / "V1_12_SOURCE_ORIGIN_BOOTSTRAP_INTERVALS.tsv"
    unit_path = output_dir / "V1_12_RESAMPLING_UNIT_SENSITIVITY.tsv"
    distributions_path = output_dir / "V1_12_EXTENDED_SENSITIVITY_REPLICATES.npz"
    source_table.to_csv(source_path, sep="\t", index=False)
    unit_table.to_csv(unit_path, sep="\t", index=False)
    np.savez_compressed(
        distributions_path,
        **source_distributions,
        **unit_distributions,
    )
    receipt = {
        "schema_version": "science_advances_v1_12_extended_sensitivity_receipt_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "PRELOCKED_EXTENDED_SENSITIVITIES_COMPLETE",
        "primary_decision_unchanged": complete["decision_branch"],
        "bootstrap_replicates": replicates,
        "minimum_valid_replicates": minimum_valid,
        "source_origin_rows": len(source_table),
        "resampling_unit_rows": len(unit_table),
        "source_origins_never_filtered_by_direction": True,
        "source_origin_intervals": "ordinary component-bootstrap percentile intervals; descriptive sensitivity, not simultaneous primary inference",
        "source_origin_resampling_warning": "Fewer than 10 source origins yields a fragile, discrete bootstrap and is labeled as such.",
        "sensitivity_code_sha256": sha256_file(Path(__file__).resolve()),
        "analysis_complete_sha256": sha256_file(analysis_complete_path),
        "outputs": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in (source_path, unit_path, distributions_path)
        },
    }
    receipt_path = output_dir / "V1_12_EXTENDED_SENSITIVITY_RECEIPT.json"
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "source_origin_rows": len(source_table),
                "resampling_unit_rows": len(unit_table),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
