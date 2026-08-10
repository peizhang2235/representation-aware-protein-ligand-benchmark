#!/usr/bin/env python3
"""Freeze v1.12 representations and predictions before authorized outcome access."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


EXPECTED_MODELS = {
    "development_exact_median_v1_8",
    "global_top5_morgan_similarity_v1_8",
    "hgb_fusion_additive_pre2018_v1_8",
    "ligand_mpnn_frozen_esm2_v1_8",
    "smiles_protein_cross_attention_v1_8",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    package_root = Path(__file__).resolve().parent.parent
    parser.add_argument("--package-root", type=Path, default=package_root)
    args = parser.parse_args()
    root = args.package_root.resolve()
    paths = {
        "protocol": root
        / "protocol"
        / "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json",
        "implementation_spec": root
        / "protocol"
        / "V1_12_ANALYSIS_IMPLEMENTATION_SPEC.json",
        "manifest": root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json",
        "source_assembler_code": root
        / "scripts"
        / "build_v1_12_blind_aggregate.py",
        "blind": root / "source_frozen" / "V1_12_BLIND_STRUCTURES.tsv",
        "join_map": root / "source_frozen" / "V1_12_OUTCOME_JOIN_MAP.tsv",
        "target_support": root / "source_frozen" / "V1_12_TARGET_SUPPORT.tsv",
        "esm2": root / "model_inputs" / "V1_12_FROZEN_ESM2_EMBEDDINGS.npz",
        "esm2_receipt": root / "model_inputs" / "V1_12_FROZEN_ESM2_RECEIPT.json",
        "esm2_equivalence": root
        / "model_inputs"
        / "V1_12_ESM2_RUNTIME_EQUIVALENCE.json",
        "esm2_code": root / "scripts" / "generate_v1_12_frozen_esm2.py",
        "predictions": root
        / "predictions_frozen"
        / "V1_12_FROZEN_PREDICTIONS.tsv",
        "prediction_start": root
        / "predictions_frozen"
        / "V1_12_PREDICTION_START_RECEIPT.json",
        "prediction_receipt": root
        / "predictions_frozen"
        / "V1_12_PREDICTIONS_FROZEN_RECEIPT.json",
        "prediction_code": root / "scripts" / "predict_v1_12_frozen_panel.py",
        "analysis_code": root / "scripts" / "run_v1_12_locked_analysis.py",
        "analysis_code_lock": root / "locks" / "V1_12_ANALYSIS_CODE_LOCK.json",
        "sensitivity_code": root / "scripts" / "run_v1_12_extended_sensitivity.py",
        "sensitivity_code_lock": root
        / "locks"
        / "V1_12_SENSITIVITY_CODE_LOCK.json",
        "execution_code_lock": root / "locks" / "V1_12_EXECUTION_CODE_LOCK.json",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Stage-2 inputs missing: {missing}")
    lock_path = root / "locks" / "V1_12_STAGE2_LOCK.json"
    if lock_path.exists():
        raise RuntimeError("Stage-2 lock already exists")

    protocol = json.loads(paths["protocol"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    esm2_receipt = json.loads(paths["esm2_receipt"].read_text(encoding="utf-8"))
    prediction_start = json.loads(
        paths["prediction_start"].read_text(encoding="utf-8")
    )
    prediction_receipt = json.loads(
        paths["prediction_receipt"].read_text(encoding="utf-8")
    )
    analysis_lock = json.loads(
        paths["analysis_code_lock"].read_text(encoding="utf-8")
    )
    sensitivity_lock = json.loads(
        paths["sensitivity_code_lock"].read_text(encoding="utf-8")
    )
    execution_lock = json.loads(
        paths["execution_code_lock"].read_text(encoding="utf-8")
    )
    checks = {
        "protocol_locked": protocol.get("status")
        == "OUTCOME_BLIND_ANALYSIS_LOCKED_AFTER_STRUCTURE_ONLY_FEASIBILITY",
        "implementation_locked": json.loads(
            paths["implementation_spec"].read_text(encoding="utf-8")
        ).get("status")
        == "LOCKED_BEFORE_OUTCOME_OR_PREDICTION_ACCESS",
        "source_membership_frozen": manifest.get("status")
        == "STRUCTURE_ONLY_MEMBERSHIP_FROZEN",
        "structural_adequacy_passed": manifest.get(
            "structural_adequacy_gate", {}
        ).get("passed")
        is True,
        "source_manifest_blind_hash_matches": manifest.get(
            "blind_structure_table_sha256"
        )
        == sha256_file(paths["blind"]),
        "source_manifest_join_hash_matches": manifest.get(
            "outcome_join_map_sha256"
        )
        == sha256_file(paths["join_map"]),
        "source_manifest_support_hash_matches": manifest.get(
            "target_support_table_sha256"
        )
        == sha256_file(paths["target_support"]),
        "source_assembler_hash_matches": manifest.get("source_assembler_sha256")
        == sha256_file(paths["source_assembler_code"]),
        "esm2_receipt_source_hash_matches": esm2_receipt.get(
            "blind_structure_table_sha256"
        )
        == sha256_file(paths["blind"]),
        "esm2_hash_matches": esm2_receipt.get("embedding_sha256")
        == sha256_file(paths["esm2"]),
        "esm2_equivalence_passed": esm2_receipt.get(
            "runtime_equivalence_passed"
        )
        is True,
        "esm2_equivalence_hash_matches": esm2_receipt.get(
            "runtime_equivalence_sha256"
        )
        == sha256_file(paths["esm2_equivalence"]),
        "esm2_code_hash_matches": esm2_receipt.get("embedding_code_sha256")
        == sha256_file(paths["esm2_code"]),
        "esm2_outcome_blind": esm2_receipt.get("outcome_rows_read") == 0,
        "prediction_started_outcome_blind": prediction_start.get(
            "future_outcome_rows_read"
        )
        == 0,
        "prediction_frozen_outcome_blind": prediction_receipt.get(
            "future_outcome_rows_read"
        )
        == 0,
        "prediction_hash_matches": prediction_receipt.get("prediction_file_sha256")
        == sha256_file(paths["predictions"]),
        "prediction_code_hash_matches": prediction_receipt.get(
            "prediction_code_sha256"
        )
        == sha256_file(paths["prediction_code"]),
        "prediction_source_hash_matches": prediction_receipt.get(
            "blind_structure_table_sha256"
        )
        == sha256_file(paths["blind"]),
        "prediction_esm2_hash_matches": prediction_receipt.get(
            "future_esm2_embeddings_sha256"
        )
        == sha256_file(paths["esm2"]),
        "calendar_prospective_false": prediction_receipt.get(
            "calendar_prospective"
        )
        is False,
        "analysis_code_locked": analysis_lock.get("analysis_code_sha256")
        == sha256_file(paths["analysis_code"]),
        "sensitivity_code_locked": sensitivity_lock.get(
            "sensitivity_code_sha256"
        )
        == sha256_file(paths["sensitivity_code"]),
        "execution_code_locked": execution_lock.get("status")
        == "EXECUTION_CODE_LOCKED_BEFORE_AGGREGATE_PREDICTION_OR_OUTCOME_ACCESS",
        "locked_source_assembler_matches": execution_lock.get(
            "source_assembler_sha256"
        )
        == sha256_file(paths["source_assembler_code"]),
        "locked_esm2_generator_matches": execution_lock.get(
            "esm2_generator_sha256"
        )
        == sha256_file(paths["esm2_code"]),
        "locked_predictor_matches": execution_lock.get("predictor_sha256")
        == sha256_file(paths["prediction_code"]),
        "locked_stage2_code_matches": execution_lock.get(
            "stage2_firewall_sha256"
        )
        == sha256_file(Path(__file__).resolve()),
        "locked_analysis_matches": execution_lock.get("analysis_code_sha256")
        == sha256_file(paths["analysis_code"]),
        "locked_sensitivity_matches": execution_lock.get(
            "sensitivity_code_sha256"
        )
        == sha256_file(paths["sensitivity_code"]),
    }

    quarantines = []
    for item in manifest.get("quarantine_sources", []):
        cohort = str(item["cohort"])
        outcome_path = Path(str(item["path"]))
        exists = outcome_path.is_file()
        checks[f"{cohort}_quarantine_exists"] = exists
        if not exists:
            continue
        mode = stat.S_IMODE(os.stat(outcome_path).st_mode)
        checks[f"{cohort}_quarantine_restricted"] = (mode & 0o377) == 0
        checks[f"{cohort}_precommitted_hash_present"] = bool(
            item.get("precommitted_sha256")
        )
        quarantines.append(
            {
                "cohort": cohort,
                "path": str(outcome_path),
                "precommitted_sha256": item["precommitted_sha256"],
                "rows": item["rows"],
                "file_mode_before_lock": format(mode, "03o"),
                "read_during_stage2": False,
            }
        )

    blind = pd.read_csv(paths["blind"], sep="\t", dtype=str)
    predictions = pd.read_csv(paths["predictions"], sep="\t", dtype=str)
    observed_models = set(predictions["model_id"])
    checks["model_membership_exact"] = observed_models == EXPECTED_MODELS
    checks["prediction_row_count_exact"] = len(predictions) == len(blind) * len(
        EXPECTED_MODELS
    )
    checks["pair_membership_exact"] = set(predictions["blind_pair_id"]) == set(
        blind["blind_pair_id"]
    )
    checks["no_duplicate_pair_model_rows"] = not predictions.duplicated(
        ["blind_pair_id", "model_id"]
    ).any()
    passed = all(checks.values())
    lock = {
        "schema_version": "science_advances_v1_12_stage2_lock_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": (
            "STAGE2_LOCKED_OUTCOME_ACCESS_AUTHORIZED"
            if passed
            else "STAGE2_LOCK_DENIED"
        ),
        "validation_design": manifest.get("validation_design"),
        "calendar_prospective": False,
        "checks": checks,
        "all_checks_passed": passed,
        "blind_pairs": len(blind),
        "prediction_rows": len(predictions),
        "model_ids": sorted(observed_models),
        "hashes": {name: sha256_file(path) for name, path in paths.items()},
        "quarantined_outcomes": quarantines,
        "prediction_rerun_allowed": False,
        "outcome_access_authorized": passed,
    }
    lock_path.write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(lock, indent=2, sort_keys=True))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
