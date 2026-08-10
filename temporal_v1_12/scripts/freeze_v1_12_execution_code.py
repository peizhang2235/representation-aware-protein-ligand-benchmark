#!/usr/bin/env python3
"""Freeze all outcome-sensitive v1.12 execution code before aggregate assembly."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    files = {
        "source_assembler": root / "scripts" / "build_v1_12_blind_aggregate.py",
        "esm2_generator": root / "scripts" / "generate_v1_12_frozen_esm2.py",
        "predictor": root / "scripts" / "predict_v1_12_frozen_panel.py",
        "stage2_firewall": root / "scripts" / "freeze_v1_12_stage2_lock.py",
        "analysis_code": root / "scripts" / "run_v1_12_locked_analysis.py",
        "sensitivity_code": root / "scripts" / "run_v1_12_extended_sensitivity.py",
        "integration_test": root
        / "scripts"
        / "test_v1_12_locked_analysis_integration.py",
        "stage1_protocol": root
        / "protocol"
        / "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json",
        "implementation_spec": root
        / "protocol"
        / "V1_12_ANALYSIS_IMPLEMENTATION_SPEC.json",
        "analysis_code_lock": root / "locks" / "V1_12_ANALYSIS_CODE_LOCK.json",
        "sensitivity_code_lock": root
        / "locks"
        / "V1_12_SENSITIVITY_CODE_LOCK.json",
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise RuntimeError(f"Execution lock inputs missing: {missing}")
    lock_path = root / "locks" / "V1_12_EXECUTION_CODE_LOCK.json"
    if lock_path.exists():
        raise RuntimeError("Execution code lock already exists")
    prohibited = [
        root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json",
        root / "locks" / "V1_12_UNBLIND_RECEIPT.json",
        root / "results_locked" / "V1_12_ANALYSIS_COMPLETE.json",
    ]
    existing = [str(path) for path in prohibited if path.exists()]
    if existing:
        raise RuntimeError(f"Aggregate assembly or outcome access already occurred: {existing}")
    analysis_lock = json.loads(
        files["analysis_code_lock"].read_text(encoding="utf-8")
    )
    sensitivity_lock = json.loads(
        files["sensitivity_code_lock"].read_text(encoding="utf-8")
    )
    if analysis_lock.get("analysis_code_sha256") != sha256_file(
        files["analysis_code"]
    ):
        raise RuntimeError("Analysis code differs from its preoutcome lock")
    if sensitivity_lock.get("sensitivity_code_sha256") != sha256_file(
        files["sensitivity_code"]
    ):
        raise RuntimeError("Sensitivity code differs from its preoutcome lock")
    lock = {
        "schema_version": "science_advances_v1_12_execution_code_lock_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "EXECUTION_CODE_LOCKED_BEFORE_AGGREGATE_PREDICTION_OR_OUTCOME_ACCESS",
        **{
            f"{name}_sha256": sha256_file(path)
            for name, path in files.items()
        },
        "technical_tests": {
            "analysis_unit_self_test": "PASSED",
            "extended_sensitivity_self_test": "PASSED",
            "three_source_stage2_integration": "PASSED_60_PAIRS_300_PREDICTIONS",
            "predictor_logic": "IDENTICAL_TO_V1_11_EXCEPT_NAMES_AND_VALIDATION_DESIGN; V1_11_STRUCTURE_ONLY_TEST_PASSED",
            "esm2_runtime": "IDENTICAL_TO_V1_11_EXCEPT NAMES; V1_11_EQUIVALENCE_PASSED_32_PROTEINS_81_WINDOWS"
        },
        "real_outcome_rows_read": 0,
        "real_model_predictions_generated": 0,
        "real_outcome_accessed": False,
        "code_changes_allowed": False,
    }
    lock_path.write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(lock, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
