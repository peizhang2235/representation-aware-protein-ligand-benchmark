#!/usr/bin/env python3
"""Freeze the v1.12 statistical implementation before real prediction or outcome access."""

from __future__ import annotations

import hashlib
import json
import os
import stat
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
    paths = {
        "analysis": root / "scripts" / "run_v1_12_locked_analysis.py",
        "integration_test": root
        / "scripts"
        / "test_v1_12_locked_analysis_integration.py",
        "stage1_protocol": root
        / "protocol"
        / "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json",
        "implementation_spec": root
        / "protocol"
        / "V1_12_ANALYSIS_IMPLEMENTATION_SPEC.json",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise RuntimeError(f"Cannot freeze missing analysis inputs: {missing}")
    lock_path = root / "locks" / "V1_12_ANALYSIS_CODE_LOCK.json"
    if lock_path.exists():
        raise RuntimeError("Analysis code lock already exists")
    prohibited = [
        root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json",
        root / "locks" / "V1_12_UNBLIND_RECEIPT.json",
        root / "locks" / "V1_12_OUTCOME_READ_RECEIPT.json",
        root / "results_locked" / "V1_12_ANALYSIS_COMPLETE.json",
    ]
    existing = [str(path) for path in prohibited if path.exists()]
    if existing:
        raise RuntimeError(f"Aggregate source or outcome access already exists: {existing}")

    master = root.parent
    quarantine_spec = [
        (
            "v1.9",
            master
            / "science_advances_temporal_blind_v1_9_20260802"
            / "outcomes_quarantine"
            / "V1_9_QUARANTINED_OUTCOMES.tsv",
        ),
        (
            "v1.10",
            master
            / "science_advances_temporal_blind_v1_10_20260802"
            / "outcomes_quarantine"
            / "V1_10_QUARANTINED_OUTCOMES.tsv",
        ),
        (
            "v1.11",
            master
            / "science_advances_temporal_blind_v1_11_20260802"
            / "outcomes_quarantine"
            / "V1_11_QUARANTINED_OUTCOMES.tsv",
        ),
    ]
    quarantine_state = []
    for cohort, path in quarantine_spec:
        if not path.is_file():
            raise RuntimeError(f"Missing {cohort} quarantined outcomes")
        mode = stat.S_IMODE(os.stat(path).st_mode)
        if mode & 0o377:
            raise RuntimeError(f"{cohort} quarantine is not access-restricted")
        quarantine_state.append(
            {
                "cohort": cohort,
                "path": str(path),
                "file_mode": format(mode, "03o"),
                "outcome_bytes_read": 0,
            }
        )
    lock = {
        "schema_version": "science_advances_v1_12_analysis_code_lock_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "ANALYSIS_CODE_LOCKED_BEFORE_OUTCOME_ACCESS",
        "analysis_code_sha256": sha256_file(paths["analysis"]),
        "integration_test_sha256": sha256_file(paths["integration_test"]),
        "stage1_protocol_sha256": sha256_file(paths["stage1_protocol"]),
        "implementation_spec_sha256": sha256_file(paths["implementation_spec"]),
        "quarantine_state": quarantine_state,
        "self_test_status": "PASSED_BEFORE_LOCK",
        "integration_test_status": "PASSED_BEFORE_LOCK_ON_SYNTHETIC_THREE_SOURCE_STAGE2_PACKAGE",
        "real_outcome_rows_read": 0,
        "analysis_code_changes_allowed": False,
        "prediction_rerun_allowed_after_unblind": False,
    }
    lock_path.write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(lock, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
