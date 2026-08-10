#!/usr/bin/env python3
"""Freeze the v1.12 additive sensitivity module before outcome access."""

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
    sensitivity = root / "scripts" / "run_v1_12_extended_sensitivity.py"
    analysis = root / "scripts" / "run_v1_12_locked_analysis.py"
    protocol = (
        root
        / "protocol"
        / "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json"
    )
    lock_path = root / "locks" / "V1_12_SENSITIVITY_CODE_LOCK.json"
    if lock_path.exists():
        raise RuntimeError("Sensitivity code lock already exists")
    if any(
        path.exists()
        for path in (
            root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json",
            root / "locks" / "V1_12_UNBLIND_RECEIPT.json",
        )
    ):
        raise RuntimeError("Aggregate source or outcome access already occurred")
    lock = {
        "schema_version": "science_advances_v1_12_sensitivity_code_lock_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "SENSITIVITY_CODE_LOCKED_BEFORE_OUTCOME_ACCESS",
        "sensitivity_code_sha256": sha256_file(sensitivity),
        "analysis_code_sha256": sha256_file(analysis),
        "stage1_protocol_sha256": sha256_file(protocol),
        "self_test_status": "PASSED_BEFORE_LOCK",
        "real_outcome_rows_read": 0,
        "changes_primary_decision": False,
        "code_changes_allowed": False,
    }
    lock_path.write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(lock, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
