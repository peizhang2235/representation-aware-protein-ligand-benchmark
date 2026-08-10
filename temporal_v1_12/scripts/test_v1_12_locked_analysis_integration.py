#!/usr/bin/env python3
"""Exercise the locked analysis on a disposable synthetic Stage 2 package."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


ALL_MODELS = [
    "development_exact_median_v1_8",
    "global_top5_morgan_similarity_v1_8",
    "hgb_fusion_additive_pre2018_v1_8",
    "ligand_mpnn_frozen_esm2_v1_8",
    "smiles_protein_cross_attention_v1_8",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    package = Path(__file__).resolve().parent.parent
    analysis_script = package / "scripts" / "run_v1_12_locked_analysis.py"
    temporary = Path(tempfile.mkdtemp(prefix="v112_analysis_integration_"))
    try:
        for directory in (
            "protocol",
            "source_frozen",
            "outcomes_quarantine",
            "predictions_frozen",
            "locks",
        ):
            (temporary / directory).mkdir(parents=True, exist_ok=True)
        protocol = json.loads(
            (package / "protocol" / "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json")
            .read_text(encoding="utf-8")
        )
        protocol["inference"]["bootstrap_replicates"] = 50
        protocol_path = temporary / "protocol" / protocol_path_name()
        write_json(protocol_path, protocol)
        implementation_path = (
            temporary / "protocol" / "V1_12_ANALYSIS_IMPLEMENTATION_SPEC.json"
        )
        shutil.copy2(
            package / "protocol" / "V1_12_ANALYSIS_IMPLEMENTATION_SPEC.json",
            implementation_path,
        )

        rng = np.random.default_rng(13)
        blind_rows = []
        outcome_rows = {"v1.9": [], "v1.10": [], "v1.11": []}
        join_rows = []
        prediction_rows = []
        sources = [
            "strict_dual_timestamp_patent",
            "database_entry_patent",
            "database_entry_chembl_import",
        ]
        for index in range(60):
            target = index // 10
            pair_id = f"SYN_{index:04d}"
            cohort = ["v1.9", "v1.10", "v1.11"][index % 3]
            origin_pair_id = f"{cohort}_SYN_{index:04d}"
            observed = 5.0 + 0.45 * target + 0.08 * (index % 10)
            blind_rows.append(
                {
                    "blind_pair_id": pair_id,
                    "pair_sha256": hashlib.sha256(pair_id.encode()).hexdigest(),
                    "canonical_smiles": "C" * (2 + index % 8),
                    "scaffold_sha256": f"scaffold_{index % 5}",
                    "protein_sequence": "ACDEFGHIKLMNPQRSTVWY" * 2,
                    "protein_sha256": f"protein_{target}",
                    "protein_length": 40,
                    "target_name": f"target {target}",
                    "target_organism": "Homo sapiens",
                    "uniprot_primary_id": f"P{target:05d}",
                    "document_component_id": f"DOCCOMP_{index % 12}",
                    "source_origin_stratum": sources[index % len(sources)],
                    "origin_cohort": cohort,
                    "origin_blind_pair_id": origin_pair_id,
                    "selection_sha256": hashlib.sha256(
                        f"selection_{index}".encode()
                    ).hexdigest(),
                    "component_selection_rank": 1 + index // 12,
                    "measurement_count": 1 + index % 2,
                    "publication_date_min": "2025-01-01",
                    "publication_date_max": "2025-01-01",
                    "bindingdb_date_min": "2025-02-01",
                    "bindingdb_date_max": "2025-02-01",
                    "source_dataset": "synthetic",
                    "temporal_design": "synthetic_firewall_test",
                }
            )
            join_rows.append(
                {
                    "blind_pair_id": pair_id,
                    "pair_sha256": hashlib.sha256(pair_id.encode()).hexdigest(),
                    "origin_cohort": cohort,
                    "origin_blind_pair_id": origin_pair_id,
                }
            )
            outcome_rows[cohort].append(
                {
                    "blind_pair_id": origin_pair_id,
                    "observed_pkd": observed,
                    "kd_nm_median": 10 ** (9 - observed),
                    "measurement_count": 1 + index % 2,
                    "is_exact_10000_nm": False,
                    "measurement_pkd_min": observed,
                    "measurement_pkd_max": observed,
                }
            )
            for model_position, model_id in enumerate(ALL_MODELS):
                if model_position == 0:
                    predicted = 6.0
                else:
                    predicted = (
                        5.0
                        + 0.40 * target
                        + (0.09 - 0.01 * model_position) * (index % 10)
                        + rng.normal(0, 0.03)
                    )
                prediction_rows.append(
                    {
                        "blind_pair_id": pair_id,
                        "model_id": model_id,
                        "family_id": f"family_{model_position}",
                        "predicted_pkd": predicted,
                        "interval_low": predicted - 1.0,
                        "interval_high": predicted + 1.0,
                        "abstained": "false",
                        "failure_reason": "",
                    }
                )
        blind = pd.DataFrame(blind_rows)
        join_map = pd.DataFrame(join_rows)
        predictions = pd.DataFrame(prediction_rows)
        blind_path = temporary / "source_frozen" / "V1_12_BLIND_STRUCTURES.tsv"
        join_map_path = temporary / "source_frozen" / "V1_12_OUTCOME_JOIN_MAP.tsv"
        prediction_path = (
            temporary / "predictions_frozen" / "V1_12_FROZEN_PREDICTIONS.tsv"
        )
        blind.to_csv(blind_path, sep="\t", index=False)
        join_map.to_csv(join_map_path, sep="\t", index=False)
        predictions.to_csv(prediction_path, sep="\t", index=False)
        quarantine_sources = []
        for cohort, rows in outcome_rows.items():
            outcome_path = temporary / "outcomes_quarantine" / f"{cohort}_OUTCOMES.tsv"
            pd.DataFrame(rows).to_csv(outcome_path, sep="\t", index=False)
            outcome_hash = sha256_file(outcome_path)
            os.chmod(outcome_path, 0)
            quarantine_sources.append(
                {
                    "cohort": cohort,
                    "path": str(outcome_path),
                    "precommitted_sha256": outcome_hash,
                    "rows": len(rows),
                }
            )
        manifest = {
            "status": "STRUCTURE_ONLY_MEMBERSHIP_FROZEN",
            "validation_design": "retrospective_outcome_blind_temporal_aggregate",
            "calendar_prospective": False,
            "blind_structure_table_sha256": sha256_file(blind_path),
            "outcome_join_map_sha256": sha256_file(join_map_path),
            "quarantine_sources": quarantine_sources,
        }
        manifest_path = temporary / "source_frozen" / "V1_12_SOURCE_MANIFEST.json"
        write_json(manifest_path, manifest)
        prediction_receipt_path = (
            temporary
            / "predictions_frozen"
            / "V1_12_PREDICTIONS_FROZEN_RECEIPT.json"
        )
        write_json(
            prediction_receipt_path,
            {
                "prediction_file_sha256": sha256_file(prediction_path),
                "blind_structure_table_sha256": sha256_file(blind_path),
            },
        )
        write_json(
            temporary / "locks" / "V1_12_STAGE2_LOCK.json",
            {
                "all_checks_passed": True,
                "outcome_access_authorized": True,
                "quarantined_outcomes": quarantine_sources,
            },
        )
        write_json(
            temporary / "locks" / "V1_12_ANALYSIS_CODE_LOCK.json",
            {
                "analysis_code_sha256": sha256_file(analysis_script),
                "stage1_protocol_sha256": sha256_file(protocol_path),
                "implementation_spec_sha256": sha256_file(implementation_path),
            },
        )
        completed = subprocess.run(
            [
                sys.executable,
                str(analysis_script),
                "--package-root",
                str(temporary),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                "Integration test failed:\n"
                + completed.stdout
                + "\n"
                + completed.stderr
            )
        complete_path = temporary / "results_locked" / "V1_12_ANALYSIS_COMPLETE.json"
        if not complete_path.is_file():
            raise RuntimeError("Analysis completion receipt was not created")
        decision = json.loads(
            (temporary / "results_locked" / "V1_12_LOCKED_DECISION.json").read_text(
                encoding="utf-8"
            )
        )
        print(
            json.dumps(
                {
                    "status": "INTEGRATION_TEST_PASSED",
                    "synthetic_pairs": len(blind),
                    "synthetic_predictions": len(predictions),
                    "decision_branch_exercised": decision["decision_branch"],
                    "real_outcome_rows_read": 0,
                },
                indent=2,
            )
        )
        return 0
    finally:
        for outcome_path in (temporary / "outcomes_quarantine").glob("*_OUTCOMES.tsv"):
            os.chmod(outcome_path, 0o600)
        shutil.rmtree(temporary, ignore_errors=True)


def protocol_path_name() -> str:
    return "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json"


if __name__ == "__main__":
    raise SystemExit(main())
