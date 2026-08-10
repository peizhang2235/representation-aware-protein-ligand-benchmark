#!/usr/bin/env python3
"""Freeze the outcome-unexposed v1.12 temporal aggregate membership."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import stat
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCHEMA_VERSION = "science_advances_v1_12_structure_only_source_manifest_v1"
SELECTION_SALT = "science_advances_v1_12_document_balance_v1"
MAX_PER_COMPONENT = 25
FORBIDDEN_COLUMNS = {
    "pkd",
    "observed_pkd",
    "experimental_pkd",
    "measured_pkd",
    "kd",
    "kd_nm",
    "kd_nm_median",
    "standard_value",
    "affinity",
    "affinity_value",
    "outcome",
    "label",
    "residual",
    "error",
    "absolute_error",
}
EXPECTED_MANIFESTS = {
    "v1.9": "14adeb6402f4fe772e8d32ee0056ae12d4d39d5f947a1ec9c0cac4838fb3591b",
    "v1.10": "ccc1a098cc839814e78739556dcca317571e891742bb3ee494709603bab6a0f8",
    "v1.11": "2896b4e12c00e50aa522b413763e3a91f91d7882713cdaf08c2527585686c653",
}
EXPECTED_OUTCOMES = {
    "v1.9": "5052d993563c236c7eb78d4e09b19f605feb7212c9c0c0de0682ceeeb0453fbb",
    "v1.10": "09910ff1cd8ea43026be5273104dfd3a935ddbb4d520c5c9e4c146c80910e587",
    "v1.11": "0e0ddfd0190cb64205c069bce4616c78911454e145f21546d8a545919663415b",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selection_hash(pair_sha256: str) -> str:
    return hashlib.sha256(
        f"{SELECTION_SALT}|{pair_sha256}".encode("utf-8")
    ).hexdigest()


def effective_units(values: pd.Series) -> float:
    counts = values.value_counts().to_numpy(dtype=float)
    return float(counts.sum() ** 2 / np.square(counts).sum())


def cohort_paths(project_root: Path) -> list[dict[str, Any]]:
    master = project_root / "MASTER_SYNTHESIS_20260730"
    return [
        {
            "cohort": "v1.9",
            "package": master / "science_advances_temporal_blind_v1_9_20260802",
            "structures": "V1_9_BLIND_STRUCTURES.tsv",
            "manifest": "V1_9_SOURCE_MANIFEST.json",
            "outcomes": "V1_9_QUARANTINED_OUTCOMES.tsv",
        },
        {
            "cohort": "v1.10",
            "package": master / "science_advances_temporal_blind_v1_10_20260802",
            "structures": "V1_10_BLIND_STRUCTURES.tsv",
            "manifest": "V1_10_SOURCE_MANIFEST.json",
            "outcomes": "V1_10_QUARANTINED_OUTCOMES.tsv",
        },
        {
            "cohort": "v1.11",
            "package": master / "science_advances_temporal_blind_v1_11_20260802",
            "structures": "V1_11_BLIND_STRUCTURES.tsv",
            "manifest": "V1_11_SOURCE_MANIFEST.json",
            "outcomes": "V1_11_QUARANTINED_OUTCOMES.tsv",
        },
    ]


def normalize_cohort(frame: pd.DataFrame, cohort: str) -> pd.DataFrame:
    required = {
        "blind_pair_id",
        "pair_sha256",
        "canonical_smiles",
        "scaffold_sha256",
        "protein_sequence",
        "protein_sha256",
        "protein_length",
        "target_name",
        "target_organism",
        "uniprot_primary_id",
        "document_component_id",
        "measurement_count",
        "publication_date_min",
        "publication_date_max",
        "bindingdb_date_min",
        "bindingdb_date_max",
        "source_dataset",
        "temporal_design",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise RuntimeError(f"{cohort} structure columns missing: {sorted(missing)}")
    forbidden = FORBIDDEN_COLUMNS.intersection(
        {column.strip().lower() for column in frame.columns}
    )
    if forbidden:
        raise RuntimeError(f"Outcome fields entered {cohort}: {sorted(forbidden)}")
    if frame["pair_sha256"].duplicated().any():
        raise RuntimeError(f"{cohort} contains duplicate pair hashes")

    output = frame.loc[:, sorted(required)].copy()
    output = output.rename(columns={"blind_pair_id": "origin_blind_pair_id"})
    output["origin_cohort"] = cohort
    if cohort == "v1.9":
        output["source_origin_stratum"] = "strict_dual_timestamp_patent"
    elif cohort == "v1.10":
        output["source_origin_stratum"] = "database_entry_patent"
    else:
        if "analysis_source_stratum" not in frame.columns:
            raise RuntimeError("v1.11 source stratum is missing")
        strata = frame["analysis_source_stratum"].astype(str).map(
            {
                "patent": "database_entry_residual_patent",
                "chembl_import": "database_entry_chembl_import",
                "pubchem_assay": "database_entry_pubchem_assay",
                "bindingdb_curated_article": "database_entry_curated_article",
                "other_documented": "database_entry_other_documented",
            }
        )
        if strata.isna().any():
            unexpected = sorted(frame.loc[strata.isna(), "analysis_source_stratum"].unique())
            raise RuntimeError(f"Unexpected v1.11 source strata: {unexpected}")
        output["source_origin_stratum"] = strata.to_numpy()
    return output


def main() -> int:
    package_root = Path(__file__).resolve().parent.parent
    project_root = package_root.parents[1]
    protocol_path = (
        package_root
        / "protocol"
        / "V1_12_OUTCOME_BLIND_TEMPORAL_AGGREGATE_PROTOCOL.json"
    )
    specification_path = (
        package_root / "protocol" / "V1_12_ANALYSIS_IMPLEMENTATION_SPEC.json"
    )
    output_dir = package_root / "source_frozen"
    structure_path = output_dir / "V1_12_BLIND_STRUCTURES.tsv"
    join_map_path = output_dir / "V1_12_OUTCOME_JOIN_MAP.tsv"
    support_path = output_dir / "V1_12_TARGET_SUPPORT.tsv"
    manifest_path = output_dir / "V1_12_SOURCE_MANIFEST.json"
    audit_path = output_dir / "V1_12_STRUCTURE_ONLY_FEASIBILITY_AUDIT.json"
    outputs = [structure_path, join_map_path, support_path, manifest_path, audit_path]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise RuntimeError(f"v1.12 frozen membership already exists: {existing}")

    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    if protocol.get("status") != "OUTCOME_BLIND_ANALYSIS_LOCKED_AFTER_STRUCTURE_ONLY_FEASIBILITY":
        raise RuntimeError("v1.12 protocol is not locked")
    if not specification_path.is_file():
        raise RuntimeError("v1.12 implementation specification is missing")

    frames: list[pd.DataFrame] = []
    source_audit: list[dict[str, Any]] = []
    for source in cohort_paths(project_root):
        cohort = source["cohort"]
        source_package = source["package"]
        structure = source_package / "source_frozen" / source["structures"]
        manifest_file = source_package / "source_frozen" / source["manifest"]
        outcome_file = source_package / "outcomes_quarantine" / source["outcomes"]
        observed_manifest_hash = sha256_file(manifest_file)
        if observed_manifest_hash != EXPECTED_MANIFESTS[cohort]:
            raise RuntimeError(f"{cohort} source manifest hash changed")
        manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
        if manifest.get("status") != "INSUFFICIENT_EVIDENCE_GATE":
            raise RuntimeError(f"{cohort} was not closed at its original gate")
        if manifest.get("outcome_values_printed") is not False:
            raise RuntimeError(f"{cohort} does not certify unopened outcomes")
        if manifest.get("quarantine", {}).get("sha256") != EXPECTED_OUTCOMES[cohort]:
            raise RuntimeError(f"{cohort} precommitted outcome hash changed")
        if manifest.get("blind_structure_table_sha256") != sha256_file(structure):
            raise RuntimeError(f"{cohort} blind structure hash changed")
        if not outcome_file.is_file():
            raise RuntimeError(f"{cohort} quarantined outcome file is missing")
        outcome_mode = stat.S_IMODE(os.stat(outcome_file).st_mode)
        if outcome_mode & 0o377:
            raise RuntimeError(
                f"{cohort} outcome quarantine has write, execute, group, or other access"
            )
        frame = pd.read_csv(structure, sep="\t", dtype=str, keep_default_na=False)
        frames.append(normalize_cohort(frame, cohort))
        source_audit.append(
            {
                "cohort": cohort,
                "source_manifest": str(manifest_file),
                "source_manifest_sha256": observed_manifest_hash,
                "blind_structure_table": str(structure),
                "blind_structure_table_sha256": sha256_file(structure),
                "structure_rows": len(frame),
                "quarantined_outcomes": str(outcome_file),
                "precommitted_outcome_sha256": EXPECTED_OUTCOMES[cohort],
                "precommitted_outcome_rows": manifest["quarantine"]["rows"],
                "outcome_file_mode": format(outcome_mode, "03o"),
                "outcome_bytes_read": 0,
            }
        )

    union = pd.concat(frames, ignore_index=True)
    union_before = len(union)
    union = union.drop_duplicates("pair_sha256", keep="first").copy()
    union["selection_sha256"] = union["pair_sha256"].map(selection_hash)
    union = union.sort_values(
        ["document_component_id", "selection_sha256", "pair_sha256"],
        kind="stable",
    )
    union["component_rank"] = union.groupby(
        "document_component_id", sort=False
    ).cumcount()
    selected = union.loc[union["component_rank"] < MAX_PER_COMPONENT].copy()
    selected["component_selection_rank"] = selected["component_rank"] + 1
    selected = selected.sort_values(
        ["selection_sha256", "pair_sha256"], kind="stable"
    ).reset_index(drop=True)
    selected.insert(
        0,
        "blind_pair_id",
        [f"BDBAGG_V112_{index:06d}" for index in range(1, len(selected) + 1)],
    )
    selected = selected.drop(columns="component_rank")

    support = (
        selected.groupby("protein_sha256", sort=True)
        .agg(
            target_name=("target_name", "first"),
            target_organism=("target_organism", "first"),
            uniprot_primary_id=("uniprot_primary_id", "first"),
            pairs=("blind_pair_id", "size"),
            unique_ligands=("canonical_smiles", "nunique"),
            unique_scaffolds=("scaffold_sha256", "nunique"),
            document_components=("document_component_id", "nunique"),
            source_origin_strata=("source_origin_stratum", "nunique"),
        )
        .reset_index()
    )
    support["eligible_within_target"] = (
        (support["pairs"] >= 5)
        & (support["unique_ligands"] >= 5)
        & (support["unique_scaffolds"] >= 3)
    )
    observed = {
        "pairs": int(len(selected)),
        "document_components": int(selected["document_component_id"].nunique()),
        "effective_document_components": effective_units(
            selected["document_component_id"]
        ),
        "exact_proteins": int(selected["protein_sha256"].nunique()),
        "eligible_within_target_groups": int(
            support["eligible_within_target"].sum()
        ),
    }
    required = {
        "pairs": 2500,
        "document_components": 250,
        "effective_document_components": 100.0,
        "exact_proteins": 300,
        "eligible_within_target_groups": 100,
    }
    checks = {key: observed[key] >= value for key, value in required.items()}
    if not all(checks.values()):
        raise RuntimeError(f"v1.12 structural adequacy gate failed: {checks}")

    output_dir.mkdir(parents=True, exist_ok=True)
    blind_columns = [
        "blind_pair_id",
        "pair_sha256",
        "canonical_smiles",
        "scaffold_sha256",
        "protein_sequence",
        "protein_sha256",
        "protein_length",
        "target_name",
        "target_organism",
        "uniprot_primary_id",
        "document_component_id",
        "measurement_count",
        "publication_date_min",
        "publication_date_max",
        "bindingdb_date_min",
        "bindingdb_date_max",
        "source_dataset",
        "temporal_design",
        "origin_cohort",
        "origin_blind_pair_id",
        "source_origin_stratum",
        "selection_sha256",
        "component_selection_rank",
    ]
    selected.loc[:, blind_columns].to_csv(structure_path, sep="\t", index=False)
    selected.loc[
        :,
        [
            "blind_pair_id",
            "pair_sha256",
            "origin_cohort",
            "origin_blind_pair_id",
        ],
    ].to_csv(join_map_path, sep="\t", index=False)
    support.to_csv(support_path, sep="\t", index=False)

    audit = {
        "schema_version": "science_advances_v1_12_structure_only_feasibility_audit_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "PASSED",
        "adaptation_disclosed": True,
        "outcome_values_accessed": False,
        "outcome_bytes_read": 0,
        "source_rows_before_deduplication": union_before,
        "unique_pairs_before_balancing": int(len(union)),
        "duplicate_pair_memberships_removed": int(union_before - len(union)),
        "document_component_cap": MAX_PER_COMPONENT,
        "observed": observed,
        "required": required,
        "checks": checks,
        "origin_rows": {
            str(key): int(value)
            for key, value in selected["origin_cohort"].value_counts().items()
        },
        "source_origin_rows": {
            str(key): int(value)
            for key, value in selected["source_origin_stratum"].value_counts().items()
        },
    }
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "STRUCTURE_ONLY_MEMBERSHIP_FROZEN",
        "validation_design": "retrospective_outcome_blind_temporal_aggregate",
        "calendar_prospective": False,
        "source_independent": False,
        "adaptation_disclosed": True,
        "protocol": str(protocol_path),
        "protocol_sha256": sha256_file(protocol_path),
        "implementation_specification": str(specification_path),
        "implementation_specification_sha256": sha256_file(specification_path),
        "source_inputs": source_audit,
        "membership": {
            "deduplication_priority": ["v1.9", "v1.10", "v1.11"],
            "selection_salt": SELECTION_SALT,
            "maximum_pairs_per_document_component": MAX_PER_COMPONENT,
            "selection_uses_outcomes": False,
        },
        "blind_structure_rows": len(selected),
        "blind_structure_table": str(structure_path),
        "blind_structure_table_sha256": sha256_file(structure_path),
        "outcome_join_map": str(join_map_path),
        "outcome_join_map_sha256": sha256_file(join_map_path),
        "target_support_table": str(support_path),
        "target_support_table_sha256": sha256_file(support_path),
        "structure_only_audit": str(audit_path),
        "structure_only_audit_sha256": sha256_file(audit_path),
        "structural_adequacy_gate": {
            "observed": observed,
            "required": required,
            "checks": checks,
            "passed": True,
        },
        "quarantine_sources": [
            {
                "cohort": item["cohort"],
                "path": item["quarantined_outcomes"],
                "precommitted_sha256": item["precommitted_outcome_sha256"],
                "rows": item["precommitted_outcome_rows"],
                "file_mode_at_membership_freeze": item["outcome_file_mode"],
            }
            for item in source_audit
        ],
        "outcome_fields_present": False,
        "outcome_values_printed": False,
        "outcome_summary_computed": False,
        "outcome_rows_read": 0,
        "source_switching_allowed": False,
        "source_membership_rerun_allowed": False,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "source_assembler_sha256": sha256_file(Path(__file__).resolve()),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "blind_pairs": len(selected),
                "structural_adequacy_gate": manifest["structural_adequacy_gate"],
                "outcome_values_printed": False,
                "outcome_rows_read": 0,
                "blind_structure_table_sha256": manifest[
                    "blind_structure_table_sha256"
                ],
                "source_manifest_sha256": sha256_file(manifest_path),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
