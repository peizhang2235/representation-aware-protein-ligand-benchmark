#!/usr/bin/env python3
"""Audit the editable v1.12 manuscript integration against locked outputs."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


MODELS = {
    "global_top5_morgan_similarity_v1_8",
    "hgb_fusion_additive_pre2018_v1_8",
    "ligand_mpnn_frozen_esm2_v1_8",
    "smiles_protein_cross_attention_v1_8",
}
ENDPOINTS = {
    "source_wide_spearman",
    "within_target_centered_spearman",
    "within_target_pairwise_concordance",
    "target_rank_attenuation",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str, checks: list[dict]) -> None:
    checks.append({"check": message, "passed": bool(condition)})
    if not condition:
        raise AssertionError(message)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    checks: list[dict] = []

    manifest_path = root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json"
    decision_path = root / "results_locked" / "V1_12_LOCKED_DECISION.json"
    primary_path = root / "results_locked" / "V1_12_PRIMARY_ENDPOINTS.tsv"
    manuscript_path = root / "manuscript" / "V1_12_MANUSCRIPT_READY_TEXT.md"
    manifest = read_json(manifest_path)
    decision = read_json(decision_path)
    primary = read_tsv(primary_path)
    manuscript = manuscript_path.read_text(encoding="utf-8")
    manuscript_flat = " ".join(manuscript.split())
    title = " ".join(manuscript.split("## Title", 1)[1].split("## Short title", 1)[0].split())
    short_title = " ".join(
        manuscript.split("## Short title", 1)[1].split("## One-sentence summary", 1)[0].split()
    )
    one_sentence_summary = " ".join(
        manuscript.split("## One-sentence summary", 1)[1].split("## Abstract", 1)[0].split()
    )
    abstract = manuscript.split("## Abstract", 1)[1].split("## Introduction", 1)[0]
    abstract_words = len(abstract.split())

    observed = manifest["structural_adequacy_gate"]["observed"]
    require(manifest["structural_adequacy_gate"]["passed"], "all structural gates passed", checks)
    require(observed["pairs"] == 2689, "pair count is 2689", checks)
    require(observed["exact_proteins"] == 340, "exact-protein count is 340", checks)
    require(observed["document_components"] == 292, "document-component count is 292", checks)
    require(observed["eligible_within_target_groups"] == 126, "eligible-target count is 126", checks)
    require(decision["decision_branch"] == "heterogeneous_or_inconclusive", "locked decision branch retained", checks)
    require(decision["counts"]["source_wide_support"] == 4, "four source-wide supports retained", checks)
    require(decision["counts"]["material_attenuation"] == 1, "one material-attenuation support retained", checks)
    require(len(primary) == 16, "primary table has 16 model-endpoint rows", checks)
    require({row["model_id"] for row in primary} == MODELS, "all four primary models retained", checks)
    require({row["endpoint"] for row in primary} == ENDPOINTS, "all four primary endpoints retained", checks)

    primary_index = {(row["model_id"], row["endpoint"]): row for row in primary}
    source_supported = sum(
        float(primary_index[(model, "source_wide_spearman")]["simultaneous_95_ci_low"]) > 0
        for model in MODELS
    )
    material = sum(
        float(primary_index[(model, "target_rank_attenuation")]["estimate"]) >= 0.10
        and float(primary_index[(model, "target_rank_attenuation")]["simultaneous_95_ci_low"]) > 0
        for model in MODELS
    )
    lower_in_all = all(
        float(primary_index[(model, "within_target_centered_spearman")]["estimate"])
        < float(primary_index[(model, "source_wide_spearman")]["estimate"])
        for model in MODELS
    )
    require(source_supported == 4, "primary table reproduces four source-wide supports", checks)
    require(material == 1, "primary table reproduces one material attenuation", checks)
    require(lower_in_all, "within-target point estimate is lower in all models", checks)

    required_manuscript_tokens = [
        "Outcome-unexposed temporal aggregate separates global transport from model-specific target ranking",
        "heterogeneous_or_inconclusive",
        "Holm-adjusted value across endpoints was P = 0.115",
        "not a calendar-prospective or source-independent replication",
        "No temporal label was used for fitting",
        "owner-readable mode 0400",
        "does not consume or replace that registry",
        "## One-sentence summary",
        "Fig. 6A and table S8",
        "figs. S12 and S13",
    ]
    for token in required_manuscript_tokens:
        require(token in manuscript_flat, f"manuscript contains boundary: {token}", checks)
    require(len(title) <= 125, "title contains no more than 125 characters", checks)
    require(len(short_title) <= 40, "short title contains no more than 40 characters", checks)
    require(len(one_sentence_summary) <= 135, "one-sentence summary contains no more than 135 characters", checks)
    require(abstract_words <= 150, "abstract contains no more than 150 words", checks)
    require("Fig. 7" not in manuscript, "temporal result is not assigned a seventh main figure", checks)

    aliases = {
        root / "tables" / "Table_S9_TEMPORAL_PRIMARY_ENDPOINTS.tsv": primary_path,
        root / "tables" / "Table_S10A_TEMPORAL_LOCKED_SENSITIVITY_POINT_ESTIMATES.tsv": root
        / "results_locked"
        / "V1_12_LOCKED_SENSITIVITY_POINT_ESTIMATES.tsv",
        root / "tables" / "Table_S10B_TEMPORAL_SOURCE_ORIGIN_INTERVALS.tsv": root
        / "results_locked"
        / "extended_sensitivity"
        / "V1_12_SOURCE_ORIGIN_BOOTSTRAP_INTERVALS.tsv",
        root / "tables" / "Table_S10C_TEMPORAL_RESAMPLING_SENSITIVITY.tsv": root
        / "results_locked"
        / "extended_sensitivity"
        / "V1_12_RESAMPLING_UNIT_SENSITIVITY.tsv",
        root / "tables" / "Table_S11_TEMPORAL_SECONDARY_METRICS.tsv": root
        / "results_locked"
        / "V1_12_SECONDARY_METRICS.tsv",
        root / "tables" / "Table_S12_TEMPORAL_RANK_COVARIANCE_DECOMPOSITION.tsv": root
        / "results_locked"
        / "V1_12_RANK_COVARIANCE_DECOMPOSITION.tsv",
    }
    for alias, source in aliases.items():
        require(alias.is_file(), f"submission table exists: {alias.name}", checks)
        require(sha256(alias) == sha256(source), f"submission table matches locked source: {alias.name}", checks)

    figure_stems = [
        root / "figures" / "main" / "Figure_T1_v1_12_temporal_architecture_test",
        root / "figures" / "main" / "Figure_6_outcome_blind_temporal_validation",
        root / "figures" / "supplementary" / "Figure_TS1_v1_12_temporal_architecture_diagnostics",
        root / "figures" / "supplementary" / "Figure_S12_temporal_architecture_diagnostics",
        root / "figures" / "supplementary" / "Figure_TS2_v1_12_resampling_and_source_heterogeneity",
        root / "figures" / "supplementary" / "Figure_S13_resampling_and_source_heterogeneity",
    ]
    figure_files: list[Path] = []
    for stem in figure_stems:
        for suffix in (".pdf", ".svg", ".png", ".tiff"):
            path = stem.with_suffix(suffix)
            require(path.is_file() and path.stat().st_size > 0, f"figure exists and is nonempty: {path.name}", checks)
            figure_files.append(path)

    figure_receipt = read_json(root / "qc" / "V1_12_FIGURE_BUILD_RECEIPT.json")
    require(figure_receipt["status"] == "FIGURES_BUILT_FROM_COMPLETE_LOCKED_RESULTS", "figure receipt is complete", checks)
    require(len(figure_receipt["outputs"]) == 24, "figure receipt contains 24 build and submission outputs", checks)

    deliverables = [
        manuscript_path,
        root / "manuscript" / "V1_12_TEMPORAL_FIGURE_LEGENDS.md",
        root / "manuscript" / "V1_12_TEMPORAL_TABLE_TITLES_AND_NOTES.md",
        root / "submission" / "V1_12_SCIENTIFIC_AND_EDITORIAL_AUDIT.md",
        root / "submission" / "V1_12_COVER_LETTER_DRAFT.md",
        root / "submission" / "V1_12_DATA_AND_CODE_AVAILABILITY_DRAFT.md",
        root / "tables" / "Table_S8_TEMPORAL_DESIGN_AND_FREEZE_AUDIT.tsv",
        *aliases.keys(),
        *figure_files,
    ]
    for path in deliverables:
        require(path.is_file() and path.stat().st_size > 0, f"deliverable exists: {path.name}", checks)

    report = {
        "schema_version": "science_advances_v1_12_submission_integration_audit_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "locked_decision_branch": decision["decision_branch"],
        "scientific_summary": {
            "source_wide_support_models": source_supported,
            "material_attenuation_models": material,
            "within_target_point_lower_in_all_models": lower_in_all,
        },
        "checks": checks,
        "deliverables": [
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in deliverables
        ],
    }
    output_json = root / "qc" / "V1_12_SUBMISSION_INTEGRATION_AUDIT.json"
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    output_md = root / "qc" / "V1_12_SUBMISSION_INTEGRATION_AUDIT.md"
    output_md.write_text(
        "# v1.12 submission integration audit\n\n"
        f"Status: **{report['status']}**\n\n"
        f"Checks passed: {sum(item['passed'] for item in checks)} / {len(checks)}\n\n"
        "Locked scientific result: four of four models support positive source-wide "
        "association; one of four supports material attenuation. The retained branch "
        "is `heterogeneous_or_inconclusive`.\n",
        encoding="utf-8",
    )
    print(f"PASS: {len(checks)} checks; {len(deliverables)} deliverables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
