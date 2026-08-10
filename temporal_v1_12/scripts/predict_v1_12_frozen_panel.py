#!/usr/bin/env python3
"""Apply the unchanged v1.8 model panel to the v1.12 blind temporal source."""

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
import torch


REQUIRED_INPUT_COLUMNS = {
    "blind_pair_id",
    "canonical_smiles",
    "protein_sequence",
    "protein_sha256",
}
FORBIDDEN_INPUT_COLUMNS = {
    "pkd",
    "observed_pkd",
    "experimental_pkd",
    "measured_pkd",
    "kd",
    "kd_nm",
    "standard_value",
    "affinity",
    "affinity_value",
    "outcome",
    "label",
    "residual",
    "error",
    "absolute_error",
}
MODEL_ORDER = [
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


def finite_q90(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    probability = min(1.0, math.ceil((len(values) + 1) * 0.90) / len(values))
    return float(np.quantile(values, probability, method="higher"))


def load_training_module(path: Path) -> Any:
    specification = importlib.util.spec_from_file_location("v18_training_module", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot import frozen v1.8 model definitions: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def output_row(
    *,
    blind_pair_id: str,
    model_id: str,
    family_id: str,
    prediction: float | None,
    q90: float | None,
    failure_reason: str = "",
) -> dict[str, Any]:
    abstained = prediction is None or not np.isfinite(prediction)
    if abstained:
        return {
            "blind_pair_id": blind_pair_id,
            "model_id": model_id,
            "family_id": family_id,
            "predicted_pkd": "",
            "interval_low": "",
            "interval_high": "",
            "abstained": "true",
            "failure_reason": failure_reason or "prediction_unavailable",
        }
    lower = prediction - q90 if q90 is not None else math.nan
    upper = prediction + q90 if q90 is not None else math.nan
    return {
        "blind_pair_id": blind_pair_id,
        "model_id": model_id,
        "family_id": family_id,
        "predicted_pkd": float(prediction),
        "interval_low": float(lower) if np.isfinite(lower) else "",
        "interval_high": float(upper) if np.isfinite(upper) else "",
        "abstained": "false",
        "failure_reason": "",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    package_root = Path(__file__).resolve().parent.parent
    project_root = package_root.parents[1]
    parser.add_argument("--package-root", type=Path, default=package_root)
    parser.add_argument("--project-root", type=Path, default=project_root)
    parser.add_argument("--blind-structures", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--future-esm2-embeddings", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=package_root / "predictions_frozen",
    )
    args = parser.parse_args()
    package_root = args.package_root.resolve()
    project_root = args.project_root.resolve()
    model_package_root = (
        project_root
        / "MASTER_SYNTHESIS_20260730"
        / "science_advances_prospective_v1_8_20260802"
    )
    output_dir = args.output_dir.resolve()

    prediction_path = output_dir / "V1_12_FROZEN_PREDICTIONS.tsv"
    start_path = output_dir / "V1_12_PREDICTION_START_RECEIPT.json"
    receipt_path = output_dir / "V1_12_PREDICTIONS_FROZEN_RECEIPT.json"
    existing = [path for path in (prediction_path, start_path, receipt_path) if path.exists()]
    if existing:
        raise RuntimeError(f"Prediction outputs already exist: {existing}")

    source_manifest = json.loads(args.source_manifest.read_text(encoding="utf-8"))
    if source_manifest.get("status") != "STRUCTURE_ONLY_MEMBERSHIP_FROZEN":
        raise RuntimeError("Source structure-only membership is not frozen")
    expected_design = (
        "retrospective_outcome_blind_temporal_aggregate"
    )
    if source_manifest.get("validation_design") != expected_design:
        raise RuntimeError("Source does not match the locked v1.12 validation design")
    if source_manifest.get("calendar_prospective") is not False:
        raise RuntimeError("v1.12 must remain explicitly non-calendar-prospective")
    if source_manifest.get("outcome_fields_present") is not False:
        raise RuntimeError("Source manifest does not certify outcome-free structures")
    if source_manifest.get("blind_structure_table_sha256") != sha256_file(
        args.blind_structures
    ):
        raise RuntimeError("Blind structure table hash differs from source manifest")

    frame = pd.read_csv(
        args.blind_structures,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    normalized_columns = {column.strip().lower() for column in frame.columns}
    missing = REQUIRED_INPUT_COLUMNS.difference(frame.columns)
    forbidden = FORBIDDEN_INPUT_COLUMNS.intersection(normalized_columns)
    if missing:
        raise RuntimeError(f"Blind structure columns missing: {sorted(missing)}")
    if forbidden:
        raise RuntimeError(f"Outcome fields entered prediction input: {sorted(forbidden)}")
    if len(frame) == 0 or frame["blind_pair_id"].duplicated().any():
        raise RuntimeError("Blind pair IDs must be nonempty and unique")

    # Import the locked environment before prepending the isolated package path.
    runtime_path = model_package_root / "runtime" / "python_packages"
    sys.path.insert(0, str(runtime_path))
    import sklearn  # noqa: PLC0415
    from rdkit import Chem, DataStructs, rdBase  # noqa: PLC0415
    from rdkit.Chem import rdFingerprintGenerator  # noqa: PLC0415

    if sklearn.__version__ != "1.4.2":
        raise RuntimeError(f"scikit-learn 1.4.2 required, got {sklearn.__version__}")
    if rdBase.rdkitVersion != "2023.09.5":
        raise RuntimeError(f"RDKit 2023.09.5 required, got {rdBase.rdkitVersion}")
    sys.path.insert(0, str(project_root / "RareMol-AI_platform"))
    import joblib  # noqa: PLC0415

    training_module_path = (
        model_package_root / "scripts" / "train_v1_8_missing_model_families.py"
    )
    training_module = load_training_module(training_module_path)
    artifacts = {
        "local": model_package_root
        / "model_artifacts"
        / "local_similarity"
        / "V1_8_LOCAL_SIMILARITY_INDEX.npz",
        "hgb": model_package_root
        / "model_artifacts"
        / "existing"
        / "hgb_fusion_additive.joblib",
        "graph": model_package_root
        / "model_artifacts"
        / "graph_sequence"
        / "V1_8_GRAPH_SEQUENCE_MODEL.pt",
        "joint": model_package_root
        / "model_artifacts"
        / "joint_interaction"
        / "V1_8_JOINT_INTERACTION_MODEL.pt",
    }
    missing_artifacts = [str(path) for path in artifacts.values() if not path.is_file()]
    if missing_artifacts:
        raise RuntimeError(f"Frozen model artifacts missing: {missing_artifacts}")

    output_dir.mkdir(parents=True, exist_ok=True)
    started_at = datetime.now(timezone.utc).isoformat()
    start_receipt = {
        "schema_version": "science_advances_v1_12_prediction_start_receipt_v1",
        "created_at": started_at,
        "status": "BLIND_PREDICTION_STARTED_OUTCOMES_UNOPENED",
        "validation_design": expected_design,
        "model_panel_source": str(model_package_root),
        "source_manifest_sha256": sha256_file(args.source_manifest),
        "blind_structure_table_sha256": sha256_file(args.blind_structures),
        "future_outcome_rows_read": 0,
        "rerun_allowed": False,
    }
    start_path.write_text(json.dumps(start_receipt, indent=2) + "\n", encoding="utf-8")

    with np.load(artifacts["local"], allow_pickle=False) as archive:
        training_bits = archive["fingerprint"].astype(np.uint8)
        training_pkd = archive["pkd"].astype(float)
        local_q90 = float(archive["conformal_q90"][0])
    median_prediction = float(np.median(training_pkd))
    constant_q90 = finite_q90(np.abs(training_pkd - median_prediction))
    training_fingerprints = []
    for bit_row in training_bits:
        fingerprint = DataStructs.ExplicitBitVect(training_bits.shape[1])
        for bit in np.flatnonzero(bit_row):
            fingerprint.SetBit(int(bit))
        training_fingerprints.append(fingerprint)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    output_rows: list[dict[str, Any]] = []
    for row in frame.itertuples(index=False):
        blind_id = str(row.blind_pair_id)
        output_rows.append(
            output_row(
                blind_pair_id=blind_id,
                model_id="development_exact_median_v1_8",
                family_id="constant_reference",
                prediction=median_prediction,
                q90=constant_q90,
            )
        )
        molecule = Chem.MolFromSmiles(str(row.canonical_smiles))
        if molecule is None:
            output_rows.append(
                output_row(
                    blind_pair_id=blind_id,
                    model_id="global_top5_morgan_similarity_v1_8",
                    family_id="local_similarity",
                    prediction=None,
                    q90=None,
                    failure_reason="invalid_smiles",
                )
            )
        else:
            query = generator.GetFingerprint(molecule)
            similarities = np.asarray(
                DataStructs.BulkTanimotoSimilarity(query, training_fingerprints),
                dtype=float,
            )
            order = np.argsort(-similarities, kind="stable")[:5]
            weights = np.maximum(similarities[order], 1e-6) ** 4.0
            prediction = float(np.average(training_pkd[order], weights=weights))
            output_rows.append(
                output_row(
                    blind_pair_id=blind_id,
                    model_id="global_top5_morgan_similarity_v1_8",
                    family_id="local_similarity",
                    prediction=prediction,
                    q90=local_q90,
                )
            )

    hgb_artifact = joblib.load(artifacts["hgb"])
    for row in frame.itertuples(index=False):
        result = hgb_artifact.predict_one(
            str(row.canonical_smiles), str(row.protein_sequence)
        )
        prediction = result.get("score") if result.get("status") != "unavailable" else None
        output_rows.append(
            output_row(
                blind_pair_id=str(row.blind_pair_id),
                model_id="hgb_fusion_additive_pre2018_v1_8",
                family_id="histogram_gradient_boosting",
                prediction=float(prediction) if prediction is not None else None,
                q90=float(hgb_artifact.conformal_q90),
                failure_reason=str(result.get("reason") or ""),
            )
        )

    with np.load(args.future_esm2_embeddings, allow_pickle=False) as archive:
        future_hashes = archive["protein_sha256"].astype(str)
        future_embeddings = archive["embedding"].astype(np.float32)
    if future_embeddings.shape != (len(future_hashes), 320):
        raise RuntimeError("Future ESM2 embedding array has an unexpected shape")
    embedding_lookup = {
        value: future_embeddings[index] for index, value in enumerate(future_hashes)
    }
    graph_valid = np.asarray(
        [
            str(row.protein_sha256) in embedding_lookup
            and Chem.MolFromSmiles(str(row.canonical_smiles)) is not None
            for row in frame.itertuples(index=False)
        ],
        dtype=bool,
    )
    graph_predictions: dict[int, float] = {}
    graph_bundle = torch.load(artifacts["graph"], map_location="cpu")
    graph_config = graph_bundle["configuration"]
    graph_model = training_module.GraphSequenceModel(
        hidden=int(graph_config["hidden"]),
        layers=int(graph_config["message_passing_layers"]),
    )
    graph_model.load_state_dict(graph_bundle["state_dict"])
    if graph_valid.any():
        graph_frame = frame.loc[graph_valid].copy().reset_index(drop=True)
        graph_frame["pkd"] = 0.0
        graph_frame["sample_weight"] = 1.0
        graph_dataset = training_module.GraphDataset(
            graph_frame,
            np.arange(len(graph_frame)),
            embedding_lookup,
            Chem,
        )
        graph_loader = training_module.make_loader(
            graph_dataset,
            batch_size=128,
            shuffle=False,
            seed=int(graph_config["seed"]),
            collate=training_module.collate_graph,
        )
        _, _, values = training_module.predict_loader(graph_model, graph_loader)
        for original_index, value in zip(np.flatnonzero(graph_valid), values):
            graph_predictions[int(original_index)] = float(value)
    for index, row in enumerate(frame.itertuples(index=False)):
        prediction = graph_predictions.get(index)
        reason = ""
        if prediction is None:
            reason = (
                "missing_frozen_esm2_embedding"
                if str(row.protein_sha256) not in embedding_lookup
                else "invalid_smiles"
            )
        output_rows.append(
            output_row(
                blind_pair_id=str(row.blind_pair_id),
                model_id="ligand_mpnn_frozen_esm2_v1_8",
                family_id="graph_sequence_neural",
                prediction=prediction,
                q90=float(graph_config["conformal_q90"]),
                failure_reason=reason,
            )
        )

    joint_bundle = torch.load(artifacts["joint"], map_location="cpu")
    joint_config = joint_bundle["configuration"]
    joint_model = training_module.JointInteractionModel(
        len(joint_bundle["smiles_vocab"]),
        len(joint_bundle["protein_vocab"]),
        dimension=int(joint_config["dimension"]),
        heads=int(joint_config["attention_heads"]),
        max_smiles=int(joint_config["max_smiles_tokens"]),
        max_protein=int(joint_config["max_protein_tokens"]),
    )
    joint_model.load_state_dict(joint_bundle["state_dict"])
    joint_valid = np.asarray(
        [
            bool(str(row.canonical_smiles).strip())
            and bool(str(row.protein_sequence).strip())
            for row in frame.itertuples(index=False)
        ],
        dtype=bool,
    )
    joint_predictions: dict[int, float] = {}
    if joint_valid.any():
        joint_frame = frame.loc[joint_valid].copy().reset_index(drop=True)
        joint_frame["pkd"] = 0.0
        joint_frame["sample_weight"] = 1.0
        joint_dataset = training_module.JointDataset(
            joint_frame,
            np.arange(len(joint_frame)),
            joint_bundle["smiles_vocab"],
            joint_bundle["protein_vocab"],
            int(joint_config["max_smiles_tokens"]),
            int(joint_config["max_protein_tokens"]),
        )
        joint_loader = training_module.make_loader(
            joint_dataset,
            batch_size=64,
            shuffle=False,
            seed=int(joint_config["seed"]),
            collate=training_module.collate_joint,
        )
        _, _, values = training_module.predict_loader(joint_model, joint_loader)
        for original_index, value in zip(np.flatnonzero(joint_valid), values):
            joint_predictions[int(original_index)] = float(value)
    for index, row in enumerate(frame.itertuples(index=False)):
        output_rows.append(
            output_row(
                blind_pair_id=str(row.blind_pair_id),
                model_id="smiles_protein_cross_attention_v1_8",
                family_id="joint_interaction_neural",
                prediction=joint_predictions.get(index),
                q90=float(joint_config["conformal_q90"]),
                failure_reason=(
                    "missing_structure_or_sequence" if index not in joint_predictions else ""
                ),
            )
        )

    predictions = pd.DataFrame(output_rows)
    model_rank = {model_id: index for index, model_id in enumerate(MODEL_ORDER)}
    predictions["model_order"] = predictions["model_id"].map(model_rank)
    predictions = predictions.sort_values(
        ["blind_pair_id", "model_order"], kind="stable"
    ).drop(columns="model_order")
    expected_rows = len(frame) * len(MODEL_ORDER)
    if len(predictions) != expected_rows:
        raise RuntimeError(f"Expected {expected_rows} predictions, got {len(predictions)}")
    observed_models = set(predictions["model_id"])
    if observed_models != set(MODEL_ORDER):
        raise RuntimeError("Frozen prediction model membership changed")
    predictions.to_csv(prediction_path, sep="\t", index=False)
    completed_at = datetime.now(timezone.utc).isoformat()
    receipt = {
        "schema_version": "science_advances_v1_12_predictions_frozen_receipt_v1",
        "created_at": completed_at,
        "status": "PREDICTIONS_FROZEN_OUTCOMES_UNOPENED",
        "validation_design": expected_design,
        "calendar_prospective": False,
        "outcome_blind_temporal_stage2_eligible": True,
        "model_panel_source": str(model_package_root),
        "blind_pairs": int(len(frame)),
        "prediction_rows": int(len(predictions)),
        "model_ids": MODEL_ORDER,
        "abstention_counts": {
            str(model_id): int(
                predictions.loc[predictions["model_id"] == model_id, "abstained"]
                .eq("true")
                .sum()
            )
            for model_id in MODEL_ORDER
        },
        "source_manifest_sha256": sha256_file(args.source_manifest),
        "blind_structure_table_sha256": sha256_file(args.blind_structures),
        "future_esm2_embeddings_sha256": sha256_file(args.future_esm2_embeddings),
        "model_artifact_sha256": {
            key: sha256_file(path) for key, path in artifacts.items()
        },
        "prediction_file": str(prediction_path),
        "prediction_file_sha256": sha256_file(prediction_path),
        "prediction_code_sha256": sha256_file(Path(__file__).resolve()),
        "runtime": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "torch": torch.__version__,
            "rdkit": rdBase.rdkitVersion,
        },
        "future_outcome_rows_read": 0,
        "outcome_access_authorized": False,
        "rerun_allowed": False,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
