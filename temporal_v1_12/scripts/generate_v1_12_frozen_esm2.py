#!/usr/bin/env python3
"""Validate ESM2 runtime equivalence and embed only v1.12 blind proteins."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch


MAX_ABSOLUTE_DIFFERENCE = 1e-5
MINIMUM_COSINE = 0.999999


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    package_root = Path(__file__).resolve().parent.parent
    project_root = package_root.parents[1]
    parser.add_argument("--package-root", type=Path, default=package_root)
    parser.add_argument("--project-root", type=Path, default=project_root)
    parser.add_argument(
        "--blind-structures",
        type=Path,
        default=package_root / "source_frozen" / "V1_12_BLIND_STRUCTURES.tsv",
    )
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=package_root / "source_frozen" / "V1_12_SOURCE_MANIFEST.json",
    )
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument(
        "--output-dir", type=Path, default=package_root / "model_inputs"
    )
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    embedding_path = output_dir / "V1_12_FROZEN_ESM2_EMBEDDINGS.npz"
    index_path = output_dir / "V1_12_FROZEN_ESM2_INDEX.tsv"
    equivalence_path = output_dir / "V1_12_ESM2_RUNTIME_EQUIVALENCE.json"
    receipt_path = output_dir / "V1_12_FROZEN_ESM2_RECEIPT.json"
    if any(
        path.exists()
        for path in (embedding_path, index_path, equivalence_path, receipt_path)
    ):
        raise RuntimeError("v1.12 ESM2 outputs already exist; rerun prohibited")

    manifest = json.loads(args.source_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "STRUCTURE_ONLY_MEMBERSHIP_FROZEN":
        raise RuntimeError("Source structure-only gate has not passed")
    if manifest.get("blind_structure_table_sha256") != sha256_file(
        args.blind_structures
    ):
        raise RuntimeError("Blind structure hash differs from source manifest")

    v18_runtime = (
        project_root
        / "MASTER_SYNTHESIS_20260730"
        / "science_advances_prospective_v1_8_20260802"
        / "runtime"
        / "python_packages"
    )
    sys.path.insert(0, str(v18_runtime))
    sys.path.insert(0, str(project_root / "RareMol-AI_platform"))
    import huggingface_hub  # noqa: PLC0415
    import rdkit  # noqa: PLC0415
    import safetensors  # noqa: PLC0415
    import sklearn  # noqa: PLC0415
    import transformers  # noqa: PLC0415
    from scripts.generate_target_aware_esm2_embeddings import (  # noqa: PLC0415
        EXPECTED_HIDDEN_SIZE,
        MODEL_ID,
        MODEL_REVISION,
        WINDOW_RESIDUES,
        embed,
        model_file_audit,
        source_proteins,
    )
    from model_services.multimodal.target_aware_baseline import (  # noqa: PLC0415
        clean_protein,
        sha256_text,
    )

    model_dir = (
        project_root
        / "RareMol-AI_platform"
        / "model_services"
        / "multimodal"
        / "pretrained"
        / "facebook_esm2_t6_8M_UR50D_c731040"
    )
    model_files = model_file_audit(model_dir)
    reference_input = (
        project_root
        / "RareMol-AI_platform"
        / "data"
        / "external_sources"
        / "chembl37_exact_kd_blind_20260729"
        / "CHEMBL37_EXACT_KD_BLIND_PAIRS.csv"
    )
    reference_embedding = (
        project_root
        / "RareMol-AI_platform"
        / "validation"
        / "target_aware_phase16_chembl37_blind_20260729"
        / "ESM2_PROTEIN_EMBEDDINGS.npz"
    )
    reference_frame = pd.read_csv(
        reference_input,
        usecols=["protein_sha256", "protein_sequence"],
        dtype=str,
        keep_default_na=False,
    )
    reference_hash_matches = reference_frame.apply(
        lambda row: bool(clean_protein(row["protein_sequence"]))
        and sha256_text(clean_protein(row["protein_sequence"]))
        == row["protein_sha256"],
        axis=1,
    )
    reference_mismatches = sorted(
        reference_frame.loc[~reference_hash_matches, "protein_sha256"].unique()
    )
    reference = source_proteins(reference_frame.loc[reference_hash_matches].copy())
    hash_sorted = reference.sort_values("protein_sha256")
    length_sorted = reference.sort_values(
        ["normalized_length", "protein_sha256"]
    )
    selected_hashes = set(hash_sorted.head(8)["protein_sha256"])
    selected_hashes.update(hash_sorted.tail(8)["protein_sha256"])
    selected_hashes.update(length_sorted.head(8)["protein_sha256"])
    selected_hashes.update(length_sorted.tail(8)["protein_sha256"])
    selected = (
        reference.loc[reference["protein_sha256"].isin(selected_hashes)]
        .sort_values("protein_sha256")
        .reset_index(drop=True)
    )
    observed, observed_chunks = embed(selected, model_dir, args.batch_size)
    with np.load(reference_embedding, allow_pickle=False) as archive:
        frozen_hashes = archive["protein_sha256"].astype(str)
        frozen_embeddings = archive["embedding"].astype(np.float32)
    frozen_lookup = {
        value: frozen_embeddings[index]
        for index, value in enumerate(frozen_hashes)
    }
    expected = np.stack(
        [frozen_lookup[value] for value in selected["protein_sha256"]]
    ).astype(np.float32)
    absolute = np.abs(observed - expected)
    cosine = np.sum(observed * expected, axis=1) / (
        np.linalg.norm(observed, axis=1) * np.linalg.norm(expected, axis=1)
    )
    maximum_absolute = float(absolute.max())
    minimum_cosine = float(cosine.min())
    equivalence_passed = bool(
        maximum_absolute <= MAX_ABSOLUTE_DIFFERENCE
        and minimum_cosine >= MINIMUM_COSINE
    )
    runtime = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "torch": torch.__version__,
        "transformers": transformers.__version__,
        "huggingface_hub": huggingface_hub.__version__,
        "safetensors": safetensors.__version__,
        "rdkit": rdkit.__version__,
        "scikit_learn": sklearn.__version__,
        "deterministic_algorithms": True,
    }
    equivalence = {
        "schema_version": "science_advances_v1_12_esm2_equivalence_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "PASSED" if equivalence_passed else "FAILED",
        "reference_input_sha256": sha256_file(reference_input),
        "reference_embedding_sha256": sha256_file(reference_embedding),
        "reference_hash_consistent_pool": len(reference),
        "reference_hash_mismatch_exclusions": len(reference_mismatches),
        "reference_hash_mismatch_sha256": hashlib.sha256(
            "|".join(reference_mismatches).encode("utf-8")
        ).hexdigest(),
        "reference_proteins": len(selected),
        "reference_windows": int(observed_chunks.sum()),
        "observed_maximum_absolute_difference": maximum_absolute,
        "required_maximum_absolute_difference": MAX_ABSOLUTE_DIFFERENCE,
        "observed_minimum_cosine": minimum_cosine,
        "required_minimum_cosine": MINIMUM_COSINE,
        "runtime": runtime,
        "candidate_embeddings_generated": False,
        "outcome_rows_read": 0,
    }
    equivalence_path.write_text(
        json.dumps(equivalence, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if not equivalence_passed:
        print(json.dumps(equivalence, indent=2, sort_keys=True))
        return 2

    blind = pd.read_csv(
        args.blind_structures, sep="\t", dtype=str, keep_default_na=False
    )
    proteins = source_proteins(blind)
    embeddings, chunk_counts = embed(proteins, model_dir, args.batch_size)
    np.savez_compressed(
        embedding_path,
        protein_sha256=proteins["protein_sha256"].to_numpy(dtype=str),
        normalized_length=proteins["normalized_length"].to_numpy(dtype=np.int32),
        chunk_count=chunk_counts,
        embedding=embeddings,
    )
    pd.DataFrame(
        {
            "protein_sha256": proteins["protein_sha256"],
            "normalized_length": proteins["normalized_length"],
            "chunk_count": chunk_counts,
            "embedding_l2_norm": np.linalg.norm(embeddings, axis=1),
        }
    ).to_csv(index_path, sep="\t", index=False)
    receipt = {
        "schema_version": "science_advances_v1_12_frozen_esm2_receipt_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "OUTCOME_FREE_FROZEN_ESM2_READY",
        "source_manifest_sha256": sha256_file(args.source_manifest),
        "blind_structure_table_sha256": sha256_file(args.blind_structures),
        "unique_proteins": len(proteins),
        "embedding_shape": list(embeddings.shape),
        "model": {
            "id": MODEL_ID,
            "revision": MODEL_REVISION,
            "files": model_files,
            "hidden_size": EXPECTED_HIDDEN_SIZE,
        },
        "embedding_policy": {
            "window_residues": WINDOW_RESIDUES,
            "pooling": (
                "last hidden-state residue mean per window; residue-count "
                "weighted across windows; final L2 normalization"
            ),
            "dtype": "float32",
            "device": "cpu",
            "batch_size": args.batch_size,
        },
        "runtime_equivalence_path": str(equivalence_path),
        "runtime_equivalence_sha256": sha256_file(equivalence_path),
        "runtime_equivalence_passed": True,
        "embedding_path": str(embedding_path),
        "embedding_sha256": sha256_file(embedding_path),
        "index_path": str(index_path),
        "index_sha256": sha256_file(index_path),
        "runtime": runtime,
        "embedding_code_sha256": sha256_file(Path(__file__).resolve()),
        "outcome_fields_read": [],
        "outcome_rows_read": 0,
        "model_predictions_generated": False,
        "rerun_allowed": False,
    }
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
