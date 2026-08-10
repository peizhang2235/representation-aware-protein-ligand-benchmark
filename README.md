# Hierarchy- and representation-aware protein-ligand benchmark validation

Status: **public-release candidate; DOI pending external archive deposit**

Public code repository: https://github.com/peizhang2235/representation-aware-protein-ligand-benchmark

This repository supports inspection of the hierarchical and temporal analyses
reported in the Journal of Cheminformatics manuscript. It contains analysis code,
locked protocols and receipts, frozen predictions, rights-bounded result
tables, figure-source tables, publication figures, and checksum records.

The archive does not imply that every upstream database record was checked
against its original document. Two planned 100-record source-review samples
remain unscored (0/100 in each sample). They did not contribute to record
selection, correction, exclusion, or evidence. Residual source-reconstruction
error remains a study limitation.

## Directory map

- `retrospective_v1_7/`: retrospective aggregate results, ChEMBL37-derived
  public rows, figure code, figure sources, figures, protocols, and environment
  records;
- `temporal_v1_12/`: frozen temporal membership, predictions, locked aggregate
  results, sensitivity outputs, figure code, figures, locks, and receipts;
- `journal_outputs/`: the six main figure PDFs and machine-readable
  supplementary tables used for the Journal of Cheminformatics manuscript;
- `release_metadata/`: data-rights boundary, dependency inventory, code-license
  license, source-review status, and reproducibility instructions;
- `PUBLIC_MANIFEST_SHA256.tsv`: file sizes and SHA-256 hashes. The manifest
  excludes itself and `STAGING_AUDIT.json`.

The clean-directory reconstruction report is intentionally distributed beside
the ZIP rather than embedded in the ZIP it audits, so its recorded archive hash
can refer to the final tested object without a self-referential change.

## Deliberate exclusions

- BioLiP/MOAD row-level reconstruction and pair-level predictions;
- quarantined pre-unblind outcome files;
- private source-review queues, notes, identities, or decisions;
- copyrighted article or patent full text;
- local credentials, caches, absolute private paths, and author personal data;
- pretrained model files and frozen model weights pending upstream and
  rights-holder release review;
- v1.14 work in progress.

## License and DOI boundary

Original project software is released under the MIT License in `LICENSE`
(`SPDX-License-Identifier: MIT`).
Database-derived files retain their source-specific terms. BioLiP/MOAD row-level
reconstruction, private source-review materials, and model weights are excluded.
The deposited archive must match the checksum manifest before a repository DOI
is cited. A persistent DOI is intentionally not claimed until the external
archive deposit has completed and its landing page resolves; the DOI will be
added at the latest before the first revision.

This archive is a research-release decision record, not legal advice.
