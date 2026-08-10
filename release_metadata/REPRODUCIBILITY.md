# Reproducibility levels and commands

Run commands from the extracted archive root.

## Level 1: checksum and file-integrity audit

```bash
python3 release_metadata/verify_manifest.py
```

This checks every file recorded in `PUBLIC_MANIFEST_SHA256.tsv`. The manifest
does not contain a DOI and does not authorize release.

## Level 2: rights-bounded retrospective figure reconstruction

```bash
Rscript retrospective_v1_7/scripts/build_v1_7_figures.R
```

The public-stage script uses aggregate fallback tables when the restricted
combined BioLiP/MOAD pair table is absent. It reconstructs the retrospective
figures without redistributing restricted rows.

## Level 3: temporal figure and figure-source reconstruction

```bash
python3 temporal_v1_12/scripts/build_v1_12_publication_figures.py \
  --package-root temporal_v1_12
```

This reconstructs the temporal figures and their source tables from frozen
membership metadata and locked aggregate results. It does not reopen or rerun
the outcome-unblinding step.

## Level 4: analysis and prediction boundaries

The v1.12 locked analysis cannot be rerun from this public candidate because
the original pre-unblind outcome files remain deliberately excluded. Frozen
predictions, locked results, receipts, and their hashes are included for
audit. Exact prediction regeneration additionally requires frozen model
weights and the pinned upstream ESM2 model, which are not included pending
rights-holder and upstream-license review.

The retrospective scientific-upgrade script likewise requires the
rights-controlled BioLiP/MOAD row-level reconstruction. Its code and aggregate
outputs are present, but the restricted rows are not.

## Recorded software

Successful source runs used R and Python. Package versions and upstream license
records are listed in `release_metadata/DEPENDENCY_LICENSE_INVENTORY.tsv`.
Figure-only reconstruction requires substantially fewer dependencies than
model prediction. Exact platform equivalence is documented in the frozen
receipts rather than inferred from the current machine.
