# V1.13 locked protocol

## Why this source

GtoPdb 2026.2 is available now, is expert curated outside the BindingDB and
ChEMBL release process, exposes pKd directly, and contains ligand structures,
target identities, species, and source references. It can therefore reproduce
the same source-wide versus within-target contrast without waiting for a future
ChEMBL release.

The database name alone is not sufficient evidence of independence. A paper or
measurement can be curated into several databases. The confirmatory claim is
therefore restricted to rows that are both exact-pair disjoint and stable-source-
reference disjoint from every dataset whose outcomes or predictions were
previously inspected in this project.

## Locked sequence

1. Hash the source registration and this protocol.
2. Acquire the fixed GtoPdb 2026.2 payload into outcome quarantine.
3. Export only interaction identifiers, ligand identifiers, target identifiers,
   species, source-reference identifiers, ligand structures, and protein
   sequences. No affinity value or derived label may enter this export.
4. Construct prior pair and source-document ledgers from already exposed data.
5. Apply the pair, document, entity, structure, sequence, support, and
   concentration gates without computing any model-performance statistic.
6. Freeze the blind membership table and its hash.
7. Score the table once with all frozen v1.8 models and freeze predictions.
8. Join outcomes once by blind pair identifier, run the locked analysis, and
   retain whichever decision branch is reached.

## Scientific endpoint

The primary contrast is not absolute error. It is the difference between the
source-wide Spearman association and target-equal within-target Spearman
association. This directly tests whether models can exploit broad target/ligand
scale differences while having weaker medicinal-chemistry ranking within a
target series.

Rows from targets with fewer than eight eligible compounds or fewer than three
Bemis-Murcko scaffolds are not confirmatory. The feasibility gate requires at
least eight eligible targets, 100 rows, and 20 resolved independent source
references after all firewalls.

## Claim boundary

This package may support an independent-curation, reference-disjoint replication.
It cannot be called a prospective external validation, and it cannot establish
that GtoPdb as a whole is independent of BindingDB or ChEMBL.

