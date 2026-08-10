# Source-specific data-rights matrix

Status: **conservative release boundary approved by the author/rights-holder representative on 2026-08-07**

| Source or material | Public basis checked 2026-08-04 | Candidate treatment | Before deposit |
|---|---|---|---|
| BindingDB-curated records | The 2024 BindingDB resource article states that data curated by BindingDB are shared under CC BY 4.0, with citation and change notation. | Include permitted BindingDB-derived temporal rows, transformations, provenance, and hashes with attribution. | Confirm each released row's origin flag and retain attribution/change notice. |
| ChEMBL-derived records | The official ChEMBL release page states CC BY-SA 3.0 Unported. The BindingDB resource article states that ChEMBL-curated records in BindingDB retain those terms. | Include the ChEMBL37-derived public ledger under source terms; do not apply the software license to it. | Add file-level CC BY-SA 3.0 metadata and the ChEMBL release/citation. |
| BioLiP | The official download page says the data are freely available, but no explicit dataset license was located. The page records removal of 24,809 PDBbind-CN affinity records in January 2025 because of a licensing issue. | Exclude row-level reconstruction and pair-level predictions; include nonrestricted aggregates, hashes, lineage, and code only. | Approved conservative boundary: no row-level BioLiP material is deposited. |
| Binding MOAD | Resource articles describe access, but no current explicit database-license statement was located and the original site has been sunset. An article license does not automatically license the database. | Exclude row-level reconstruction and pair-level predictions; include nonrestricted aggregates, hashes, lineage, and code only. | Approved conservative boundary: no row-level Binding MOAD material is deposited. |
| Primary literature and patent text | Copyright remains with the relevant rights holders. | Do not include full text or substantial excerpts; include identifiers, citations, and derived summaries only. | Confirm that no private notes or source text entered the archive. |
| Original analysis code | MIT selected for original project software; upstream dependencies retain their own licenses. | Include original code with the MIT LICENSE; do not apply MIT to data, publication text, weights, or dependencies. | Approved by the author/rights-holder representative; retain dependency inventory. |
| Model weights and pretrained artifacts | Upstream and author-generated artifacts may have distinct terms. | Exclude weights; include frozen predictions, receipts, and artifact hashes. | Decide whether author-generated weights may be deposited and record upstream model terms. |

Primary pages:

- BindingDB resource article: https://www.bindingdb.org/rwd/bind/gkae1075.pdf
- ChEMBL official release page: https://www.ebi.ac.uk/chembldb/
- BioLiP official download page: https://zhanggroup.org/BioLiP/download.html
- Binding MOAD sunset article: https://pmc.ncbi.nlm.nih.gov/articles/PMC9944886/

This matrix documents a conservative release boundary and is not legal advice.
