# Gene annotation review

Checks whether a protein that the pipeline cannot give one official gene can safely be assigned one, using
the DIA-NN peptides and current public records. It is read-only towards the pipeline: it reads DIA-NN outputs
and writes only to `--out`. Approved results go into `data/metadata/manual_gene_annotation_overrides.csv`
by hand.

## What it asks

For each accession:

1. **Which genes could it come from?** These are the Gene IDs that UniProt links to the accession, plus those the
   installed `org.Mm.eg.db` links to it (the source of the pipeline's gene annotation).
2. **Do the observed peptides tell those genes apart?** Every peptide DIA-NN reported is searched in every
   RefSeq protein of every candidate gene, with I and L treated as the same residue. With `--fasta`, it is also
   searched in every entry of a whole-proteome FASTA.
   - A peptide is **gene-unique** only if it occurs in one candidate gene's proteins and in no other gene's FASTA
     entry. Entries of the same gene under older names or aliases do not count against it.
   - Any other peptide is **shared**.

   DIA-NN's own `Proteotypic` flag is not used for this:
   - it only means unique among the search FASTA entries;
   - a copy gene without its own entry never counts against it;
   - an I/L variant in another gene does not count either.

   Always supply a whole-proteome FASTA. Without it, a peptide is only compared with the candidate genes. For
   example, Cbx3's IIGATDSSGELMFLMK occurs unchanged in Cbx1 (P83917, residues 122–137, same flanking
   residues). DIA-NN listed only Cbx3 for it, and only the proteome check reveals the match.
3. **How much of the measured amount is gene-specific?** The protein-group value in `report.pg_matrix.tsv` is a
   MaxLFQ over all precursors of the group. The script recomputes MaxLFQ from all, gene-unique-only and
   shared-only precursors. It reports:
   - the gene-unique share of the signal;
   - how closely each version follows the delivered value;
   - at animal level, how much including shared peptides shifts each group contrast.

Read the outcome per accession (`peptide_specificity_summary.csv`, column `verdict`) like this:

| Verdict | Meaning |
|---|---|
| one candidate gene; every peptide in its proteins | A database annotation gap; an override is defensible |
| all peptides shared | The gene of origin cannot be resolved from these data; leave the protein out of gene-level analyses |
| gene-unique and shared peptides | Identity supported, but the quantity is mixed; check step 4 before assigning |
| no candidate gene | No NCBI Gene ID; it cannot enter GO or EWCE |

## How to run

Run from the repository root:

```
Rscript tools/gene_annotation_review/01_extract_precursors.R --pr-matrix <report.pr_matrix.tsv> --accessions <accessions.csv> --out <dir>
Rscript tools/gene_annotation_review/02_fetch_records.R --out <dir> [--cache <dir>] [--extra-gene-ids <csv>]
Rscript tools/gene_annotation_review/03_peptide_gene_specificity.R --out <dir> [--fasta <search FASTA>[,<contaminants FASTA>]]
Rscript tools/gene_annotation_review/04_quantity_attribution.R --out <dir> --pr-matrix <...> --pg-matrix <...> --run-map data/metadata/pg_matrix_rerun_run_map.csv
```

- `accessions.csv` needs one column, `accession`.
- Step 2 sends only public accession numbers and Gene IDs. It caches every response and logs the URL and retrieval time
  in `fetch_log.csv`; delete the cache to refresh the records.
- For `--fasta`, give the FASTA the search used if you have it. Otherwise give the current complete UniProt mouse
  proteome (Swiss-Prot plus TrEMBL, with isoforms). Several files can be given separated by commas.
- Keep `--out` short: R on this host cannot open paths of 260 characters or more.

## Limits

- **Today's records:** RefSeq and UniProt are current records, not necessarily what the search FASTA held.
- **Approximate MaxLFQ:** the re-implementation follows DIA-NN's protein-group values closely (r about 0.93–0.995 for
  the 2026-11-16 check) but not exactly. DIA-NN's gene-unique quantity (`Genes.MaxLFQ.Unique` in the main report)
  is the better reference when available.
- **Animal-level tests:** these have little power with 3/3/2 animals. Read the bias sizes, not only the P values.
