# ChEMBL fragment space analysis

UMAP of Morgan fingerprints for the fragment-like subset of ChEMBL, plus scaffold
representation, clustering and target analysis. Everything runs from `workflow.ipynb`.

Run against the **ChEMBL 37** SQLite release (2026-05-29). The release is read-only
throughout (`file:...?mode=ro`), lives at `Data/chembl_37/chembl_37_sqlite/chembl_37.db`
and is gitignored (~30 GB).

Environment: `conda env create -f env.yaml` (the kernel is named `general`).

## Layout

- `workflow.ipynb` — the pipeline.
- `Data/` — fragments file, UMAP coordinates, cluster and target CSVs, `diagnostics/`.
- `Artifacts/` — figures, SVG + PNG of each.
- `Saves/` — frozen runs (see below).
- `old_pipeline.py` — the script this notebook replaced, kept for reference.

## Notebook structure

1. **Imports and settings** — paths and the constants used everywhere below.
2. **Extract the fragments** — `QUERY` against `compound_properties`, so only
   fragment-like rows ever reach RDKit. Prints the subset size.
3. **Clean and write** — parse, drop the counter-ion (keeping the piece nearest
   `mw_freebase`), canonicalise, dedup on the canonical SMILES. Batched read,
   streaming write to `Data/chembl_37_fragments.csv`.
4. **Filter diagnostics** — read-only side checks, one CSV each in `Data/diagnostics/`:
   rows inside the MW limit rejected by another property, stereoisomers that end up
   as separate rows, and fragments too small to mean anything.
5. **Morgan fingerprints** — bit matrix, with unparsable rows dropped from the frame.
6. **Fit UMAP** — Jaccard metric, coordinates saved to `Data/chembl_37_fragments_umap.csv`.
7. **Plots** — density (hexbin, log scale) and scatter of the whole map.
8. **Scaffold representation** — Murcko scaffold per fragment, then one end of the
   representation ranking highlighted in red over the map.
9. **Clusters within the selection** — HDBSCAN on the selected points, or one cluster
   per scaffold; the largest few coloured, the rest of the selection grey. Saves the
   per-fragment cluster table and a grid of example molecules.
10. **Targets** — maps the selection back through `activities` to targets and ChEMBL
    protein classes, two histograms, three CSVs.

## Parameters

Cell 3 (settings):

- `DB`, `FRAGMENTS`, `EMBEDDING`, `FIGURE_STEM` — input and output paths.
- `BATCH` — rows per fetch while streaming the release.
- `FLATTEN` — strip stereochemistry before dedup, so stereoisomers collapse into one
  row instead of one row (and one UMAP point) each. The stereo diagnostic still runs
  and finds nothing.
- `RADIUS`, `N_BITS` — Morgan fingerprint.
- `N_NEIGHBORS`, `MIN_DIST`, `N_JOBS` — UMAP. No `random_state`, since seeding forces
  UMAP single-threaded.

`QUERY` (cell 5) — the property filter itself; edit it to change the subset. MW is the
freebase weight, so salts are judged on the parent. The diagnostics read their
thresholds back out of `QUERY`, so they cannot drift from it.

Scaffold representation:

- `SCAFFOLD_END` — `"bottom"` for the least represented scaffolds, `"top"` for the most.
- `SCAFFOLD_NUMBER` — take exactly this many scaffolds; `None` to cut by fraction.
- `SCAFFOLD_FRACTION` — take n% (`0.05`, `0.01`, ...); `None` when cutting by count.
  Setting both raises.
- `SCAFFOLD_BASIS` — n% *of what*: `"fragments"` or `"scaffolds"`. Ignored (with a
  warning) when `SCAFFOLD_NUMBER` is set. At the top end use `"scaffolds"`, since
  benzene alone carries more than 5% of the fragments.

Clustering:

- `MIN_CLUSTER_SIZE`, `MIN_SAMPLES` — HDBSCAN.
- `CLUSTER_SELECTION_METHOD` — `"leaf"` (tight local clusters), `"eom"` (fewer, much
  larger), or `"scaffold"` (no HDBSCAN, one cluster per scaffold in the cut; warns).
- `N_CLUSTERS_SHOWN`, `CLUSTER_COLOURS` — how many clusters get a colour.

Targets:

- `N_TARGETS_SHOWN` — bars per histogram.
- `PROTEIN_CLASS_LEVEL` — class depth: 1 = Enzyme / Membrane receptor, 2 = Kinase /
  Family A GPCR, 3+ narrower. The cell walks `parent_id` up to this level rather than
  filtering on it, or whole families would drop out.
- `TARGET_TYPES` — e.g. `("SINGLE PROTEIN",)`; `None` keeps everything.
- `TARGET_LABEL` — bar labels: `"gene"` symbol, `"accession"` (UniProt) or `"name"`
  (ChEMBL pref_name). All-or-nothing: if any bar drawn lacks that label, the whole
  histogram falls back one level, gene -> accession -> name.

Output names carry the cut, so runs don't overwrite each other, e.g.
`Scaffolds_top5_...`, `Clusters_bottom5pct_...`, `Clusters_top5_by_scaffold_...`.

## Saved runs

Both under `Saves/`, each with `config.txt` (the exact settings used), `data/` and
`artefacts/` (split into `scaffold_top` and `scaffold_bottom`):

- `ro3` — rule of three: `mw_freebase <= 300`, `alogp <= 3`, `hbd <= 3`, `hba <= 3`, `rtb <= 3`.
- `ro3_ex_rtb` — the same, with rotatable bonds relaxed to `rtb <= 6`.

Each run covers two cuts: the top 5 scaffolds (clustered by scaffold) and the bottom
5% by fragments (HDBSCAN, `eom`).
