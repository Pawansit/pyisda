# pyisda

A local, pip-installable Python toolkit for protein structural and mutational analysis built
around the **IBDC ISDA REST API** (`https://ibdc.dbt.gov.in/isda/api/documentation/`).

> **Domain migration note:** IBDC has moved from `ibdc.dbtindia.gov.in` to
> `ibdc.dbt.gov.in` and asks users to switch to the new domain. 


- Protein/sequence lookup (`pyisda.protein`)
- Structural (PDB) coverage details + ChimeraX mutation scripts (`pyisda.structure`)
- UniProt ⇄ PDB residue-numbering mapping (`pyisda.mapper`)
- Clinical + computational (AlphaMissense-style) mutation retrieval and analysis (`pyisda.mutation`)
- Per-residue SASA calculation via biotite (`pyisda.sasa`)

---

## Installation

Requires Python ≥ 3.9.

```bash
# From the extracted/cloned package directory:
cd pyisda
pip install -e .
```

This installs the core dependencies (`requests`, `pandas`, `numpy`,
`matplotlib`). The SASA module depends on `biotite`, which is kept optional
since it's heavier and not everyone needs it:

```bash
pip install -e ".[sasa]"
# or just:
pip install biotite
```

To install without `pyproject.toml`/editable mode (e.g. into a plain venv):

```bash
pip install -r requirements.txt
pip install .
```

---

## Quick start

```python
import pyisda as ibdc

# Protein-level metadata
details = ibdc.get_protein_details("P00533")   # EGFR
seq = ibdc.get_sequence("P00533")

# Structural coverage
coverage = ibdc.get_structure_details("P00533", include_pdb_records=True)

# Entry-level PDB summary
pdb_info = ibdc.get_pdb_summary("6q0j", additional_outputs="micromolecular_data")

# UniProt <-> PDB residue mapping, filtered to one chain directly
residue_map = ibdc.get_residue_map("P00533", "1m17", auth_chain_id="A")

# Mutation retrieval + analysis
clinvar_df = ibdc.get_mutation_table("P00533", source="clinvar")
annotated = ibdc.analyze_mutation_properties(clinvar_df)
summary = ibdc.get_mutation_summary(clinvar_df)

# Calculate per residue Solvent Accessible Surface Area (SASA)
result = ibdc.calculate_sasa_for_chain("6gel", "A")
```

### Worked example: focusing on a binding-site residue set

Combine experimental (ClinVar) and computational (AlphaMissense-style)
mutation sources into one table, then narrow it down to a specific
residue set — e.g. a ligand-binding site:

```python
import pyisda as ibdc

uniprot_id = "P00533"

# 1. Structural coverage -> pick the best-covered PDB + its chain
coverage = ibdc.get_structure_details(uniprot_id, include_pdb_records=True)
best_pdb = coverage["highest_coverage_pdb_id"]
chain_id = coverage["highest_coverage_auth_chain"]

# 2. Residue map, filtered to that chain directly
residue_map = ibdc.get_residue_map(uniprot_id, best_pdb, auth_chain_id=chain_id)

# 3. Experimental + computational mutation sources
clinvar_df = ibdc.get_mutation_table(uniprot_id, source="clinvar", record_type="full")
annotated_mutations = ibdc.analyze_mutation_properties(clinvar_df)
comp_predictions = ibdc.get_computational_mutations(uniprot_id)

# 4. Attach a computational pathogenicity prediction to each experimental
#    mutation, matched on the exact substitution (position + alt AA)
combined = ibdc.merge_experimental_with_computational(annotated_mutations, comp_predictions)

# 5. Narrow down to a residue set of interest, e.g. ligand-binding-site residues
binding_site_residues = [790, 793, 797, 855, 858]  # example UniProt positions
focus_set = ibdc.subset_by_positions(combined, binding_site_residues)
```

---

## Module reference

### `pyisda.protein`

| Function | Description |
|---|---|
| `get_protein_details(uniprot_id)` | Organism, gene name(s)/synonyms, structural record count, sequence length |
| `get_sequence(uniprot_id)` | Raw sequence string, length, molecular weight |
| `mutate_sequence(sequence_record, mutations)` | Apply a `Position`/`Alt_AA` (3-letter code) mutation table to a sequence, returns the mutated sequence string |


```python
seq_record = ibdc.get_sequence("P00533")
muts = pd.DataFrame({"Alt_AA": ["Gly"], "Position": [12]})
mutated_seq = ibdc.mutate_sequence(seq_record, muts)
```

### `pyisda.structure`

| Function | Description |
|---|---|
| `get_structure_details(uniprot_id, include_pdb_records=False)` | PDB/chain counts and the highest-coverage structure; pass `include_pdb_records=True` for the full per-PDB record list |
| `get_pdb_summary(pdb_id, additional_outputs=None)` | Entry-level summary for a single PDB ID (`{ISDA_BASE_URL}/pdb_summary/<pdb_id>/`)`.`additional_outputs`  also accepts a list, joined with commas and cover the following `experiment, publication, micromolecular_data, interactions` |
| `generate_mutated_structure_script(mutation_table, pdb_id, auth_chain_id, output_dir=".")` | Writes a ChimeraX `.cxc` script that opens `<pdb_id>.cif`, applies `swapaa` mutations for the given chain, and saves `<pdb_id>_mutated_structure.cif` |

### `pyisda.mapper`

| Function | Description |
|---|---|
| `get_residue_map(uniprot_id, pdb_id, auth_chain_id=None)` | List of `{unp_residue, pdb_residue, pdb_auth_chain, pdb_chain}` dicts mapping UniProt numbering to PDB numbering for one PDB entry. Pass `auth_chain_id` (e.g. `"A"`) to filter to that chain directly — records for other chains are skipped before the map is even built, equivalent to but more efficient than fetching everything and then `df[df['pdb_auth_chain'] == chain_id]` |

### `pyisda.mutation`

| Function | Description |
|---|---|
| `get_mutation_table(uniprot_id, source="clinvar", record_type="default", selection=None)` | Clinical mutation records with `Ref_AA`/`Position`/`Alt_AA` parsed from `Variant AA Change` |
| `get_computational_mutations(uniprot_id, ranges=None, significance_type=None)` | AlphaMissense-style computational predictions (replaces the two divergent `AFMutationTable` copies) |
| `filter_by_significance(df, pattern, case=False)` | Regex filter on `Clinical Significance` |
| `analyze_mutation_properties(df, ...)` | Adds hydrophobic/hydrophilic property columns + change-type label |
| `get_mutation_summary(df, top_n=10)` | Count, significance distribution, consequence-type distribution, top mutation hotspots |
| `merge_mutations_with_pdb_map(mutation_df, residue_map_df)` | Joins a mutation table to `get_residue_map` output on residue position, to attach PDB-numbered residue coordinates. Works even when `residue_map_df` is a filtered slice (e.g. one chain). **Not** for joining two mutation tables to each other — pass a `get_residue_map`-shaped DataFrame (must contain `unp_residue`) or it raises a clear `ValueError` pointing you at `merge_experimental_with_computational` instead |
| `merge_experimental_with_computational(experimental_df, computational_df, position_col="Position", alt_aa_col="Alt_AA", how="left")` | Attaches AlphaMissense-style computational pathogenicity predictions to each experimental/clinical mutation, joined on the exact substitution (position + alt AA). Handles the 1-letter vs. 3-letter `Alt_AA` mismatch between `get_mutation_table` and `get_computational_mutations` automatically |
| `subset_by_positions(df, positions, position_col="Position")` | Filters any mutation/merged table down to a specific set of residue positions — e.g. a ligand-binding-site residue list — to focus downstream analysis on that region |
| `plot_mutation_matrix(df, ax=None)` | Position × Alt_AA heatmap of mutation counts |

### `pyisda.sasa`

| Function | Description |
|---|---|
| `calculate_sasa_for_chain(pdb_id, chain_id, probe_radius=1.4, local_pdb_path=None, source="isda")` | Per-residue SASA via biotite. By default (`source="isda"`) downloads the structure as `.cif` from the IBDC ISDA download endpoint (`{ISDA_BASE_URL}/download.<pdb_id>.cif`), which replaces the previous RCSB-only fetch; pass `source="rcsb"` to fall back to fetching a `.pdb` file from RCSB instead, or `local_pdb_path` to use a file you already have (`.cif` or `.pdb`, detected by extension). Returns a `SASAResult(residue_ids, residue_names, sasa_per_residue)` namedtuple, or `None` on error |
| `fetch_structure_cif(pdb_id, output_dir=".")` | Just the download step: fetches `{ISDA_BASE_URL}/download.<pdb_id>.cif` and saves it locally, returning the file path. Used internally by `calculate_sasa_for_chain`, but also handy standalone (e.g. as the input structure for `generate_mutated_structure_script`) |

### `pyisda.visualize`

Turns the outputs of the other modules into residue-position tracks —
the kind of coverage/lollipop/profile views you'd see on UniProt or PDBe.

| Function | Description |
|---|---|
| `plot_structural_coverage(pdb_records, sequence_length=None, max_entries=None, ax=None)` | Genome-browser-style track: one horizontal bar per (PDB, chain) spanning its `UNIPROT_START`..`UNIPROT_END` range. Input is `get_structure_details(uniprot_id, include_pdb_records=True)["pdb_records"]` |
| `plot_mutation_lollipop(mutation_df, position_col="Position", significance_col="Clinical Significance", sequence_length=None, ax=None)` | Stem plot of mutation counts per position, colored red where any record at that position is flagged pathogenic. Input is `get_mutation_table(...)` or `get_computational_mutations(...)` |
| `plot_sasa_profile(sasa_result, sequence_length=None, ax=None)` | Line/area plot of per-residue SASA. Input is `calculate_sasa_for_chain(...)` |
| `plot_protein_overview(sequence_length=None, pdb_records=None, mutation_df=None, sasa_result=None, title=None)` | Stacks any subset of the three plots above on one shared residue-position x-axis — a single combined figure |

```python
import pyisda as ibdc

uniprot_id = "P00533"  # EGFR
seq = ibdc.get_sequence(uniprot_id)
coverage = ibdc.get_structure_details(uniprot_id, include_pdb_records=True)
mutations = ibdc.get_mutation_table(uniprot_id)

fig = ibdc.plot_protein_overview(
    sequence_length=seq["length"],
    pdb_records=coverage["pdb_records"],
    mutation_df=mutations,
    title=f"{uniprot_id} structural + mutation overview",
)
fig.savefig("egfr_overview.png", dpi=150)
```

Each function also works standalone and accepts an `ax=` so you can drop
it into your own subplot grid instead of the combined view:

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 4))
ibdc.plot_structural_coverage(coverage["pdb_records"], sequence_length=seq["length"], ax=ax)
```

**Numbering caveat:** structural coverage and clinical mutation positions
are UniProt-numbered, while a SASA profile from
`calculate_sasa_for_chain` is numbered per the PDB chain itself — these
usually track closely but aren't guaranteed identical. For an
exact residue-by-residue overlay, translate one numbering scheme to the
other first with `get_residue_map`.

### `pyisda.interactive` — HTML/Plotly versions

Same data inputs as `visualize`, rendered as zoomable/pannable HTML with
hover tooltips instead of a static image. Reach for these when a protein
has many structures or mutations and the static plots' y-axis labels or
coverage ranges get too cramped to read — you can zoom into a region,
pan around, and hover any bar/point/marker for its exact values instead
of squinting at overlapping labels.

Requires the `interactive` extra:
```bash
pip install -e ".[interactive]"   # installs plotly
```

| Function | Description |
|---|---|
| `plot_structural_coverage_html(pdb_records, sequence_length=None, title=None)` | Same coverage track as `plot_structural_coverage`, but zoomable/pannable with an x-axis range slider and hover tooltips giving the exact PDB ID, chain, and UniProt range |
| `plot_mutation_lollipop_html(mutation_df, ...)` | Same as `plot_mutation_lollipop`; hovering a point shows position, mutation count, and the set of significance labels at that position |
| `plot_sasa_profile_html(sasa_result, ...)` | Same as `plot_sasa_profile`; hover shows exact residue id, name, and SASA value |
| `plot_protein_overview_html(...)` | Combined figure with synced x-axes across panels — zoom/pan one panel and the others follow |
| `save_html(fig, path, standalone=True)` | Writes any of the above to a self-contained `.html` file. `standalone=True` (default) embeds the Plotly JS so the file works fully offline (~4-5 MB); `standalone=False` loads it from a CDN instead (much smaller file, needs internet to view) |

```python
import pyisda as ibdc

uniprot_id = "P00533"
seq = ibdc.get_sequence(uniprot_id)
coverage = ibdc.get_structure_details(uniprot_id, include_pdb_records=True)
mutations = ibdc.get_mutation_table(uniprot_id)

fig = ibdc.plot_protein_overview_html(
    sequence_length=seq["length"],
    pdb_records=coverage["pdb_records"],
    mutation_df=mutations,
    title=f"{uniprot_id} interactive overview",
)
ibdc.save_html(fig, "egfr_overview.html")   # open in any browser
# or, in a notebook:
# fig.show()
```

If `plotly` isn't installed, importing `pyisda` still works fine —
only calling one of the `*_html` functions raises a clear `ImportError`
telling you to install the `interactive` extra.

---

## Configuration

The ISDA base URL defaults to `https://ibdc.dbt.gov.in/isda/api` (the
current domain as of July 2026) and is centralized and overridable, e.g. to
point at a staging environment or roll back to the old host:

```python
import pyisda as ibdc
ibdc.config.ISDA_BASE_URL = "https://staging.example.org/isda/api"
```

## Logging

All modules log through the `"pyisda"` logger instead of `print()`.
Enable it in your own scripts/notebooks with:

```python
import logging
logging.basicConfig(level=logging.INFO)  # or logging.WARNING
```

By default (no `basicConfig` call), the package stays silent — it attaches
a `NullHandler` so it won't print anything unless you opt in.

## Error handling

Every network-backed function follows the same contract: on request
failure (timeout, connection error, HTTP error, bad JSON) it logs a warning
and returns `None` — it never raises. Internally, request failures are
raised as `pyisda.ISDARequestError` and caught at each public
function's boundary, so if you're extending the package yourself you can
catch that exception directly instead of parsing log output.

---

## Testing

```bash
pip install pytest
pytest tests/
```

`tests/test_offline.py` covers all pure-Python logic (sequence mutation,
mutation property analysis, summary stats, residue-map merging, ChimeraX
script generation) plus URL-construction checks for the network-backed
functions (via monkeypatching, no real request made). Functions that call
the live ISDA API (`get_protein_details`, `get_structure_details`,
`get_pdb_summary`, `get_residue_map`, `get_mutation_table`,
`get_computational_mutations`, `calculate_sasa_for_chain`) aren't covered
end-to-end here since that requires network access to `ibdc.dbt.gov.in`;
verify those against the live API in an environment with access to it.

---

## Project layout

```
pyisda/
├── pyproject.toml
├── requirements.txt
├── README.md
├── LICENSE
├── .gitignore
├── pyisda/
│   ├── __init__.py        # public API
│   ├── _client.py         # shared HTTP helper, ISDA_BASE_URL, logger
│   ├── protein.py
│   ├── structure.py
│   ├── mapper.py
│   ├── mutation.py
│   ├── sasa.py
│   ├── visualize.py
│   └── interactive.py
└── tests/
    └── test_offline.py
```

## License

MIT — see [LICENSE](LICENSE).
