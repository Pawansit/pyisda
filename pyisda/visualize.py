"""
Visualizations built on top of the other modules' outputs:

- `plot_structural_coverage`   — track of PDB coverage ranges along a
  sequence, from `structure.get_structure_details(..., include_pdb_records=True)["pdb_records"]`
- `plot_mutation_lollipop`     — mutation positions/counts along a
  sequence, from `mutation.get_mutation_table` / `get_computational_mutations`
- `plot_sasa_profile`          — per-residue SASA, from `sasa.calculate_sasa_for_chain`
- `plot_protein_overview`      — stacks the three above on a shared
  residue-position axis for a single combined view

Note on numbering: PDB structural coverage and clinical mutation
positions are typically UniProt-numbered, while a SASA profile from
`sasa.calculate_sasa_for_chain` is numbered according to the PDB chain
itself. These are usually close but not always identical — for a
residue-exact overlay, first translate one numbering scheme to the other
with `mapper.get_residue_map`.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ._client import logger

_DEFAULT_FIGSIZE_TRACK = (12, 4)


def plot_structural_coverage(
    pdb_records: List[Dict[str, Any]],
    sequence_length: Optional[int] = None,
    max_entries: Optional[int] = None,
    ax=None,
):
    """
    Plot a genome-browser-style track of PDB structural coverage along a
    UniProt sequence: one horizontal bar per (PDB, chain) record spanning
    its UNIPROT_START..UNIPROT_END range.

    Args:
        pdb_records: The `pdb_records` list from
            `structure.get_structure_details(uniprot_id, include_pdb_records=True)`.
            Each record is expected to have `PDB`, `UNIPROT_START`,
            `UNIPROT_END`, and (optionally) `PDB_CHAIN`/`AUTH_CHAIN`.
        sequence_length: Full sequence length, to fix the x-axis extent
            (e.g. from `protein.get_sequence(...)["length"]`). If omitted,
            the axis is scaled to the covered range only.
        max_entries: If set, only plot the top N records by coverage span
            (useful when a protein has hundreds of structures).
        ax: Optional matplotlib Axes to draw on. If None, creates a new
            figure and calls `plt.show()`.

    Returns:
        The matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    if not pdb_records:
        logger.warning("plot_structural_coverage: no pdb_records to plot.")
        if ax is None:
            _, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE_TRACK)
        ax.text(0.5, 0.5, "No structural coverage records", ha="center", va="center")
        return ax

    rows = []
    for record in pdb_records:
        try:
            start = int(record["UNIPROT_START"])
            end = int(record["UNIPROT_END"])
        except (KeyError, TypeError, ValueError):
            continue
        rows.append({
            "pdb_id": record.get("PDB", "?"),
            "chain": record.get("PDB_CHAIN") or record.get("AUTH_CHAIN") or "?",
            "start": start,
            "end": end,
            "span": end - start,
        })

    if not rows:
        logger.warning("plot_structural_coverage: no records had usable UNIPROT_START/END.")
        if ax is None:
            _, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE_TRACK)
        ax.text(0.5, 0.5, "No usable coverage records", ha="center", va="center")
        return ax

    df = pd.DataFrame(rows).sort_values("start").reset_index(drop=True)
    if max_entries:
        df = df.sort_values("span", ascending=False).head(max_entries).sort_values("start").reset_index(drop=True)

    created_fig = ax is None
    if created_fig:
        height = max(3, 0.25 * len(df) + 1)
        _, ax = plt.subplots(figsize=(_DEFAULT_FIGSIZE_TRACK[0], height))

    colors = plt.cm.tab20(np.linspace(0, 1, max(len(df), 1)))
    for i, row in df.iterrows():
        ax.broken_barh([(row["start"], row["span"] or 1)], (i - 0.4, 0.8), facecolors=colors[i % len(colors)])

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"{r.pdb_id}:{r.chain}" for r in df.itertuples()], fontsize=8)
    ax.set_xlabel("UniProt residue position")
    ax.set_title(f"Structural coverage ({len(df)} record{'s' if len(df) != 1 else ''})")
    if sequence_length:
        ax.set_xlim(0, sequence_length)
    ax.invert_yaxis()

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax


def plot_mutation_lollipop(
    mutation_df: pd.DataFrame,
    position_col: str = "Position",
    significance_col: Optional[str] = "Clinical Significance",
    pathogenic_pattern: str = "pathogenic",
    sequence_length: Optional[int] = None,
    ax=None,
):
    """
    Plot mutation counts per residue position as a lollipop/stem plot,
    optionally colored by whether the position has any pathogenic-flagged
    record.

    Args:
        mutation_df: Output of `mutation.get_mutation_table` or
            `mutation.get_computational_mutations` (must have a numeric
            `position_col`).
        position_col: Column holding the (1-indexed) residue position.
        significance_col: Column to check for a pathogenic flag, or None
            to skip coloring and plot all stems in one color.
        pathogenic_pattern: Case-insensitive substring marking a record as
            pathogenic (matched against `significance_col`).
        sequence_length: Full sequence length, to fix the x-axis extent.
        ax: Optional matplotlib Axes to draw on. If None, creates a new
            figure and calls `plt.show()`.

    Returns:
        The matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    if mutation_df is None or mutation_df.empty:
        logger.warning("plot_mutation_lollipop: mutation_df is empty or None.")
        if ax is None:
            _, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE_TRACK)
        ax.text(0.5, 0.5, "No mutation records", ha="center", va="center")
        return ax

    df = mutation_df.copy()
    df[position_col] = pd.to_numeric(df[position_col], errors="coerce")
    df = df.dropna(subset=[position_col])

    if significance_col and significance_col in df.columns:
        is_pathogenic = df.groupby(position_col)[significance_col].apply(
            lambda vals: vals.str.contains(pathogenic_pattern, case=False, na=False).any()
        )
        counts = df.groupby(position_col).size()
        colors = ["#c0392b" if is_pathogenic.get(pos, False) else "#2980b9" for pos in counts.index]
    else:
        counts = df.groupby(position_col).size()
        colors = "#2980b9"

    created_fig = ax is None
    if created_fig:
        _, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE_TRACK)

    markerline, stemlines, baseline = ax.stem(counts.index, counts.values, basefmt=" ")
    plt.setp(stemlines, color="#95a5a6", linewidth=1)
    plt.setp(markerline, marker="o", markersize=5)
    if isinstance(colors, list):
        markerline.set_markerfacecolor("none")
        ax.scatter(counts.index, counts.values, c=colors, s=25, zorder=3)
    else:
        markerline.set_markerfacecolor(colors)

    ax.set_xlabel("Residue position")
    ax.set_ylabel("Mutation count")
    ax.set_title(f"Mutation positions ({int(counts.sum())} total)")
    if sequence_length:
        ax.set_xlim(0, sequence_length)

    if significance_col and significance_col in df.columns:
        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#c0392b", markersize=6, label="Pathogenic"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#2980b9", markersize=6, label="Other"),
        ]
        ax.legend(handles=legend_handles, loc="upper right", fontsize=8)

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax


def plot_sasa_profile(sasa_result, sequence_length: Optional[int] = None, ax=None):
    """
    Plot per-residue SASA as a line/area profile.

    Args:
        sasa_result: A `sasa.SASAResult` (or any object with
            `residue_ids` and `sasa_per_residue` attributes/items) from
            `sasa.calculate_sasa_for_chain`.
        sequence_length: Full sequence length, to fix the x-axis extent.
        ax: Optional matplotlib Axes to draw on. If None, creates a new
            figure and calls `plt.show()`.

    Returns:
        The matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    if sasa_result is None:
        logger.warning("plot_sasa_profile: sasa_result is None.")
        if ax is None:
            _, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE_TRACK)
        ax.text(0.5, 0.5, "No SASA data", ha="center", va="center")
        return ax

    residue_ids = np.asarray(sasa_result.residue_ids)
    values = np.asarray(sasa_result.sasa_per_residue)

    created_fig = ax is None
    if created_fig:
        _, ax = plt.subplots(figsize=_DEFAULT_FIGSIZE_TRACK)

    ax.fill_between(residue_ids, values, color="#27ae60", alpha=0.3)
    ax.plot(residue_ids, values, color="#27ae60", linewidth=1)
    ax.set_xlabel("Residue position (chain numbering)")
    ax.set_ylabel("SASA (\u00c5\u00b2)")
    ax.set_title("Per-residue solvent-accessible surface area")
    if sequence_length:
        ax.set_xlim(0, sequence_length)

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax


def plot_protein_overview(
    sequence_length: Optional[int] = None,
    pdb_records: Optional[List[Dict[str, Any]]] = None,
    mutation_df: Optional[pd.DataFrame] = None,
    sasa_result=None,
    title: Optional[str] = None,
):
    """
    Combine structural coverage, mutation lollipop, and SASA profile into
    one figure with a shared residue-position x-axis. Any subset of the
    three data sources can be provided — only matching panels are drawn.

    Args:
        sequence_length: Full sequence length; recommended so all panels
            share the same x-axis extent.
        pdb_records: See `plot_structural_coverage`.
        mutation_df: See `plot_mutation_lollipop`.
        sasa_result: See `plot_sasa_profile`.
        title: Optional overall figure title (e.g. a UniProt ID/gene name).

    Returns:
        The matplotlib Figure.
    """
    import matplotlib.pyplot as plt

    panels = []
    if pdb_records is not None:
        panels.append(("coverage", pdb_records))
    if mutation_df is not None:
        panels.append(("mutations", mutation_df))
    if sasa_result is not None:
        panels.append(("sasa", sasa_result))

    if not panels:
        raise ValueError("Provide at least one of pdb_records, mutation_df, or sasa_result.")

    fig, axes = plt.subplots(
        nrows=len(panels), ncols=1, sharex=True,
        figsize=(12, 3 * len(panels)),
        gridspec_kw={"hspace": 0.4},
    )
    if len(panels) == 1:
        axes = [axes]

    for (kind, data), ax in zip(panels, axes):
        if kind == "coverage":
            plot_structural_coverage(data, sequence_length=sequence_length, ax=ax)
        elif kind == "mutations":
            plot_mutation_lollipop(data, sequence_length=sequence_length, ax=ax)
        elif kind == "sasa":
            plot_sasa_profile(data, sequence_length=sequence_length, ax=ax)

    if title:
        fig.suptitle(title, fontsize=13, y=1.02)

    return fig
