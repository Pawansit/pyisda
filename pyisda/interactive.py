"""
Interactive (Plotly) counterparts to `pyisda.visualize` — same data
inputs, but rendered as zoomable/pannable HTML with hover tooltips instead
of static matplotlib images. Aimed at cases where a protein has many
structures/mutations and a static plot's y-axis labels or coverage ranges
become too cramped to read.

Requires the `interactive` extra: ``pip install -e ".[interactive]"`` (installs
`plotly`).

Every `plot_*_html` function returns a `plotly.graph_objects.Figure` that
you can display in a notebook (`fig.show()`) or save with `save_html`.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from ._client import logger


def _require_plotly():
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        return go, make_subplots
    except ImportError as exc:
        raise ImportError(
            "Interactive HTML plots require plotly. Install it with: "
            'pip install -e ".[interactive]"  (or: pip install plotly)'
        ) from exc


def save_html(fig, path: str, standalone: bool = True, auto_open: bool = False) -> str:
    """
    Save a Plotly figure to a standalone HTML file.

    Args:
        fig: A `plotly.graph_objects.Figure`, e.g. from one of the
            `plot_*_html` functions in this module.
        path: Output file path, e.g. "coverage.html".
        standalone: If True (default), embeds the Plotly JS library
            directly in the file (~3-4 MB, works fully offline). If
            False, loads Plotly JS from a CDN instead (much smaller file,
            requires internet access to view).
        auto_open: If True, opens the file in a browser after saving.

    Returns:
        The path the file was written to.
    """
    fig.write_html(path, include_plotlyjs=(True if standalone else "cdn"), auto_open=auto_open)
    logger.info("Saved interactive HTML to %s", path)
    return path


def plot_structural_coverage_html(
    pdb_records: List[Dict[str, Any]],
    sequence_length: Optional[int] = None,
    title: Optional[str] = None,
):
    """
    Interactive structural-coverage track: one horizontal bar per (PDB,
    chain) record. Hover shows the exact PDB ID, chain, and UniProt
    range; drag to zoom/pan on either axis (much more legible than a
    static plot when there are dozens of structures).

    Args:
        pdb_records: `pdb_records` from
            `structure.get_structure_details(uniprot_id, include_pdb_records=True)`.
        sequence_length: Full sequence length, to fix the x-axis extent.
        title: Optional figure title.

    Returns:
        A `plotly.graph_objects.Figure`.
    """
    go, _ = _require_plotly()

    if not pdb_records:
        fig = go.Figure()
        fig.add_annotation(text="No structural coverage records", showarrow=False)
        return fig

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
            "span": max(end - start, 1),
        })

    if not rows:
        fig = go.Figure()
        fig.add_annotation(text="No usable coverage records", showarrow=False)
        return fig

    df = pd.DataFrame(rows).sort_values("start").reset_index(drop=True)
    df["label"] = df["pdb_id"] + ":" + df["chain"]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=df["span"],
        y=df["label"],
        base=df["start"],
        orientation="h",
        marker=dict(color=df.index, colorscale="Turbo", showscale=False),
        customdata=np.stack([df["pdb_id"], df["chain"], df["start"], df["end"]], axis=-1),
        hovertemplate=(
            "PDB %{customdata[0]} chain %{customdata[1]}<br>"
            "UniProt %{customdata[2]}\u2013%{customdata[3]}"
            "<extra></extra>"
        ),
    ))

    fig.update_layout(
        title=title or f"Structural coverage ({len(df)} records)",
        xaxis_title="UniProt residue position",
        yaxis=dict(automargin=True, fixedrange=False),
        xaxis=dict(
            range=[0, sequence_length] if sequence_length else None,
            rangeslider=dict(visible=True, thickness=0.08),
        ),
        height=max(350, 24 * len(df) + 150),
        margin=dict(l=120),
        bargap=0.25,
    )
    return fig


def plot_mutation_lollipop_html(
    mutation_df: pd.DataFrame,
    position_col: str = "Position",
    significance_col: Optional[str] = "Clinical Significance",
    pathogenic_pattern: str = "pathogenic",
    sequence_length: Optional[int] = None,
    title: Optional[str] = None,
):
    """
    Interactive mutation lollipop: hover a point to see its exact
    position, mutation count, and (if available) the breakdown of
    clinical significance labels at that position.

    Args:
        mutation_df: Output of `mutation.get_mutation_table` /
            `mutation.get_computational_mutations`.
        position_col: Column holding the (1-indexed) residue position.
        significance_col: Column to check for a pathogenic flag, or None
            to skip coloring.
        pathogenic_pattern: Case-insensitive substring marking a record
            as pathogenic.
        sequence_length: Full sequence length, to fix the x-axis extent.
        title: Optional figure title.

    Returns:
        A `plotly.graph_objects.Figure`.
    """
    go, _ = _require_plotly()

    if mutation_df is None or mutation_df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No mutation records", showarrow=False)
        return fig

    df = mutation_df.copy()
    df[position_col] = pd.to_numeric(df[position_col], errors="coerce")
    df = df.dropna(subset=[position_col])

    has_sig = significance_col and significance_col in df.columns
    if has_sig:
        grouped = df.groupby(position_col)[significance_col].apply(
            lambda vals: ", ".join(sorted(set(vals.dropna())))
        )
        is_pathogenic = df.groupby(position_col)[significance_col].apply(
            lambda vals: vals.str.contains(pathogenic_pattern, case=False, na=False).any()
        )
        counts = df.groupby(position_col).size()
        colors = np.where(is_pathogenic.reindex(counts.index, fill_value=False), "#c0392b", "#2980b9")
        hover_extra = grouped.reindex(counts.index).fillna("")
    else:
        counts = df.groupby(position_col).size()
        colors = "#2980b9"
        hover_extra = None

    fig = go.Figure()
    # stems
    for pos, count in counts.items():
        fig.add_trace(go.Scatter(
            x=[pos, pos], y=[0, count], mode="lines",
            line=dict(color="#bdc3c7", width=1), showlegend=False, hoverinfo="skip",
        ))

    customdata = None
    hovertemplate = "Position %{x}<br>Mutations: %{y}<extra></extra>"
    if hover_extra is not None:
        customdata = hover_extra.values
        hovertemplate = "Position %{x}<br>Mutations: %{y}<br>%{customdata}<extra></extra>"

    fig.add_trace(go.Scatter(
        x=counts.index, y=counts.values, mode="markers",
        marker=dict(size=10, color=colors, line=dict(width=1, color="white")),
        customdata=customdata,
        hovertemplate=hovertemplate,
        showlegend=False,
    ))

    if has_sig:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                  marker=dict(size=10, color="#c0392b"), name="Pathogenic"))
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                  marker=dict(size=10, color="#2980b9"), name="Other"))

    fig.update_layout(
        title=title or f"Mutation positions ({int(counts.sum())} total)",
        xaxis_title="Residue position",
        yaxis_title="Mutation count",
        xaxis=dict(
            range=[0, sequence_length] if sequence_length else None,
            rangeslider=dict(visible=True, thickness=0.08),
        ),
        height=450,
    )
    return fig


def plot_sasa_profile_html(
    sasa_result,
    sequence_length: Optional[int] = None,
    title: Optional[str] = None,
):
    """
    Interactive per-residue SASA profile. Hover shows the exact residue
    id, name, and SASA value.

    Args:
        sasa_result: A `sasa.SASAResult` from `sasa.calculate_sasa_for_chain`.
        sequence_length: Full sequence length, to fix the x-axis extent.
        title: Optional figure title.

    Returns:
        A `plotly.graph_objects.Figure`.
    """
    go, _ = _require_plotly()

    if sasa_result is None:
        fig = go.Figure()
        fig.add_annotation(text="No SASA data", showarrow=False)
        return fig

    residue_ids = np.asarray(sasa_result.residue_ids)
    values = np.asarray(sasa_result.sasa_per_residue)
    names = np.asarray(sasa_result.residue_names) if getattr(sasa_result, "residue_names", None) is not None else None

    customdata = names.reshape(-1, 1) if names is not None else None
    hovertemplate = "Residue %{x}<br>SASA: %{y:.1f} \u00c5\u00b2"
    if customdata is not None:
        hovertemplate = "Residue %{x} (%{customdata[0]})<br>SASA: %{y:.1f} \u00c5\u00b2"
    hovertemplate += "<extra></extra>"

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=residue_ids, y=values, mode="lines",
        line=dict(color="#27ae60", width=1.5),
        fill="tozeroy", fillcolor="rgba(39, 174, 96, 0.2)",
        customdata=customdata,
        hovertemplate=hovertemplate,
    ))

    fig.update_layout(
        title=title or "Per-residue solvent-accessible surface area",
        xaxis_title="Residue position (chain numbering)",
        yaxis_title="SASA (\u00c5\u00b2)",
        xaxis=dict(
            range=[0, sequence_length] if sequence_length else None,
            rangeslider=dict(visible=True, thickness=0.08),
        ),
        height=400,
    )
    return fig


def plot_protein_overview_html(
    sequence_length: Optional[int] = None,
    pdb_records: Optional[List[Dict[str, Any]]] = None,
    mutation_df: Optional[pd.DataFrame] = None,
    sasa_result=None,
    title: Optional[str] = None,
):
    """
    Combine structural coverage, mutation lollipop, and SASA profile into
    one interactive figure with a shared, synced residue-position x-axis
    (zooming/panning one panel moves the others together). Any subset of
    the three data sources can be provided.

    Args:
        sequence_length: Full sequence length; recommended so all panels
            share the same x-axis extent.
        pdb_records: See `plot_structural_coverage_html`.
        mutation_df: See `plot_mutation_lollipop_html`.
        sasa_result: See `plot_sasa_profile_html`.
        title: Optional overall figure title.

    Returns:
        A `plotly.graph_objects.Figure`.
    """
    go, make_subplots = _require_plotly()

    panel_specs = []
    if pdb_records is not None:
        panel_specs.append(("coverage", "Structural coverage"))
    if mutation_df is not None:
        panel_specs.append(("mutations", "Mutations"))
    if sasa_result is not None:
        panel_specs.append(("sasa", "SASA"))

    if not panel_specs:
        raise ValueError("Provide at least one of pdb_records, mutation_df, or sasa_result.")

    n = len(panel_specs)
    row_heights = [0.5 if kind == "coverage" else 0.25 for kind, _ in panel_specs]
    fig = make_subplots(
        rows=n, cols=1, shared_xaxes=True, vertical_spacing=0.08,
        subplot_titles=[label for _, label in panel_specs],
        row_heights=row_heights,
    )

    for i, (kind, _) in enumerate(panel_specs, start=1):
        if kind == "coverage":
            sub = plot_structural_coverage_html(pdb_records, sequence_length=sequence_length)
        elif kind == "mutations":
            sub = plot_mutation_lollipop_html(mutation_df, sequence_length=sequence_length)
        else:
            sub = plot_sasa_profile_html(sasa_result, sequence_length=sequence_length)
        for trace in sub.data:
            fig.add_trace(trace, row=i, col=1)
        # carry over the y-axis title/config for this row
        fig.update_yaxes(title_text=sub.layout.yaxis.title.text, row=i, col=1)

    fig.update_xaxes(
        title_text="Residue position",
        range=[0, sequence_length] if sequence_length else None,
        rangeslider=dict(visible=(n == 1), thickness=0.06),
        row=n, col=1,
    )
    fig.update_layout(
        title=title,
        height=300 * n,
        showlegend=False,
    )
    return fig
