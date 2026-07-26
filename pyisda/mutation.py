"""
Mutation retrieval (clinical + AlphaMissense-style computational
predictions), property analysis, and mutation-matrix plotting.

Note: the original codebase had two divergent copies of the computational
mutation fetcher (`Mutation.py::AFMutationTable` and a near-duplicate in
`Afmissense.py` with a malformed URL). This module keeps a single,
corrected implementation: `get_computational_mutations`.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ._client import ISDA_BASE_URL, ISDARequestError, get_json, logger

_CLINICAL_AA_PATTERN = r"p\.(?P<Ref_AA>[A-Z][a-z]{2})(?P<Position>\d+)(?P<Alt_AA>[A-Z][a-z]{2}|\*)"
_COMPUTATIONAL_AA_PATTERN = r"^([A-Z])(\d+)([A-Z])$"

DEFAULT_HYDROPHOBIC_AAS = ["Ala", "Val", "Leu", "Ile", "Pro", "Phe", "Trp", "Met", "Cys"]
DEFAULT_HYDROPHILIC_AAS = ["Gly", "Ser", "Thr", "Tyr", "Asn", "Gln", "Asp", "Glu", "Lys", "Arg", "His", "*"]


def get_mutation_table(
    uniprot_id: str,
    source: str = "clinvar",
    record_type: str = "default",
    selection: Optional[list] = None,
) -> Optional[pd.DataFrame]:
    """
    Fetch clinical mutation records (e.g. from ClinVar) for a UniProt
    accession and parse the protein-level AA change into Ref_AA / Position
    / Alt_AA columns.

    Args:
        uniprot_id: UniProt accession.
        source: Mutation data source, e.g. "clinvar".
        record_type: Record type/output shape requested from the API.
        selection: Optional list of columns to restrict the result to.

    Returns:
        A DataFrame of mutation records, or None on failure.
    """
    url = f"{ISDA_BASE_URL}/mutation/?uniprot_id={uniprot_id}&sources={source}&output_type={record_type}"

    try:
        data = get_json(url)
    except ISDARequestError as exc:
        logger.warning(str(exc))
        return None

    if "mutation" not in data or not data["mutation"]:
        logger.info("No mutation records found for %s", uniprot_id)
        return None

    df = pd.DataFrame(data["mutation"])
    df[["Ref_AA", "Position", "Alt_AA"]] = df["Variant AA Change"].str.extract(_CLINICAL_AA_PATTERN)

    if selection is None:
        return df

    missing = [col for col in selection if col not in df.columns]
    if missing:
        logger.warning("Missing columns in selection: %s", missing)
        return None
    return df[selection]


def get_computational_mutations(
    uniprot_id: str,
    ranges: Optional[str] = None,
    significance_type: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """
    Fetch computational (e.g. AlphaMissense) mutation predictions for a
    UniProt accession and parse Ref_AA / Position / Alt_AA columns.

    Args:
        uniprot_id: UniProt accession.
        ranges: Optional residue range filter, API-specific format.
        significance_type: Optional significance-category filter.

    Returns:
        A DataFrame of predicted mutations, or None on failure.
    """
    url = f"{ISDA_BASE_URL}/computational_mutations/{uniprot_id}"

    if ranges and significance_type:
        url += f"/?ranges={ranges}&significance_type={significance_type}"
    elif ranges:
        url += f"/?ranges={ranges}"
    elif significance_type:
        url += f"?significance_type={significance_type}"

    try:
        data = get_json(url)
    except ISDARequestError as exc:
        logger.warning(str(exc))
        return None

    if data is None:
        logger.info("No computational mutation data returned for %s", uniprot_id)
        return None

    df = pd.DataFrame(data)
    df[["Ref_AA", "Position", "Alt_AA"]] = df["protein_variant"].str.extract(_COMPUTATIONAL_AA_PATTERN)
    return df


def filter_by_significance(df: pd.DataFrame, pattern: str, case: bool = False) -> pd.DataFrame:
    """
    Filter a mutation DataFrame by a regex pattern against the
    `Clinical Significance` column.

    Args:
        df: DataFrame containing a `Clinical Significance` column.
        pattern: Regex pattern to search for.
        case: Case-sensitive match if True (default: False).

    Returns:
        A filtered copy of `df` with a reset index.
    """
    if "Clinical Significance" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'Clinical Significance' column.")
    if not isinstance(pattern, str) or not pattern.strip():
        raise ValueError("A valid search pattern (string) must be provided.")

    return df[
        df["Clinical Significance"].str.contains(pattern, case=case, na=False)
    ].reset_index(drop=True)


def _classify_change(row: pd.Series) -> str:
    ref, alt = row["Ref_Prop"], row["Alt_Prop"]
    if pd.isna(ref) or pd.isna(alt):
        return "Unknown/N/A"
    return f"No change ({ref})" if ref == alt else f"Change from {ref} to {alt}"


def analyze_mutation_properties(
    df: pd.DataFrame,
    hydrophobic_aas: list = DEFAULT_HYDROPHOBIC_AAS,
    hydrophilic_aas: list = DEFAULT_HYDROPHILIC_AAS,
    category_labels: list = ("Hydrophobic", "Hydrophilic"),
) -> Optional[pd.DataFrame]:
    """
    Annotate a mutation DataFrame with hydrophobic/hydrophilic property
    classes for Ref_AA and Alt_AA, plus a combined change-type label.

    Returns:
        `df` with `Ref_Prop`, `Alt_Prop`, `Mutation_Effect_Type` columns
        added, or None if `df` is empty/None.
    """
    if df is None or df.empty:
        logger.warning("Cannot analyze properties: DataFrame is empty or None.")
        return None

    prop_map = {aa: category_labels[0] for aa in hydrophobic_aas}
    prop_map.update({aa: category_labels[1] for aa in hydrophilic_aas})

    df = df.copy()
    df["Ref_Prop"] = df["Ref_AA"].map(prop_map).fillna("Unknown")
    df["Alt_Prop"] = df["Alt_AA"].map(prop_map).fillna("Unknown")
    df["Mutation_Effect_Type"] = df.apply(_classify_change, axis=1)
    return df


def get_mutation_summary(df: pd.DataFrame, top_n: int = 10) -> Optional[Dict]:
    """
    Summarize a mutation DataFrame: total count, clinical significance
    distribution, consequence-type distribution, and top mutated positions.

    Args:
        df: DataFrame with `Clinical Significance` and `Position` columns.
        top_n: Number of top mutation-hotspot positions to report.

    Returns:
        A summary dict, or None if `df` is invalid/missing required columns.
    """
    if df is None or df.empty:
        logger.warning("Cannot summarize: DataFrame is empty or None.")
        return None

    required = ["Clinical Significance", "Position"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        logger.warning("DataFrame is missing required columns: %s", missing)
        return None

    return {
        "mutation_count": df.shape[0],
        "significance_status": df["Clinical Significance"].value_counts().to_dict(),
        "mutation_type": df["Consequence Type"].value_counts(sort=True).to_dict(),
        "maximum_mutation_spot": df["Position"].value_counts(sort=True).head(top_n).to_dict(),
    }


def merge_mutations_with_pdb_map(
    mutation_df: pd.DataFrame, residue_map_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Join a mutation table to a UniProt<->PDB residue map on residue
    position, so each mutation gets its corresponding PDB-numbered residue.

    Args:
        mutation_df: Output of `get_mutation_table` / `get_computational_mutations`
            (must contain a `Position` column).
        residue_map_df: Output of `mapper.get_residue_map` as a DataFrame,
            whose first column will be renamed to `Position` for the join.
            Works whether `residue_map_df` is the full mapping or a
            filtered/sliced subset of it (e.g. one chain).

    Returns:
        The merged DataFrame.
    """
    residue_map_df = residue_map_df.rename(columns={residue_map_df.columns[0]: "Position"})
    return pd.merge(mutation_df, residue_map_df, on="Position")


def _normalize_alt_aa(series: pd.Series) -> pd.Series:
    """Map a mixed 1-/3-letter Alt_AA column to single-letter codes for joining."""
    from .protein import AA_3TO1

    def _convert(value):
        if pd.isna(value):
            return value
        value = str(value)
        if len(value) == 1:
            return value.upper()
        return AA_3TO1.get(value, value)

    return series.apply(_convert)


def merge_experimental_with_computational(
    experimental_df: pd.DataFrame,
    computational_df: pd.DataFrame,
    position_col: str = "Position",
    alt_aa_col: str = "Alt_AA",
    how: str = "left",
    suffixes: Tuple[str, str] = ("_experimental", "_computational"),
) -> pd.DataFrame:
    """
    Attach computational (e.g. AlphaMissense) pathogenicity predictions to
    each experimentally-observed mutation, so both sources are visible
    for the same substitution side by side.

    Joins on residue position + substituted amino acid. Since
    `get_mutation_table` reports `Alt_AA` as a 3-letter code (e.g. "Gly")
    and `get_computational_mutations` reports it as 1-letter (e.g. "G"),
    both are normalized to 1-letter internally before joining — the
    original `Position`/`Alt_AA` values from each side are dropped in
    favor of one canonical `Position`/`Alt_AA` pair in the result (they
    always agree by construction of the join), while every other
    overlapping column (e.g. `Ref_AA`) is kept from both sides, suffixed
    per `suffixes`.

    Args:
        experimental_df: Output of `get_mutation_table` (or any DataFrame
            with `position_col`/`alt_aa_col` columns).
        computational_df: Output of `get_computational_mutations` (or any
            DataFrame with `position_col`/`alt_aa_col` columns).
        position_col: Residue-position column name, present in both inputs.
        alt_aa_col: Substituted-amino-acid column name, present in both
            inputs (mixed 1-/3-letter codes are handled automatically).
        how: Join type, passed to `pandas.merge` — "left" (default) keeps
            every experimental mutation and attaches a computational
            prediction where one exists (NaN otherwise); use "inner" to
            keep only mutations with both experimental and computational
            support.
        suffixes: Suffixes applied to other overlapping columns from the
            two inputs.

    Returns:
        The merged DataFrame with canonical `Position`/`Alt_AA` columns.
    """
    if experimental_df is None or experimental_df.empty:
        logger.warning("merge_experimental_with_computational: experimental_df is empty or None.")
        return experimental_df
    if computational_df is None or computational_df.empty:
        logger.warning("merge_experimental_with_computational: computational_df is empty or None.")
        return experimental_df

    for df_, name in ((experimental_df, "experimental_df"), (computational_df, "computational_df")):
        missing = [c for c in (position_col, alt_aa_col) if c not in df_.columns]
        if missing:
            raise ValueError(f"{name} is missing required column(s): {missing}")

    left = experimental_df.copy()
    right = computational_df.copy()

    left["_join_position"] = left[position_col].astype(str)
    right["_join_position"] = right[position_col].astype(str)
    left["_join_alt_aa"] = _normalize_alt_aa(left[alt_aa_col])
    right["_join_alt_aa"] = _normalize_alt_aa(right[alt_aa_col])

    merged = pd.merge(left, right, on=["_join_position", "_join_alt_aa"], how=how, suffixes=suffixes)

    # Drop the original (now-suffixed, redundant-by-construction) position/
    # Alt_AA columns from both sides in favor of the canonical join columns.
    redundant = [f"{position_col}{suffixes[0]}", f"{position_col}{suffixes[1]}",
                 f"{alt_aa_col}{suffixes[0]}", f"{alt_aa_col}{suffixes[1]}"]
    merged = merged.drop(columns=[c for c in redundant if c in merged.columns])
    merged = merged.rename(columns={"_join_position": position_col, "_join_alt_aa": alt_aa_col})

    return merged


def subset_by_positions(
    df: pd.DataFrame,
    positions: Sequence[Union[int, str]],
    position_col: str = "Position",
) -> pd.DataFrame:
    """
    Filter a mutation (or merged mutation) DataFrame down to a specific
    set of residue positions — e.g. a ligand-binding-site or active-site
    residue list — to focus downstream analysis on that region.

    Args:
        df: Any DataFrame with a residue-position column, such as
            `get_mutation_table`, `get_computational_mutations`, or
            `merge_experimental_with_computational` output.
        positions: Residue positions to keep (ints or strings; compared
            as strings so "12" and 12 both match).
        position_col: Column to filter on. Note that after
            `merge_mutations_with_pdb_map` this is still `"Position"`,
            but you may want `"pdb_residue"` instead if filtering by
            PDB-numbered binding-site residues.

    Returns:
        A filtered copy of `df` with a reset index.
    """
    if position_col not in df.columns:
        raise ValueError(f"'{position_col}' column not found in DataFrame.")

    wanted = {str(p) for p in positions}
    return df[df[position_col].astype(str).isin(wanted)].reset_index(drop=True)


def plot_mutation_matrix(df: pd.DataFrame, ax=None):
    """
    Plot a heatmap of mutation counts per (Alt_AA, Position).

    Args:
        df: DataFrame with `Position`, `Alt_AA`, and `AC` columns.
        ax: Optional matplotlib Axes to draw on. If None, creates a new
            figure and calls `plt.show()`.
    """
    import matplotlib.pyplot as plt

    if df is None or df.empty:
        logger.warning("Cannot plot; DataFrame is empty or None.")
        return

    matrix_df = pd.crosstab(index=df["Alt_AA"], columns=df["Position"])

    created_fig = ax is None
    if created_fig:
        fig, ax = plt.subplots(figsize=(15, 8))

    im = ax.imshow(matrix_df.values, cmap="viridis", aspect="auto")
    ax.set_xticks(np.arange(len(matrix_df.columns)))
    ax.set_xticklabels(matrix_df.columns, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(matrix_df.index)))
    ax.set_yticklabels(matrix_df.index, fontsize=10)
    ax.set_xlabel("Mutation Position (AA)", fontsize=12)
    ax.set_ylabel("Alternative Amino Acid (Alt_AA)", fontsize=12)
    ax.set_title(f"Mutation Matrix for UniProt ID: {df['AC'].iloc[0]}")
    plt.colorbar(im, ax=ax, label="Mutation Count")

    if created_fig:
        plt.tight_layout()
        plt.show()
