"""
Mutation retrieval (clinical + AlphaMissense-style computational
predictions), property analysis, and mutation-matrix plotting.

Note: the original codebase had two divergent copies of the computational
mutation fetcher (`Mutation.py::AFMutationTable` and a near-duplicate in
`Afmissense.py` with a malformed URL). This module keeps a single,
corrected implementation: `get_computational_mutations`.
"""

from __future__ import annotations

from typing import Dict, Optional

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
