"""
PDB structural coverage details for a UniProt accession, and generation of
ChimeraX scripts that apply point mutations to a downloaded structure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from ._client import ISDA_BASE_URL, ISDARequestError, get_json, logger


def get_structure_details(
    uniprot_id: str, include_pdb_records: bool = False
) -> Optional[Dict[str, Any]]:
    """
    Fetch PDB structural coverage summary for a UniProt accession.

    Args:
        uniprot_id: UniProt accession, e.g. "P00533".
        include_pdb_records: If True, also include the full list of raw
            per-PDB mapping records under "pdb_records".

    Returns:
        A dict of extracted fields, or None on request failure.
    """
    url = f"{ISDA_BASE_URL}/protein_detail/{uniprot_id}/?additional_outputs=structuralMapping"

    try:
        data = get_json(url)
    except ISDARequestError as exc:
        logger.warning(str(exc))
        return None

    structural_mapping = data.get("structuralMapping", {}) or {}
    highest_coverage = structural_mapping.get("highest_coverage", {}) or {}

    extracted = {
        "pdb_count": structural_mapping.get("pdb_count"),
        "chain_count": structural_mapping.get("chains_count"),
        "highest_coverage_pdb_id": highest_coverage.get("pdb_id", "NA"),
        "highest_coverage": highest_coverage.get("value", "NA"),
        "highest_coverage_pdb_chain": highest_coverage.get("pdb_chain", "NA"),
        "highest_coverage_auth_chain": highest_coverage.get("auth_chain", "NA"),
    }

    if include_pdb_records:
        extracted["pdb_records"] = structural_mapping.get("pdb_data", [])

    return extracted


def get_pdb_summary(
    pdb_id: str,
    additional_outputs: Optional[Union[str, List[str]]] = None,
) -> Optional[Dict[str, Any]]:
    """
    Fetch PDB-entry-level summary information from the ISDA API, e.g.
    ``{ISDA_BASE_URL}/pdb_summary/6q0j/?additional_outputs=micromolecular_data``.

    Args:
        pdb_id: PDB accession, e.g. "6q0j" (case-insensitive; sent lowercase).
        additional_outputs: Optional extra output block(s) to request from
            the API, e.g. "micromolecular_data", or a list of several
            (joined with commas), matching the `additional_outputs` query
            parameter used elsewhere in the ISDA API.

    Returns:
        The raw parsed JSON response as a dict, or None on request
        failure. This endpoint's response schema isn't otherwise
        documented, so inspect the returned dict directly (e.g. via
        `.keys()`) rather than assuming fixed fields.
    """
    pdb_id_lower = pdb_id.lower()
    url = f"{ISDA_BASE_URL}/pdb_summary/{pdb_id_lower}/"

    params = None
    if additional_outputs:
        if isinstance(additional_outputs, (list, tuple)):
            params = {"additional_outputs": ",".join(additional_outputs)}
        else:
            params = {"additional_outputs": additional_outputs}

    try:
        return get_json(url, params=params)
    except ISDARequestError as exc:
        logger.warning(str(exc))
        return None


def generate_mutated_structure_script(
    mutation_table: pd.DataFrame,
    pdb_id: str,
    auth_chain_id: str,
    output_dir: Union[str, Path] = ".",
) -> Path:
    """
    Write a ChimeraX (.cxc) script that opens a structure, applies the given
    mutations via `swapaa`, and saves the mutated structure.

    Args:
        mutation_table: DataFrame with `pdb_auth_chain`, `Alt_AA`
            (three-letter code), and `pdb_residue` columns.
        pdb_id: PDB ID of the structure to mutate (used to build filenames;
            the corresponding `<pdb_id>.cif` is expected in `output_dir`
            when the script is later run in ChimeraX).
        auth_chain_id: Author (deposited) chain ID to restrict mutations to.
        output_dir: Directory to write the script into. Defaults to the
            current working directory.

    Returns:
        Path to the generated .cxc script.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    script_path = output_dir / f"ChimeraXSession_{pdb_id}.cxc"
    input_cif = f"{pdb_id}.cif"
    output_cif = f"{pdb_id}_mutated_structure.cif"

    chain_rows = mutation_table[mutation_table["pdb_auth_chain"] == auth_chain_id]

    lines = [f"open {input_cif}"]
    for _, row in chain_rows.iterrows():
        alt_aa = row["Alt_AA"]
        position = row["pdb_residue"]
        if pd.isna(alt_aa) or alt_aa in ("", "nan"):
            continue
        lines.append(f"swapaa /{auth_chain_id}:{position} {alt_aa}")
    lines.append(f"save {output_cif}")
    lines.append("cli quit")

    script_path.write_text("\n".join(lines) + "\n")
    logger.info("ChimeraX script written to %s", script_path)
    return script_path
