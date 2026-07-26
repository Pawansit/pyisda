"""
UniProt-numbering <-> PDB-numbering residue mapping, built from the
structural mapping records exposed by the ISDA API.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from ._client import ISDA_BASE_URL, ISDARequestError, get_json, logger


def get_residue_map(
    uniprot_id: str,
    pdb_id: str,
    auth_chain_id: Optional[str] = None,
) -> Optional[List[Dict[str, str]]]:
    """
    Build a residue-by-residue mapping between UniProt numbering and PDB
    numbering for a given (uniprot_id, pdb_id) pair.

    Args:
        uniprot_id: UniProt accession, e.g. "P00533".
        pdb_id: 4-character PDB ID, e.g. "1xyz" (case-insensitive).
        auth_chain_id: Optional author (deposited) chain ID, e.g. "A". If
            given, only mapping records whose `pdb_auth_chain` matches are
            returned — equivalent to but more efficient than fetching the
            full map and then filtering with
            `df[df['pdb_auth_chain'] == chain_id]`, since non-matching
            chains are skipped before the residue map is even built.

    Returns:
        A list of dicts, each with keys `unp_residue`, `pdb_residue`,
        `pdb_auth_chain`, `pdb_chain`; an empty list if the PDB ID (or
        PDB ID + chain) has no matching records; or None if the API
        request itself failed.
    """
    pdb_id_upper = pdb_id.upper()
    url = f"{ISDA_BASE_URL}/protein_detail/{uniprot_id}/?additional_outputs=structuralMapping"

    try:
        data = get_json(url)
    except ISDARequestError as exc:
        logger.warning(str(exc))
        return None

    try:
        pdb_records = data["structuralMapping"]["pdb_data"]
    except KeyError:
        logger.warning("'structuralMapping.pdb_data' missing from API response for %s", uniprot_id)
        return None

    residue_map: List[Dict[str, str]] = []

    for record in pdb_records:
        if record.get("PDB") != pdb_id_upper:
            continue

        if auth_chain_id is not None and str(record.get("AUTH_CHAIN", "")) != str(auth_chain_id):
            continue

        unp_start = int(record["UNIPROT_START"])
        unp_end = int(record["UNIPROT_END"])

        try:
            pdb_start = int(record["PDB_START"])
        except (ValueError, TypeError):
            logger.warning(
                "Non-integer PDB_START (%s) for %s/%s; skipping this range",
                record.get("PDB_START"), uniprot_id, pdb_id_upper,
            )
            continue

        unp_range = range(unp_start, unp_end + 1)
        pdb_range = range(pdb_start, pdb_start + len(unp_range))

        for unp_residue, pdb_residue in zip(unp_range, pdb_range):
            residue_map.append({
                "unp_residue": str(unp_residue),
                "pdb_residue": str(pdb_residue),
                "pdb_auth_chain": str(record.get("AUTH_CHAIN", "N/A")),
                "pdb_chain": str(record.get("PDB_CHAIN", "N/A")),
            })

    if not residue_map and pdb_records:
        logger.info(
            "No mapping found for PDB ID %s%s within available records for %s",
            pdb_id_upper,
            f" chain {auth_chain_id}" if auth_chain_id else "",
            uniprot_id,
        )

    return residue_map
