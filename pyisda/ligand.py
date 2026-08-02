"""
Bound-ligand lookup and structural binding-site residue detection.

`get_bound_ligands` uses the ISDA `pdb_summary` endpoint (no biotite
required). `get_binding_site_residues` needs an actual 3D structure, so
it reuses `sasa.fetch_structure_array` (biotite) to fetch/parse one and
computes residues within a distance cutoff of the ligand.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

import numpy as np

from ._client import ISDARequestError, logger

# Candidate key names to look for the ligand/micromolecule record list in
# the pdb_summary response, and candidate chain-identifying fields within
# each record. The exact schema of `additional_outputs=micromolecular_data`
# isn't otherwise documented, so this is deliberately permissive.
_LIGAND_LIST_KEYS = (
    "micromolecular_data", "micromolecules", "ligands",
    "hetero_compounds", "non_polymer_entities",
)
_CHAIN_KEYS = (
    "auth_asym_id", "chain_id", "auth_chain", "pdb_chain",
    "chain", "AUTH_CHAIN", "PDB_CHAIN",
)


def get_bound_ligands(pdb_id: str, chain_id: Optional[str] = None) -> Optional[Any]:
    """
    Fetch bound small-molecule (ligand/HET) information for a PDB entry,
    via the ISDA `pdb_summary` endpoint's `micromolecular_data` output.

    Args:
        pdb_id: PDB accession, e.g. "6q0j".
        chain_id: Optional chain ID (`auth_asym_id`) to filter results to.
            If omitted, all ligand records for the entry are returned.

    Returns:
        A list of ligand-record dicts, or None on request failure. If the
        response's ligand list uses a recognizable chain field (checked
        among `auth_asym_id`, `chain_id`, `auth_chain`, `pdb_chain`,
        `chain`), results are filtered to `chain_id`. If no ligand list or
        no recognizable chain field can be found — this endpoint's exact
        schema isn't otherwise documented — the raw `pdb_summary` response
        (or its unfiltered ligand list) is returned instead, with a
        logged warning, so inspect its keys directly in that case.
    """
    from .structure import get_pdb_summary  # local import: avoid a circular import at module load

    data = get_pdb_summary(pdb_id, additional_outputs="micromolecular_data")
    if data is None:
        return None

    ligand_records = None
    for key in _LIGAND_LIST_KEYS:
        value = data.get(key) if isinstance(data, dict) else None
        if isinstance(value, list):
            ligand_records = value
            break

    if ligand_records is None:
        logger.warning(
            "Could not locate a ligand/micromolecule list in the pdb_summary "
            "response for %s (checked keys: %s). Returning the raw response "
            "instead — inspect its keys directly.",
            pdb_id, ", ".join(_LIGAND_LIST_KEYS),
        )
        return data

    if not chain_id:
        return ligand_records

    filterable = [r for r in ligand_records if isinstance(r, dict) and any(k in r for k in _CHAIN_KEYS)]
    if not filterable:
        logger.warning(
            "Ligand records for %s don't have a recognizable chain field "
            "(checked: %s); returning all %d record(s) unfiltered.",
            pdb_id, ", ".join(_CHAIN_KEYS), len(ligand_records),
        )
        return ligand_records

    def _matches_chain(record: Dict[str, Any]) -> bool:
        return any(
            key in record and str(record[key]).upper() == chain_id.upper()
            for key in _CHAIN_KEYS
        )

    return [r for r in ligand_records if _matches_chain(r)]


def get_binding_site_residues(
    pdb_id: str,
    chain_id: str,
    ligand_id: str,
    cutoff: float = 5.0,
    local_pdb_path: Optional[str] = None,
    source: str = "isda",
) -> Optional[List[Dict[str, Any]]]:
    """
    Identify the residues surrounding a bound ligand within a distance
    cutoff — i.e. the ligand's binding-site/pocket residues.

    Args:
        pdb_id: PDB accession, e.g. "6q0j".
        chain_id: The ligand instance's chain ID (`auth_asym_id` in the
            structure file). If unsure, check with `get_bound_ligands`
            first for the ligand's chain.
        ligand_id: 3-letter ligand/HET code, e.g. "TPO" (case-insensitive)
            — matched against each residue's name.
        cutoff: Distance cutoff in Angstrom for "surrounding" (default
            5.0, a common binding-site definition).
        local_pdb_path: Optional local structure file (.cif or .pdb) to
            use instead of fetching one.
        source: Where to fetch the structure from when `local_pdb_path`
            is not given: "isda" (default, IBDC ISDA download endpoint)
            or "rcsb". See `sasa.calculate_sasa_for_chain`.

    Returns:
        A list of dicts, one per binding-site residue:
        `{chain_id, residue_id, residue_name, min_distance}`, sorted by
        chain then residue ID. Returns an empty list if the ligand isn't
        found in the given chain. Returns None on a structure fetch/read
        error.
    """
    # Imported lazily: biotite is an optional/heavier dependency only
    # needed for structure-consuming functions.
    import biotite.structure as struc
    from .sasa import fetch_structure_array

    try:
        array = fetch_structure_array(pdb_id, local_pdb_path=local_pdb_path, source=source)
    except ISDARequestError as exc:
        logger.warning("Error fetching structure %s from %s: %s", pdb_id, source, exc)
        return None
    except Exception as exc:
        logger.warning("Error reading structure file for %s: %s", pdb_id, exc)
        return None

    ligand_mask = (array.chain_id == chain_id) & (array.res_name == ligand_id.upper())
    ligand_atoms = array[ligand_mask]

    if ligand_atoms.array_length() == 0:
        logger.warning(
            "Ligand '%s' not found in chain '%s' of structure %s.",
            ligand_id, chain_id, pdb_id,
        )
        return []

    protein_atoms = array[struc.filter_amino_acids(array)]
    if protein_atoms.array_length() == 0:
        logger.warning("No amino-acid residues found in structure %s.", pdb_id)
        return []

    # Minimum distance from each protein atom to any ligand atom.
    diff = protein_atoms.coord[:, None, :] - ligand_atoms.coord[None, :, :]
    min_dist_per_atom = np.linalg.norm(diff, axis=-1).min(axis=1)

    within_cutoff = min_dist_per_atom <= cutoff
    nearby_atoms = protein_atoms[within_cutoff]
    nearby_dist = min_dist_per_atom[within_cutoff]

    residues: Dict[tuple, Dict[str, Any]] = {}
    for c, res_id, res_name, dist in zip(
        nearby_atoms.chain_id, nearby_atoms.res_id, nearby_atoms.res_name, nearby_dist
    ):
        key = (str(c), int(res_id))
        if key not in residues or dist < residues[key]["min_distance"]:
            residues[key] = {
                "chain_id": str(c),
                "residue_id": int(res_id),
                "residue_name": str(res_name),
                "min_distance": round(float(dist), 2),
            }

    return sorted(residues.values(), key=lambda r: (r["chain_id"], r["residue_id"]))
