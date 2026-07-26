"""
Per-residue solvent-accessible surface area (SASA) calculation for a PDB
chain, using biotite.

Structure files are fetched from the IBDC ISDA download endpoint by
default (``{ISDA_BASE_URL}/download.<pdb_id>.cif``), which replaces the
previous RCSB-only fetch. RCSB and local files are still supported as
fallbacks.
"""

from __future__ import annotations

import os
from typing import NamedTuple, Optional

import numpy as np

from ._client import ISDA_BASE_URL, ISDARequestError, download_file, logger


class SASAResult(NamedTuple):
    residue_ids: np.ndarray
    residue_names: np.ndarray
    sasa_per_residue: np.ndarray


def fetch_structure_cif(pdb_id: str, output_dir: str = ".") -> str:
    """
    Download a structure's CIF file from the IBDC ISDA download endpoint:
    ``{ISDA_BASE_URL}/download.<pdb_id>.cif``.

    Args:
        pdb_id: 4-character PDB accession, e.g. "6gel" (case-insensitive).
        output_dir: Directory to save the downloaded file into. Defaults
            to the current working directory.

    Returns:
        Path to the downloaded `.cif` file.

    Raises:
        ISDARequestError: on timeout, connection error, HTTP error, or an
            empty response body.
    """
    pdb_id_lower = pdb_id.lower()
    url = f"{ISDA_BASE_URL}/download.{pdb_id_lower}.cif"
    os.makedirs(output_dir, exist_ok=True)
    dest_path = os.path.join(output_dir, f"{pdb_id_lower}.cif")
    return download_file(url, dest_path)


def _read_structure_file(file_path: str, model: int = 1):
    """Read a .cif/.pdb file into a biotite AtomArray, dispatching on extension."""
    ext = os.path.splitext(file_path)[1].lower()

    if ext in (".cif", ".mmcif"):
        from biotite.structure.io.pdbx import CIFFile, get_structure
        cif_file = CIFFile.read(file_path)
        return get_structure(cif_file, model=model)

    from biotite.structure.io import pdb as pdb_io
    pdb_file = pdb_io.PDBFile.read(file_path)
    return pdb_file.get_structure(model=model)


def calculate_sasa_for_chain(
    pdb_id: str,
    chain_id: str,
    probe_radius: float = 1.4,
    local_pdb_path: Optional[str] = None,
    source: str = "isda",
) -> Optional[SASAResult]:
    """
    Calculate per-residue solvent-accessible surface area (SASA) for a
    specific chain of a PDB structure.

    Args:
        pdb_id: 4-character PDB accession, e.g. "6gel".
        chain_id: Chain identifier, e.g. "A".
        probe_radius: Solvent probe radius in Angstrom (default 1.4, water).
        local_pdb_path: Optional path to a local structure file (.cif or
            .pdb). If given, this is used and `source` is ignored.
        source: Where to fetch the structure from when `local_pdb_path`
            is not given: "isda" (default) downloads a `.cif` from the
            IBDC ISDA API; "rcsb" falls back to fetching a `.pdb` file
            from RCSB.

    Returns:
        A `SASAResult(residue_ids, residue_names, sasa_per_residue)`
        namedtuple of aligned NumPy arrays, or None on error.
    """
    # Imported lazily: biotite is an optional/heavier dependency only
    # needed for this module.
    import biotite.structure as struc

    try:
        if local_pdb_path:
            file_path = local_pdb_path
        elif source == "rcsb":
            import biotite.database.rcsb as rcsb
            file_path = rcsb.fetch(pdb_id, "pdb", os.getcwd())
        elif source == "isda":
            file_path = fetch_structure_cif(pdb_id, os.getcwd())
        else:
            raise ValueError(f"Unknown source '{source}'; expected 'isda' or 'rcsb'.")

        array = _read_structure_file(file_path)
    except ISDARequestError as exc:
        logger.warning("Error fetching structure %s from %s: %s", pdb_id, source, exc)
        return None
    except Exception as exc:
        logger.warning("Error reading structure file for %s: %s", pdb_id, exc)
        return None

    if chain_id not in array.chain_id:
        logger.warning("Chain ID '%s' not found in PDB structure %s.", chain_id, pdb_id)
        return None

    chain = array[array.chain_id == chain_id]
    protein_chain = chain[struc.filter_amino_acids(chain)]

    sasa_per_atom = struc.sasa(protein_chain, vdw_radii="ProtOr", probe_radius=probe_radius)
    sasa_per_residue = struc.apply_residue_wise(protein_chain, sasa_per_atom, np.sum)
    residue_ids, residue_names = struc.get_residues(protein_chain)

    return SASAResult(residue_ids, residue_names, sasa_per_residue)
