"""
Protein-level metadata: gene names, organism info, sequence retrieval,
and sequence/single-letter <-> three-letter mutation helpers.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import pandas as pd

from ._client import ISDA_BASE_URL, ISDARequestError, get_json, logger

AA_1TO3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "*",
}
AA_3TO1 = {v: k for k, v in AA_1TO3.items()}


def _extract_gene_names(data: Dict[str, Any]) -> List[str]:
    """Extract all `geneName` values from the `genes` list of a protein record."""
    names: List[str] = []
    for gene_entry in data.get("genes", []) or []:
        try:
            value = gene_entry.get("geneName", {}).get("value")
        except AttributeError:
            continue
        if value:
            names.append(value)
    return names


def _extract_gene_synonyms(data: Dict[str, Any]) -> List[str]:
    """Extract all gene synonym values from the `genes` list of a protein record."""
    synonyms: List[str] = []
    for gene_entry in data.get("genes", []) or []:
        try:
            for synonym in gene_entry.get("synonyms", []) or []:
                value = synonym.get("value")
                if value:
                    synonyms.append(value)
        except AttributeError:
            continue
    return synonyms


def get_protein_details(uniprot_id: str) -> Optional[Dict[str, Any]]:
    """
    Fetch summary protein details (organism, gene names, structural coverage)
    for a UniProt accession from the ISDA API.

    Args:
        uniprot_id: UniProt accession, e.g. "P00533".

    Returns:
        A dict of extracted fields, or None on request failure.
    """
    url = f"{ISDA_BASE_URL}/protein_detail/{uniprot_id}/?additional_outputs=proteinFeatures,structuralMapping"

    try:
        data = get_json(url)
    except ISDARequestError as exc:
        logger.warning(str(exc))
        return None

    recommended_name = (
        data.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value")
    )

    return {
        "uniprot_id": data.get("uniprot_id"),
        "scientific_name": data.get("organism", {}).get("scientificName"),
        "taxon_id": data.get("organism", {}).get("taxonId"),
        "protein_name": recommended_name,
        "gene_name": ",".join(_extract_gene_names(data)),
        "gene_synonyms": ",".join(_extract_gene_synonyms(data)),
        "structural_records": data.get("structuralMapping", {}).get("pdb_count", 0),
        "sequence_length": data.get("sequence", {}).get("length"),
    }


def get_sequence(uniprot_id: str) -> Optional[Dict[str, Any]]:
    """
    Fetch the amino-acid sequence, length, and molecular weight for a
    UniProt accession.

    Args:
        uniprot_id: UniProt accession, e.g. "P00533".

    Returns:
        A dict with `sequence`, `length`, `mol_weight`, or None on failure.
    """
    url = f"{ISDA_BASE_URL}/protein_detail/{uniprot_id}/?additional_outputs=proteinFeatures"

    try:
        data = get_json(url)
    except ISDARequestError as exc:
        logger.warning(str(exc))
        return None

    sequence_block = data.get("sequence") or {}
    return {
        "sequence": sequence_block.get("value"),
        "length": sequence_block.get("length"),
        "mol_weight": sequence_block.get("molWeight"),
    }


def mutate_sequence(sequence_record: Dict[str, Any], mutations: pd.DataFrame) -> str:
    """
    Apply a table of point mutations to a sequence and return the mutated
    sequence string (does not modify `sequence_record`).

    Args:
        sequence_record: Output of `get_sequence`, or any dict with a
            "sequence" key holding the wild-type sequence string.
        mutations: DataFrame with `Alt_AA` (three-letter code, e.g. "Ala")
            and `Position` (1-indexed residue number) columns.

    Returns:
        The mutated sequence as a string.
    """
    residues = list(sequence_record.get("sequence", ""))
    if not residues:
        raise ValueError("sequence_record has no 'sequence' value to mutate")

    records = (
        mutations[["Alt_AA", "Position"]]
        .dropna(how="any")
        .to_dict("records")
    )

    for record in records:
        try:
            position = int(record["Position"]) - 1
            residues[position] = AA_3TO1[record["Alt_AA"]]
        except (KeyError, IndexError, ValueError):
            logger.warning("Skipping unmappable mutation record: %s", record)
            continue

    return "".join(residues)


def one_letter_to_three(amino_acid: str) -> str:
    """Convert a single-letter amino acid code (e.g. 'A') to three-letter ('Ala')."""
    return AA_1TO3[amino_acid]


def three_letter_to_one(amino_acid: str) -> str:
    """Convert a three-letter amino acid code (e.g. 'Ala') to single-letter ('A')."""
    return AA_3TO1[amino_acid]
