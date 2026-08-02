"""
pyisda
======

A local toolkit for protein structural analysis built around the IBDC
ISDA REST API: protein/sequence lookup, UniProt<->PDB residue mapping,
structural coverage, mutation retrieval & analysis, ChimeraX mutation
script generation, and per-residue SASA calculation.

Quick start
-----------
>>> import pyisda
>>> pyisda.get_protein_details("P00533")
>>> pyisda.get_residue_map("P00533", "1m17")
>>> pyisda.get_mutation_table("P00533")
"""

from . import _client as config  # noqa: F401  (exposed for ISDA_BASE_URL overrides)
from ._client import ISDARequestError

from .protein import (
    get_protein_details,
    get_sequence,
    mutate_sequence,
    one_letter_to_three,
    three_letter_to_one,
)
from .structure import (
    get_structure_details,
    get_pdb_summary,
    generate_mutated_structure_script,
)
from .mapper import get_residue_map
from .mutation import (
    get_mutation_table,
    get_computational_mutations,
    filter_by_significance,
    analyze_mutation_properties,
    get_mutation_summary,
    merge_mutations_with_pdb_map,
    merge_experimental_with_computational,
    subset_by_positions,
    plot_mutation_matrix,
)
from .sasa import calculate_sasa_for_chain, fetch_structure_cif, fetch_structure_array, SASAResult
from .ligand import get_bound_ligands, get_binding_site_residues
from .visualize import (
    plot_structural_coverage,
    plot_mutation_lollipop,
    plot_sasa_profile,
    plot_protein_overview,
)
from .interactive import (
    plot_structural_coverage_html,
    plot_mutation_lollipop_html,
    plot_sasa_profile_html,
    plot_protein_overview_html,
    save_html,
)
from .viewer import (
    build_chain_view,
    show_structure_chains,
    save_structure_html,
)

__version__ = "0.1.0"

__all__ = [
    "ISDARequestError",
    "get_protein_details",
    "get_sequence",
    "mutate_sequence",
    "one_letter_to_three",
    "three_letter_to_one",
    "get_structure_details",
    "get_pdb_summary",
    "generate_mutated_structure_script",
    "get_residue_map",
    "get_mutation_table",
    "get_computational_mutations",
    "filter_by_significance",
    "analyze_mutation_properties",
    "get_mutation_summary",
    "merge_mutations_with_pdb_map",
    "merge_experimental_with_computational",
    "subset_by_positions",
    "plot_mutation_matrix",
    "calculate_sasa_for_chain",
    "fetch_structure_cif",
    "fetch_structure_array",
    "SASAResult",
    "get_bound_ligands",
    "get_binding_site_residues",
    "plot_structural_coverage",
    "plot_mutation_lollipop",
    "plot_sasa_profile",
    "plot_protein_overview",
    "plot_structural_coverage_html",
    "plot_mutation_lollipop_html",
    "plot_sasa_profile_html",
    "plot_protein_overview_html",
    "save_html",
    "build_chain_view",
    "show_structure_chains",
    "save_structure_html",
]
