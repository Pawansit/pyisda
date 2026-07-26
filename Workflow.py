import pyisda as ibdc 
import pandas as pd

# Protein-level metadata
details = ibdc.get_protein_details("P00533")   # EGFR
seq = ibdc.get_sequence("P00533")

# Structural coverage
coverage = ibdc.get_structure_details("P00533", include_pdb_records=True)

# UniProt <-> PDB residue mapping
residue_map = ibdc.get_residue_map("P00533", "1m17")

# Mutation retrieval + analysis
clinvar_df = ibdc.get_mutation_table("P00533", source="clinvar",record_type="Full")
annotated = ibdc.analyze_mutation_properties(clinvar_df)
summary = ibdc.get_mutation_summary(clinvar_df)

# Entry-level PDB summary
pdb_info = ibdc.get_pdb_summary("6q0j", additional_outputs="micromolecular_data")
print(pd.DataFrame(pdb_info['MicroMolecule']))

