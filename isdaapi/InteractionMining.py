import io
import pandas as pd
from arpeggio.core import InteractionComplex 
from collections import defaultdict
from pdbeccdutils.core.boundmolecule import infer_bound_molecules
import glob
import json
import requests
import gzip
import os
from typing import List, Dict, Union, Optional

def analyze_protein_interaction(struc_id: str, bound_molecule_obj) -> list:
    """
    Analyzes protein interactions for a given structure ID and a specific 
    bound molecule using the arpeggio library.

    Args:
        struc_id: PDB ID or path to the structure file.
        bound_molecule_obj: An object that can be converted to arpeggio selection format.

    Returns:
        A list of contact dictionaries returned by Complex.get_contacts().
    """
    # Use standard Python naming conventions
    complex_obj = InteractionComplex(struc_id) 
    
    # Run setup steps
    # complex_obj.write_hydrogenated("./", struc_id) # Commented out as in original
    complex_obj.structure_checks()
    complex_obj.address_ambiguities()
    complex_obj.initialize()

    # Define interaction parameters (can be made input parameters if needed)
    interaction_cutoff = 5.0
    compensation_factor = 0.1
    include_neighbors = False

    # Convert the bound molecule object to the format Arpeggio expects
    selection = bound_molecule_obj.to_arpeggio()

    complex_obj.run_arpeggio(
        selection, 
        interaction_cutoff, 
        compensation_factor, 
        include_neighbors
    )

    return complex_obj.get_contacts()


def analyze_bound_comp_interaction(pdb_id_no_h: str, pdb_id_with_h: str) -> Optional[pd.DataFrame]:
    """
    Infers bound molecules for a PDB ID in mmCIF format (without Hydrogens) and analyzes 
    their interactions using the corresponding PDB ID (with Hydrogens) structure.

    Args:
        pdb_id_no_h: PDB ID without hydrogens (used for inferring bound molecules).
        pdb_id_with_h: PDB ID with hydrogens (used for interaction calculation in Arpeggio).

    Returns:
        A pandas DataFrame of filtered interaction data, or None if no data found.
    """
    
    # Use standard Python naming conventions
    discarded_ligands = ["HOH"]
    intx_fields = ["auth_asym_id", "auth_seq_id", "label_comp_id", "auth_atom_id"]

    try:
        # Infer bound molecules from the H-less structure
        bound_molecules = infer_bound_molecules(pdb_id_no_h, discarded_ligands, assembly=True)
    except Exception as e:
        print(f"Error inferring bound molecules for {pdb_id_no_h}: {e}")
        return None

    if not bound_molecules:
        print(f"No bound molecules identified for {pdb_id_no_h} (excluding HOH).")
        return None

    print(f"Identified {len(bound_molecules)} bound molecule(s). Processing interactions...")
    
    pdb_interactions = defaultdict(list)
    
    for i, bound_mol in enumerate(bound_molecules, start=1):
        bm_residues_set = set(
            f"{r['auth_asym_id']}_{r['auth_seq_id']}_{r['label_comp_id']}" 
            for r in bound_mol.to_dict().get('residues', [])
        )
        
        # Get interactions using the H-containing structure
        protein_ligand_complex_contacts = analyze_protein_interaction(pdb_id_with_h, bound_mol)
        
        for intx in protein_ligand_complex_contacts:
            intx_entities = intx["interacting_entities"]
            
            if intx_entities in ("INTER", "INTRA_SELECTION", "SELECTION_WATER"):
                
                begin_res_id = f"{intx['bgn']['auth_asym_id']}_{intx['bgn']['auth_seq_id']}_{intx['bgn']['label_comp_id']}"
                ordered_intx_labels = ('bgn', 'end') if begin_res_id in bm_residues_set else ('end', 'bgn')

                pdb_interactions["interacting_entities"].append(intx_entities)
                pdb_interactions["interaction_type"].append(intx["type"])
                pdb_interactions["distance"].append(intx["distance"])
                pdb_interactions["bm_id"].append(f"bm_{i}")
                pdb_interactions["contact_type"].append(intx["contact"]) # Note: 'contact' in intx is a list/set of types
                
                for field in intx_fields:
                    pdb_interactions[f"{field}_1"].append(intx[ordered_intx_labels[0]][field])
                    pdb_interactions[f"{field}_2"].append(intx[ordered_intx_labels[1]][field])

    # Convert results to DataFrame
    pdb_interactions_df = pd.DataFrame.from_dict(pdb_interactions)
    
    if pdb_interactions_df.empty:
        print("No interaction data found after processing.")
        return None

    # Filter out 'proximal' contacts (often just van der Waals) and keep only 'INTER' actions
    
    # Filter out rows where 'proximal' is in the list of contact types
    pdb_interactions_df = pdb_interactions_df[
        ~pdb_interactions_df['contact_type'].apply(lambda types: 'proximal' in types)
    ]
    
    # Filter to keep only "INTER" (inter-molecular) interactions
    pdb_interactions_df = pdb_interactions_df[
        pdb_interactions_df['interacting_entities'] == "INTER"
    ].reset_index(drop=True) # Reset index once at the very end

    return pdb_interactions_df

def analyze_residue_interactions(
    structure_path: str, 
    selection_list: List[str], 
    cutoff_distance: float = 5.0
) -> Optional[pd.DataFrame]:
    """
    Analyzes protein-ligand or specific residue interactions within a PDB structure 
    using Arpeggio, returning a filtered DataFrame of significant, inter-entity contacts.

    Args:
        structure_path: Path to the local PDB file (e.g., './data/1abc.pdb').
        selection_list: A list of Arpeggio selection strings (e.g., ['/A/447/']).
        cutoff_distance: The maximum distance in Angstroms for 'proximal' interactions.

    Returns:
        A pandas DataFrame containing filtered interaction details, or None if errors occur 
        or no significant interactions are found.
    """
    
    # --- Input Validation ---
    if not os.path.exists(structure_path):
        print(f"Error: Structure file not found at {structure_path}")
        return None
    
    if not selection_list:
        print("Error: Selection list cannot be empty.")
        return None

    try:
        # --- Initialize Arpeggio Complex ---
        complex_obj = InteractionComplex(structure_path)
        complex_obj.structure_checks()
        complex_obj.address_ambiguities()
        complex_obj.initialize()

        # --- Run Arpeggio Interaction Analysis ---
        # The 'interacting_cutoff' defines the maximum distance for a 'proximal' interaction
        complex_obj.run_arpeggio(
            selection_list, 
            interacting_cutoff=cutoff_distance, 
            vdw_comp=0.1, 
            include_sequence_adjacent=False
        )

        contacts_data = complex_obj.get_contacts()

        if not contacts_data:
            print(f"No contacts found for the selection: {selection_list}")
            return None
        
        # --- Process and Filter Results using Pandas ---
        df_contact = pd.DataFrame(contacts_data)
        
        # Filter 1: Exclude 'proximal' contacts (often just weak van der Waals)
        df_contact = df_contact[
            ~df_contact['contact'].apply(lambda types: 'proximal' in types)
        ]

        # Filter 2: Keep only "INTER" (inter-entity/inter-chain) interactions
        df_contact = df_contact[
            df_contact['interacting_entities'] == "INTER"
        ].reset_index(drop=True) # Reset index once at the end

        if df_contact.empty:
            print("No significant INTER-entity interactions found after filtering proximal contacts.")
            return None
        

        return df_contact

    except Exception as e:
        # General exception handling for Arpeggio-specific errors
        print(f"An error occurred during Arpeggio processing: {e}")
        return None

# Example Usage:
# df_interactions = analyze_residue_interactions(
#     structure_path="./path/to/my_structure.pdb", 
#     selection_list=['/A/447/']
# )
# print(df_interactions.head())







































# import io
# import pandas as pd
# from arpeggio.core import InteractionComplex
# from collections import defaultdict
# from pdbeccdutils.core.boundmolecule import infer_bound_molecules
# import glob
# import json
# import requests
# import gzip

# def ProteinInteraction(StruID, boundSelect):
#     Complex = InteractionComplex(StruID)
#     #Complex.write_hydrogenated("./", StruID)
#     Complex.structure_checks()
#     Complex.address_ambiguities()
#     Complex.initialize()

#     interaction_cutoff=5.0
#     compensation_factor=0.1
#     include_neighbors=False

#     selection = boundSelect.to_arpeggio()

#     Complex.run_arpeggio(selection, interaction_cutoff, compensation_factor, include_neighbors)

#     return(Complex.get_contacts())




# def BoundCompInteraction(WoHID, WHID):
#     pdb_idwoh = f"{WoHID}"
#     pdb_idwh  = f"{WHID}"
#     discarded_ligands = ["HOH"]

#     bound_molecules = infer_bound_molecules(pdb_idwoh, discarded_ligands, assembly=True)
#     print(f"List of Bound Molecules Identified")
#     print(bound_molecules[0])
#     pdb_interactions = defaultdict(list)
#     intx_fields = ["auth_asym_id", "auth_seq_id", "label_comp_id", "auth_atom_id"]
#     for i, boundMol in enumerate(bound_molecules, start=1):
#        print(boundMol.to_dict())
#        bm_residues = set(f"{residue['auth_asym_id']}_{residue['auth_seq_id']}_{residue['label_comp_id']}" for residue in boundMol.to_dict().get('residues'))
#        print(bm_residues)
#        protein_lignad_complex= ProteinInteraction(pdb_idwh, boundMol)
#        for intx in protein_lignad_complex:
#           if intx["interacting_entities"] in ("INTER", "INTRA_SELECTION", "SELECTION_WATER"):
#              for contact in intx["contact"]:
#                 pdb_interactions["interacting_entities"].append(intx["interacting_entities"])
#                 pdb_interactions["interaction_type"].append(intx["type"])
#                 pdb_interactions["distance"].append(intx["distance"])
#                 pdb_interactions["bm_id"].append(f"bm_{i}")
#                 pdb_interactions["contact_type"].append(contact)
#                 begin_res_id = f"{intx['bgn']['auth_asym_id']}_{intx['bgn']['auth_seq_id']}_{intx['bgn']['label_comp_id']}"
#                 if begin_res_id in bm_residues:
#                    ordered_intx_residues = {1:'bgn', 2: 'end'} #keeping bound-molecule as the first residue in the table
#                 else:
#                    ordered_intx_residues = {1:'end', 2: 'bgn'}
#                 for field in intx_fields:
#                    for order in ordered_intx_residues:
#                       pdb_interactions[f"{field}_{order}"].append(intx[ordered_intx_residues[order]][field])

#     PDBID = pdb_idwoh.split('/')[-1].split('.')[0]
#     pdb_interactions_df = pd.DataFrame.from_dict(pdb_interactions)
#     pdb_interactions_df = pdb_interactions_df[pdb_interactions_df['contact'].apply(lambda types: 'proximal' not in types)].reset_index()    # Df_contact[Df_contact['contact'] == ["proximal"]]
#     pdb_interactions_df = pdb_interactions_df[pdb_interactions_df['interacting_entities'] == "INTER"].reset_index()
#     return(pdb_interactions_df)



# def ResidueInteraction(pdb_id:str, selection:list):
   
#    structure_file = pdb_id
#    ##selection = ['/A/447/']
#    complex = InteractionComplex(structure_file)
#    complex.structure_checks()
#    complex.address_ambiguities()
#    complex.initialize()

#    # 4. Run Arpeggio for the specific selection
#    # The 'interacting_cutoff' defines the maximum distance for a 'proximal' interaction

#    complex.run_arpeggio(selection, interacting_cutoff=5.0, vdw_comp=0.1, include_sequence_adjacent=False)

#    # 5. Get the results
#    contacts = complex.get_contacts()

#    # 6. Analyze the results (example: count and list interaction types)
#    print(f"Total contacts found for residue A:200: {len(contacts)}")

#    # Displaying some details of the interactions

#    if contacts:
#     Df_contact = pd.DataFrame(contacts)
#     Df_contact = Df_contact[Df_contact['contact'].apply(lambda types: 'proximal' not in types)].reset_index()    # Df_contact[Df_contact['contact'] == ["proximal"]]
#     Df_contact = Df_contact[Df_contact['interacting_entities'] == "INTER"].reset_index()
#     return(Df_contact)


# #BoundCompInteraction("1uwh.cif","1uwh.cif")

# #ResidueInteraction("1uwh.cif", Selection=['/A/1723/BAX'])


# ResidueInteraction("1uwhH.cif",selection=['/A/472/'] )