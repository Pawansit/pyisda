#import numpy as np
#from biotite.structure.io import pdb
#from biotite.structure import sasa
#import biotite.structure as struc
#import os




import numpy as np
import os
from typing import Optional, Tuple
from biotite.structure.io import pdb
import biotite.structure as struc
import biotite.database.rcsb as rcsb

def calculate_sasa_for_chain(
    pdb_id: str, 
    chain_id: str, 
    probe_radius: float =  1.4,
    local_pdb_path: Optional[str] = None
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Calculates the per-residue solvent-accessible surface area (SASA) 
    for a specific PDB ID and chain.

    Args:
        pdb_id: The 4-character PDB accession ID (e.g., '1l2y').
        chain_id: The specific chain identifier (e.g., 'A').
        local_pdb_path: Optional path to a local PDB file. If None, the 
                        function fetches the file from the RCSB PDB database.

    Returns:
        A tuple of (residue_ids, residue_names, sasa_per_residue) as NumPy arrays,
        or None if an error occurs.
    """
    try:
        if local_pdb_path:
            # Use local file if provided
            file_path = local_pdb_path
        else:
            # Fetch from RCSB PDB database dynamically
            file_path = rcsb.fetch(pdb_id, "pdb", os.getcwd()) # Saves to current working directory

        file = pdb.PDBFile.read(file_path)
        
        array = file.get_structure(model=1) 

    except Exception as e:
        print(f"Error reading PDB file {pdb_id}: {e}")
        return None

    if chain_id not in array.chain_id:
        print(f"Error: Chain ID '{chain_id}' not found in PDB structure {pdb_id}.")
        return None
        
    chain = array[array.chain_id == chain_id]
    is_amino_acid = struc.filter_amino_acids(chain)
    protein_chain = chain[is_amino_acid]
    
    # 3. Calculate atom-wise SASA
    sasa_per_atom = struc.sasa(protein_chain, vdw_radii='ProtOr', probe_radius=probe_radius) 
    
    # 4. Sum up SASA for each residue (residue-wise calculation)
    sasa_per_residue = struc.apply_residue_wise(protein_chain, sasa_per_atom, np.sum)
    
    # 5. Extract residue identifiers for mapping results
    ids, names = struc.get_residues(protein_chain)

    return ids, names, sasa_per_residue
