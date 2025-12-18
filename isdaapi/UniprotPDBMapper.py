import requests
import json
from typing import List, Dict, Optional

def GetResidueMap(uniprot_id: str, pdb_id: str) -> Optional[List[Dict[str, str]]]:
    """
    Retrieves a mapping between UniProt residue numbers and PDB residue numbers 
    for a specific PDB structure associated with a UniProt ID.

    Args:
        uniprot_id: The UniProt accession ID (e.g., 'P00533').
        pdb_id: The 4-character PDB ID (e.g., '1XYZ').

    Returns:
        A list of dictionaries containing the residue mappings, 
        or None if data cannot be retrieved or processed.
    """
    
    # Ensure PDB ID is correctly formatted for lookup
    pdb_id_upper = pdb_id.upper()
    #url = f"https://ibdc.dbtindia.gov.in/isda/api/protein_detail/{uniprot_id}/?additional_outputs=structuralMapping"
    url = f"http://10.74.0.47:8000/api/protein_detail/{uniprot_id}/?additional_outputs=structuralMapping"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        
        data = response.json()
        
        # --- Robustness: Check for nested keys explicitly ---
        try:
            pdb_records = data['structuralMapping']['pdb_data']
        except KeyError:
            print("Error: 'structuralMapping' or 'pdb_data' key missing from API response.")
            return None

        residue_map = []

        for mapping in pdb_records:
            if mapping.get('PDB') == pdb_id_upper:
                # Type conversion for range parameters is handled once
                unp_start = int(mapping["UNIPROT_START"])
                unp_end = int(mapping["UNIPROT_END"])
                
                # Try converting PDB_START to integer early, handle potential errors immediately
                try:
                    pdb_start = int(mapping["PDB_START"])
                except (ValueError, TypeError):
                    print(f"Warning: Non-integer PDB_START value found for mapping: {mapping['PDB_START']}. Skipping range calculation.")
                    continue # Skip this mapping block

                # Use zip for clearer iteration over simultaneous ranges
                unp_range = range(unp_start, unp_end + 1)
                pdb_range = range(pdb_start, pdb_start + len(unp_range))

                for unp_residue, pdb_residue in zip(unp_range, pdb_range):
                    map_residue = {
                        "unp_residue": str(unp_residue), 
                        "pdb_residue": str(pdb_residue), 
                        "pdb_auth_chain": str(mapping.get('AUTH_CHAIN', 'N/A')),
                        "pdb_chain" : str(mapping.get('PDB_CHAIN', 'N/A'))
                    }
                    residue_map.append(map_residue)
        
        if not residue_map and pdb_records:
            print(f"No mapping found for PDB ID {pdb_id_upper} within the available records.")

        return residue_map

    # --- Comprehensive Error Handling ---
    except requests.exceptions.Timeout:
        print(f"Request timed out while fetching data for {uniprot_id}.")
    except requests.exceptions.ConnectionError:
        print(f"A connection error occurred for URL: {url}")
    except requests.exceptions.HTTPError as e:
        print(f"HTTP error occurred: {e.response.status_code} {e.response.reason}")
    except json.JSONDecodeError:
        print(f"Failed to decode JSON response for {uniprot_id}.")
    except Exception as e:
        # Catch any other unexpected exceptions
        print(f"An unexpected error occurred: {e}")
        
    # Ensure the function always returns None if an exception occurs
    return None













