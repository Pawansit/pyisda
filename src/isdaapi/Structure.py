import requests
import pandas as pd

def StructureDetails(uniprot_id, PDBDetails=False):
    
    #url = f"https://ibdc.dbtindia.gov.in/isda/api/protein_detail/{uniprot_id}/?additional_outputs=structuralMapping"
    url = f"http://10.74.0.47:8000/api/protein_detail/{uniprot_id}/?additional_outputs=structuralMapping"

    try:
        # Make the GET request to the API
        response = requests.get(url, timeout=10) # Set a timeout for the request
        response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)
        
        # Parse the JSON response into a Python dictionary
        data = response.json()
        ProtLength =  data.get("sequence", {}).get("length")
        if PDBDetails == True:
            extracted_data = {
                "PDB Count": data.get("structuralMapping",{}).get("pdb_count",{}),
                "Chain Count": data.get("structuralMapping", {}).get("chains_count",{}),
                "PDBRecords" : data.get("structuralMapping", {}).get("pdb_data",{}),
                "HighestCoverageID" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("pdb_id","NA"),
                "HighestCoverage" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("value","NA"),
                "HighestCoveragePDB_Chain" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("pdb_chain","NA"),
                "HighestCoverageAuth_Chain" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("auth_chain","NA")
            }
            return extracted_data
        else:
            extracted_data = {
                "PDB Count": data.get("structuralMapping",{}).get("pdb_count",{}),
                "Chain Count": data.get("structuralMapping", {}).get("chains_count",{}),
                "HighestCoverageID" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("pdb_id","NA"),
                "HighestCoverage" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("value","NA"),
                "HighestCoveragePDB_Chain" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("pdb_chain","NA"),
                "HighestCoverageAuth_Chain" : data.get("structuralMapping", {}).get("highest_coverage",{}).get("auth_chain","NA")
            }
            return extracted_data
    except requests. exceptions.RequestException as e:
        return {"error": f"An error occurred during the request: {e}"}
    except ValueError as e:
        return {"error": f"Error decoding JSON: {e}"}
    

def Mutated_Stru(MutationTable_df: pd.DataFrame, pdb_id: str, auth_chain_id: str):
    """
    Generates a ChimeraX script to apply multiple mutations to a structure 
    and saves the script for execution.

    Args:
        MutationTable_df: DataFrame with 'Alt_AA' (new residue) and 'Position' (residue number).
        pdb_id: PDB ID of the structure (used for filenames).
        auth_chain_id: Use the Auth chain ID The specific chain ID to mutate.
    """

    script_filename = f"ChimeraXSession_{pdb_id}.cxc"
    input_cif_filename = f"{pdb_id}.cif"
    output_cif_filename = f"{pdb_id}_mutated_structure.cif"

    chain_df = MutationTable_df[MutationTable_df["pdb_auth_chain"] == auth_chain_id]

    with open(script_filename, "w") as P:
        print(f"open {input_cif_filename}", file=P)
        
        for index, row in chain_df.iterrows():
            alt_aa = row['Alt_AA']
            position = row['pdb_residue']

            if pd.isna(alt_aa) or alt_aa == "" or alt_aa == "nan":
                continue
            
            # swapaa NEW_RESIDUE :RES_NUM.CHAIN_ID
            print(f"swapaa /{auth_chain_id}: {position}  {alt_aa}", file=P)

        # 3. Save the result
        print(f"save {output_cif_filename}", file=P)
        # Add a quit command so it exits when run headlessly
        print("cli quit", file=P)
    
    print(f"ChimeraX script '{script_filename}' generated successfully.")
