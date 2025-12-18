import requests
import os 
import pandas as pd


def Genes_Details(data_structure: dict):
    """
    Extracts all gene values from a nested dictionary structure, 
    handling multiple entries in the 'genes' list.

    Args:
        data_structure (dict): The input dictionary.

    Returns:
        list: A list of extracted gene values.
    """
    all_gene_values = []
    
    # Safely get the 'genes' list, default to an empty list if key is missing
    genes_list = data_structure.get('genes', [])
    
    # Iterate over each item in the genes list
    for gene_entry in genes_list:
        try:
            # Safely navigate the nested structure within each entry
            gene_value = gene_entry.get('geneName', {}).get('value')
            
            # If a value was found (is not None), add it to our list
            if gene_value:
                all_gene_values.append(gene_value)
                
        except (TypeError, AttributeError):
            # This handles cases where 'geneName' isn't a dict, etc.
            continue # Skip to the next item in the loop if an error occurs

    return all_gene_values

def Genes_Synonyms_Details(data_structure: dict):
    """
    Extracts all gene values from a nested dictionary structure, 
    handling multiple entries in the 'genes' list.

    Args:
        data_structure (dict): The input dictionary.

    Returns:
        list: A list of extracted gene values.
    """
    all_alter_values = []
    
    # Safely get the 'genes' list, default to an empty list if key is missing
    genes_list = data_structure.get('genes', [])
    
    
    # Iterate over each item in the genes list
    for gene_entry in genes_list:
        try:
            #print(gene_entry.get("synonyms"))
            # Safely navigate the nested structure within each entry
            gene_value = gene_entry.get('synonyms', {})
            for GValue in gene_value:
                all_alter_values.append(GValue.get("value",{}))        
        except (TypeError, AttributeError):
            # This handles cases where 'geneName' isn't a dict, etc.
            continue # Skip to the next item in the loop if an error occurs

    return all_alter_values



def get_protein_details(uniprot_id:str):
    """
    Fetches protein details from the IBDC ISDA API and extracts specific information.

    Args:
        uniprot_id (str): The UniProt ID of the protein.

    Returns:
        dict: A dictionary containing the extracted details, or an error message.
    """
    #url = f"https://ibdc.dbtindia.gov.in/isda/api/protein_detail/{uniprot_id}/?additional_outputs=proteinFeatures,structuralMapping"
    url = f"http://10.74.0.47:8000/api/protein_detail/{uniprot_id}/?additional_outputs=proteinFeatures,structuralMapping"
    
    try:
        # Make the GET request to the API
        response = requests.get(url, timeout=10) # Set a timeout for the request
        response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)
        
        # Parse the JSON response into a Python dictionary
        data = response.json()
        
        # Extract the required information using dictionary keys
        extracted_data = {
            "uniprot_id": data.get("uniprot_id"),
            "scientificName": data.get("organism", {}).get("scientificName"),
            "taxonId": data.get("organism", {}).get("taxonId"),
            "Protein recommended name": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName",{}).get("value",{}),
            "Gene Name": (",").join(Genes_Details(data)),
            "Gene synonyms": (",").join(Genes_Synonyms_Details(data)),
            "StructuralRecords": data.get("structuralMapping",{}).get("pdb_count",0),
            "SequenceLength": data.get("sequence",{}).get("length",{})
        }

        return extracted_data

    except requests. exceptions.RequestException as e:
        return {"error": f"An error occurred during the request: {e}"}
    except ValueError as e:
        return {"error": f"Error decoding JSON: {e}"} 
    
def Sequence(uniprot_id:str):
    """
    Fetches protein details from the IBDC ISDA API and extracts specific information.

    Args:
        uniprot_id (str): The UniProt ID of the protein.

    Returns:
        dict: A dictionary containing the extracted details, or an error message.
    """
    #url = f"https://ibdc.dbtindia.gov.in/isda/api/protein_detail/{uniprot_id}/?additional_outputs=proteinFeatures"
    url  = f"http://10.74.0.47:8000/api/protein_detail/{uniprot_id}/?additional_outputs=proteinFeatures"

    try:
        # Make the GET request to the API
        response = requests.get(url, timeout=10) # Set a timeout for the request
        response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)
        
        # Parse the JSON response into a Python dictionary
        data = response.json()

        extracted_data = {
            "Sequence": data.get("sequence").get("value"),
            "Length": data.get("sequence", {}).get("length"),
            "Mol Weight": data.get("sequence", {}).get("molWeight")
        }
        return extracted_data

    except requests. exceptions.RequestException as e:
        return {"error": f"An error occurred during the request: {e}"}
    except ValueError as e:
        return {"error": f"Error decoding JSON: {e}"} 
    
def MutateSequence(Sequence:str, Mutations: pd.DataFrame):

    AA_CODES = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    '*': '*'}
  
    Seq_list = list(Sequence.get("Sequence", "Provide the Sequence"))
    print(f"Original Sequence: {''.join(Seq_list)}")

    Mutations = Mutations[['Alt_AA','Position']].dropna(how='any').reset_index()
    Mutation_list = Mutations.to_dict("records")


    for Record in Mutation_list:
        try: 
            Position = int(Record['Position']) - 1
            AA = AA_CODES[Record['Alt_AA']]
            Seq_list[Position] = AA
        except KeyError:
            continue

    print(f"Mutated Sequence: {''.join(Seq_list)}")


def Codes_1to3(AA:str):
    AA_CODES = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': '*'}

    return(AA_CODES[AA])