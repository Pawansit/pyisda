import requests
import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np




def AFMutationTable(uniprot_id: str, ranges: str|None, significance_type : str|None) -> pd.DataFrame | None:
    """
    Retrieves computational mutation data for a given Uniprot ID from the ISDA API, 
    processes it into a DataFrame.

    Args:
        uniprot_id: The Uniprot ID for the query.
        source: The data source (e.g., 'clinvar').
        record_type: The type of record (e.g., 'default').
        selection: A list of columns to select. If None, returns all columns.

    Returns:
        A pandas DataFrame with mutation data, or None if an error occurs.
    """
    
    url = f"http://10.74.0.47:8000/api/computational_mutations/{uniprot_id}"


    if ranges and significance_type:
        url = url+f"/?ranges={ranges}&significance_type={significance_type}"
    elif ranges:
        url =url+f"/?ranges={ranges}"
    elif significance_type:
        url = url+f"?significance_type={significance_type}"

    AApattern = r"^([A-Z])(\d+)([A-Z])$"

    try:
        response = requests.get(url, timeout=10) 
        response.raise_for_status() 

        data = response.json()


        
        if data is None:
            print(f"Could not retrieve or process data")
            return None
            
        
        df = pd.DataFrame(data)

        df[['Ref_AA', 'Position', 'Alt_AA']] = df['protein_variant'].str.extract(AApattern)

        return(df)

    except requests.exceptions.Timeout:
        print("The request timed out. The server took too long to respond.")
    except requests.exceptions.ConnectionError:
        print("A connection error occurred. Check your network or the URL.")
    except requests.exceptions.HTTPError as e:
        print(f"HTTP error occurred: {e.response.status_code} {e.response.reason}")
    except requests.exceptions.JSONDecodeError:
        print("Failed to decode JSON response. The API might not have returned valid JSON.")
    except KeyError as e:
        print(f"KeyError: Missing expected key in the API response: {e}") 
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        
    return None