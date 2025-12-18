import requests
import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np

def MutationTable(uniprot_id: str, source: str = 'clinvar', record_type: str = "default", selection: list = None) -> pd.DataFrame | None:
    """
    Retrieves mutation data for a given Uniprot ID from the ISDA API, 
    processes it into a DataFrame, and filters columns based on selection.

    Args:
        uniprot_id: The Uniprot ID for the query.
        source: The data source (e.g., 'clinvar').
        record_type: The type of record (e.g., 'default').
        selection: A list of columns to select. If None, returns all columns.

    Returns:
        A pandas DataFrame with mutation data, or None if an error occurs.
    """
    
    #url = f"https://ibdc.dbtindia.gov.in/isda/api/mutation/?uniprot_id={uniprot_id}&sources={source}&output_type={record_type}"
    url = f"http://10.74.0.47:8000/api/mutation/?uniprot_id={uniprot_id}&sources={source}&output_type={record_type}"
    AApattern = r"p\.(?P<Ref_AA>[A-Z][a-z]{2})(?P<Position>\d+)(?P<Alt_AA>[A-Z][a-z]{2}|\*)"
    

    try:
        response = requests.get(url, timeout=10) 
        response.raise_for_status() 

        data = response.json()
        
        if 'mutation' not in data or not data['mutation']:
            print("No mutation records found in the API response data.")
            return None
            
        df = pd.DataFrame(data['mutation'])
        
        df[['Ref_AA', 'Position', 'Alt_AA']] = df['Variant AA Change'].str.extract(AApattern)
        
        if selection is None:
            return df
        elif isinstance(selection, list):
            df = df[selection]
            return df
        else:
            if all(col in df.columns for col in selection):
                return df[selection]
            else:
                missing_cols = [col for col in selection if col not in df.columns]
                print(f"Error: Missing columns in selection: {missing_cols}")
                return None

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


def MatchSignificance(df: pd.DataFrame, pattern: str, case: bool = False) -> pd.DataFrame:
    """
    Filters a pandas DataFrame based on whether the 'Clinical_Significance' column 
    contains a specific regex pattern.

    Args:
        df: The input DataFrame containing a 'Clinical_Significance' column.
        pattern: The regex pattern (string) to search for.
        case: If True, the pattern matching is case-sensitive; if False (default), it is case-insensitive.

    Returns:
        A new DataFrame containing only the rows that match the pattern, 
        with the index reset.
    
    Raises:
        ValueError: If the input DataFrame is missing the required 'Clinical_Significance' column.
    """

    if 'Clinical Significance' not in df.columns:
        raise ValueError("Input DataFrame must contain a 'Clinical Significance' column.")

    if not isinstance(pattern, str) or not pattern.strip():
         raise ValueError("A valid search pattern (string) must be provided.")

    filtered_df = df.copy()

    filtered_df = filtered_df[
        filtered_df['Clinical Significance'].str.contains(pattern, case=case, na=False)
    ].reset_index(drop=True)

    return filtered_df


def get_change_type(row: pd.Series) -> str:
    """Helper function to determine the type of mutation effect based on properties."""
    ref = row['Ref_Prop']
    alt = row['Alt_Prop']
    
    # Use pandas/numpy checks for missing values
    if pd.isna(ref) or pd.isna(alt):
        return "Unknown/N/A"
    
    if ref == alt:
        return f"No change ({ref})"
    else:
        return f"Change from {ref} to {alt}"
    

def AnalyzeMutationProperties(
    df: pd.DataFrame, 
    hydrophobic_aas: list = ['Ala', 'Val', 'Leu', 'Ile', 'Pro', 'Phe', 'Trp', 'Met', 'Cys'], 
    hydrophilic_aas: list = ['Gly', 'Ser', 'Thr', 'Tyr', 'Asn', 'Gln', 'Asp', 'Glu', 'Lys', 'Arg', 'His', '*'], 
    category_labels: list = ["Hydrophobic", "Hydrophilic"]
) -> pd.DataFrame | None:
    """
    Analyzes mutation properties (e.g., hydrophobic/hydrophilic change) for a given protein ID.

    Args:
        ProteinID: The UniProt ID used to fetch data via MutationTable.
        hydrophobic_aas: List of amino acids considered hydrophobic.
        hydrophilic_aas: List of amino acids considered hydrophilic or 'other'.
        category_labels: Labels corresponding to the two AA lists.
    
    Returns:
        A filtered pandas DataFrame with new property columns, or None if data fetching fails.
    """
    
    if df is None or df.empty:
        print(f"Could not retrieve or process data Frame")
        return None
    
    prop_map = {aa: category_labels[0] for aa in hydrophobic_aas}
    prop_map.update({aa: category_labels[1] for aa in hydrophilic_aas})
    
    df['Ref_Prop'] = df['Ref_AA'].map(prop_map).fillna('Unknown')
    df['Alt_Prop'] = df['Alt_AA'].map(prop_map).fillna('Unknown')
    
    df['Mutation_Effect_Type'] = df.apply(get_change_type, axis=1)

    return df



def GetMutationSummary(
    df: pd.DataFrame, 
    top_n: int = 10
) -> dict | None:
    """
    Analyzes a DataFrame of mutation data to provide summary statistics.

    Args:
        df: The input DataFrame expected to have 'Clinical_Significance' and 'Position' columns.
        top_n: The number of top mutation spots to return in the summary.

    Returns:
        A dictionary summarizing mutation counts, significance distribution, 
        and top mutation positions, or None if the DataFrame is invalid.
    """
    if df is None or df.empty:
        print("Warning: Input DataFrame is empty or None.")
        return None
        
    required_cols = ['Clinical Significance', 'Position']
    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        print(f"Error: DataFrame is missing required columns: {missing}")
        return None

    extracted_data = {
        "mutation_count": df.shape[0],
        "significance_status": [df['Clinical Significance'].value_counts().to_dict()],
        "Mutation_type": [df['Consequence Type'].value_counts(sort=True).to_dict()],
        "maximum_mutation_spot": [df['Position'].value_counts(sort=True).head(top_n).to_dict()]
    }
    
    return extracted_data



def plot_mutation_matrix(df: pd.DataFrame):
    """
    Generates a heatmap/matrix visualization of mutation counts per position and
    alternative amino acid.

    Args:
        df: A pandas DataFrame containing 'Position', 'Ref_AA', and 'Alt_AA' columns.
    """
    if df is None or df.empty:
        print("Cannot plot; DataFrame is empty or None.")
        return

    # Pivot the data to create a matrix of counts
    # Rows: Alternative Amino Acids, Columns: Positions
    matrix_df = pd.crosstab(index=df['Alt_AA'], columns=df['Position'])

    # Ensure consistent ordering for common AAs in the plot if possible
    # (The order will primarily be determined by the data present)
    
    plt.figure(figsize=(15, 8))
    plt.imshow(matrix_df.values, cmap='viridis', aspect='auto')
    
    # Set up the axis labels and ticks
    plt.xticks(ticks=np.arange(len(matrix_df.columns)), labels=matrix_df.columns, rotation=90, fontsize=8)
    plt.yticks(ticks=np.arange(len(matrix_df.index)), labels=matrix_df.index, fontsize=10)
    
    plt.xlabel('Mutation Position (AA)', fontsize=12)
    plt.ylabel('Alternative Amino Acid (Alt_AA)', fontsize=12)
    plt.title(f'Mutation Matrix for Uniprot ID: {df["AC"].iloc[0]}')      ### Check here 
    plt.colorbar(label='Mutation Count')
    
    # Add count annotations to the cells (optional, can be messy with lots of data)
    # for i in range(len(matrix_df.index)):
    #     for j in range(len(matrix_df.columns)):
    #         count = matrix_df.iloc[i, j]
    #         if count > 0:
    #             plt.text(j, i, count, ha='center', va='center', color='white', fontsize=6)

    plt.tight_layout()
    plt.show()


# def DfMerge(Mutationtable_df,UniprotPDBMapper_df):

#     """
#     Maps mutation sites from Mutationtable_df to UniprotPDBMapper_df 
#     using pandas vectorization.
#     """
    
#     Strutmutation = {}
#     for M in Mutationtable_df['Position'].to_list():
#         Strutmutation[M] = "Known Mutation Site"
    
#     # Use the map function directly. Pandas handles alignment efficiently.
#     UniprotPDBMapper_df['Mutation Site'] = UniprotPDBMapper_df['unp_residue'].map(Strutmutation)
#     UniprotPDBMapper_df = UniprotPDBMapper_df.dropna(how='any').reset_index()
    
#     return UniprotPDBMapper_df


def DfMerge(Mutationtable_df,UniprotPDBMapper_df):

    """
    Maps mutation sites from Mutationtable_df to UniprotPDBMapper_df 
    using pandas vectorization.
    """

    UniprotPDBMapper_df.columns.values[0] = "Position"

    merged_df = pd.merge(Mutationtable_df, UniprotPDBMapper_df, on='Position')
    return(merged_df)

