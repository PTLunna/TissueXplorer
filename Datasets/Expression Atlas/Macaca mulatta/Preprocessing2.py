#  --------------------- PACKAGES ---------------------
import pandas as pd
import os

# For APIs
import requests
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# -----------------Normalizing gene symbols------------

# API to transform gene IDs into gene symbols
BATCH_SIZE = 100
MAX_WORKERS = 10  # Número de chamadas simultâneas


def get_gene_symbol(ensembl_id):
    "Calls API Ensembl to obtain the Gene Symbol of a single Gene ID."
    url = (
        f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
    )
    response = requests.get(url, timeout=120)
    if response.status_code == 200:
        data = response.json()
        return ensembl_id, data.get("display_name", "Unknown")
    else:
        return ensembl_id, "Not Found"


def get_gene_symbol_batch(ensembl_ids):
    "Executes múltiple API calls in paralell using ThreadPoolExecutor."
    gene_symbols = {}
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        results = executor.map(get_gene_symbol, ensembl_ids)
    for gene_id, symbol in results:
        gene_symbols[gene_id] = symbol
    return gene_symbols


def add_gene_names(file_path, species):
    df = pd.read_csv(file_path)
    tqdm.pandas()
    # Filters lines where GeneSymbol is empty or Na
    missing_mask = df["GeneSymbol"].isna() | (
        df["GeneSymbol"].astype(str).str.strip() == ""
    )
    missing_symbols = df[missing_mask]

    # Removes duplicates before calling the API
    unique_gene_ids = missing_symbols["Gene ID"].drop_duplicates().tolist()
    gene_symbols = {}  # Dictionary to store unique results
    # Processing in batches
    for i in tqdm(
        range(0, len(unique_gene_ids), BATCH_SIZE), desc="Processing in batches"
    ):
        batch_ids = unique_gene_ids[i : i + BATCH_SIZE]
        batch_results = get_gene_symbol_batch(batch_ids)
        gene_symbols.update(batch_results)  # Save in cache

    df["GeneSymbol"] = df["GeneSymbol"].astype("object")

    # Substituting only the gene symbols that we did not know
    df.loc[missing_mask, "GeneSymbol"] = df.loc[missing_mask, "Gene ID"].map(
        gene_symbols
    )
    df.to_csv(f"{species}_with_gene_symbols.csv", index=False)


# --------------------Final preprocessing ------


# Transforming "Unknown" or "Not Found" in GeneSymbol collumn in Gene IDs and delete the collumn Gene ID
def back_geneid():
    # Obtains the list of files in the current folder that end with "_with_gene_symbols.csv"
    files = [f for f in os.listdir() if f.endswith("_with_gene_symbols.csv")]

    for file in files:
        df = pd.read_csv(file)

        if "GeneSymbol" in df.columns and "Gene ID" in df.columns:
            # Replace "Unknown" or "Not Found" with Gene ID
            mask = df["GeneSymbol"].isin(["Unknown", "Not Found"])
            df.loc[mask, "GeneSymbol"] = df.loc[mask, "Gene ID"]

            # Drop the 'Gene ID' column
            df.drop(columns=["Gene ID"], inplace=True)
            # Naming and saving final .csv file
            new_filename = f"final_{file}"
            df.to_csv(new_filename, index=False)
            print(f"File saved: {new_filename}")
        else:
            print(f"Error: {file} Does not contain the necessary collumns")


# -------------------- TESTING -----------
# Uncomment the following functions to run the code
# --


# add_gene_names(
#     "PP_ExpAtlas_Macaca mulatta_E-MTAB-2799-query-results.tpms.csv",
#     "Macaca mulatta",
# )

# back_geneid()
