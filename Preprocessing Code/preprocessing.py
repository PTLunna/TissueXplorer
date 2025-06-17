#  --------------------- PACKAGES ---------------------
# For preprocessing
import pandas as pd
import os
import numpy as np

# -----------------------Preprocessing------------------------------


# Function to preprocess BGee dataframes
def preprocess_BGee():
    # For all the species in the Bgee folder
    for species in os.listdir("..//Datasets/Bgee"):
        species_path = os.path.join("..//Datasets/Bgee", species)
        for file in os.listdir(species_path):
            # Takes all .tsv files, which is the data from Bgee
            if file.endswith(".tsv"):
                file_path = os.path.join(species_path, file)
                # Read .tsv file and transform into a pandas dataframe
                df = pd.read_csv(file_path, sep="\t")
                # Eliminate unecessary collumns
                df = df.drop(
                    columns=[
                        "Experiment ID",
                        "Library ID",
                        "Library type",
                        "Anatomical entity ID",
                        "Stage ID",
                        "Stage name",
                        "Sex",
                        "Strain",
                        "Read count",
                        "Rank",
                        "Detection flag",
                        "pValue",
                        "State in Bgee",
                    ]
                )

                # Creates a new dataframe called df_pivot where each tissue in the collumn "Anatomical entity name" of the df dataframe becomes a new collumn filled by the corresponding "TPM" values.
                df_pivot = df.pivot(columns="Anatomical entity name", values="TPM")
                # Joins the two dataframes in a new one called df_final
                df_final = pd.concat([df, df_pivot], axis=1)
                # Removes "Anatomical entity name" and "TPM" collumns
                df_final = df_final.drop(columns=["Anatomical entity name", "TPM"])
                # PREPROCESSING
                # Substituting missing values for 0
                df_final.fillna(0, inplace=True)
                # Removing lines where all collumns except "Gene ID" are 0
                df_final = df_final.loc[
                    ~(df_final.drop(columns=["Gene ID"], errors="ignore") == 0).all(
                        axis=1
                    )
                ]
                # Adding a collumn called "GeneSymbol"
                df_final.insert(1, "GeneSymbol", None)
                # Saving df_final in a .csv file
                output_file = os.path.join(
                    species_path, f"PP_Bgee_{file.replace('.tsv', '.csv')}"
                )
                df_final.to_csv(output_file, index=False)
                print(f"Complete for {species}")


# Function to preprocess Expression Atlas dataframes
def preprocess_ExpAtlas():
    # For all the species in the Expression Atlas folder
    for species in os.listdir("..//Datasets/Expression Atlas"):
        species_path = os.path.join("..//Datasets/Expression Atlas", species)
        for file in os.listdir(species_path):
            # Takes all .tsv files, which is the data from Expression Atlas
            if file.endswith(".tsv"):
                file_path = os.path.join(species_path, file)
                # Read .tsv file and transform into a pandas dataframe
                df = pd.read_csv(file_path, sep="\t", header=4)
                # Substite missing values with 0
                df.fillna(0, inplace=True)
                # Change the name of the collumn Gene Name to "GeneSymbol"
                df.rename(columns={"Gene Name": "GeneSymbol"}, inplace=True)
                # Delete geneIDs in the Gene Symbol collumn
                df.loc[
                    df["GeneSymbol"].astype(str).str.contains("00000"), "GeneSymbol"
                ] = np.nan
                # Save in a csv file
                output_file = os.path.join(
                    species_path,
                    f"PP_ExpAtlas_{species}_{file.replace('.tsv', '.csv')}",
                )
                df.to_csv(output_file, index=False)
                print(f"Complete for {species}")


# -------------------- TESTING -----------

# preprocess_ExpAtlas()
# preprocess_BGee()
