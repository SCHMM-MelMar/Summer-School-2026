from pathlib import Path

import numpy as np
import pandas as pd


def normalize_negative_control_values(bioactivity_df):
    result = bioactivity_df.copy()
    negative_control_mask = (
        result["value"].fillna("").astype(str).str.strip().str.lower() == "negative control"
    )

    if negative_control_mask.any():
        result.loc[negative_control_mask, "moa"] = "Negative control"
        result.loc[negative_control_mask, "value"] = np.nan

    return result


def generate_tsvs():
    input_file = "EUbOPEN compounds_cv_approved.csv"
    print(f"Reading data from {input_file}...")
    
    # Read the CSV file once
    df = pd.read_csv(input_file, low_memory=False)
    
    # Handle slight variations in user's column naming (InChi Key vs InChiKey)
    inchikey_col = 'Compound InChi Key' if 'Compound InChi Key' in df.columns else 'Compound InChiKey'
    
    # ==========================================
    # 1. target.tsv
    # ==========================================
    target_df = pd.DataFrame()
    target_df['target_key'] = df['Target ID'] 
    target_df['name'] = np.nan
    target_df['type'] = np.nan
    
    target_df = target_df.drop_duplicates(subset=['target_key']).dropna(subset=['target_key'])
    target_df.to_csv("target.tsv", sep='\t', index=False)
    print("Successfully generated target.tsv")
    
    # ==========================================
    # 2. uniprot.tsv
    # ==========================================
    uniprot_df = pd.DataFrame()
    uniprot_df['uniprot_id'] = df['UniProt ID']
    uniprot_df['target_key'] = df['Target ID']
    uniprot_df['HGNC'] = np.nan
    uniprot_df['species'] = np.nan
    
    uniprot_df = uniprot_df.drop_duplicates(subset=['uniprot_id', 'target_key'])
    uniprot_df = uniprot_df.dropna(subset=['uniprot_id', 'target_key'], how='all')
    uniprot_df.to_csv("uniprot.tsv", sep='\t', index=False)
    print("Successfully generated uniprot.tsv")

    # ==========================================
    # 3. compound.tsv
    # ==========================================
    compound_df = pd.DataFrame()
    compound_df['inchikey'] = df[inchikey_col]
    compound_df['smiles'] = df['Compound SMILES']
    compound_df['chembl_id'] = np.nan
    compound_df['name'] = df['Virtual Compound Preferred Name']
    
    compound_df = compound_df.drop_duplicates(subset=['inchikey'])
    compound_df = compound_df.dropna(subset=['inchikey', 'smiles'], how='all')
    compound_df.to_csv("compound.tsv", sep='\t', index=False)
    print("Successfully generated compound.tsv")

    # ==========================================
    # 4. bioactivity.tsv
    # ==========================================
    # Extract only the integer from Recommended Concentration
    extracted_concentration = df['Recommended Concentration'].astype(str).str.extract(r'(\d+)', expand=False)
    
    # --- Create Biochemical DataFrame ---
    bio1 = pd.DataFrame()
    bio1['inchikey'] = df[inchikey_col]
    bio1['target_key'] = df['Target ID']
    bio1['moa'] = df['Mode of Action']
    bio1['cell_line'] = np.nan
    bio1['unit'] = "nM"
    bio1['concentration'] = extracted_concentration
    bio1['concentration_unit'] = "uM"
    bio1['source_db'] = "EUbOPEN"
    bio1['source'] = np.nan
    bio1['xref_id'] = np.nan
    
    bio1['assay_type'] = "biochemical"
    bio1['value'] = df['Affinity Biochemical (nM)']
    bio1['bioactivity_type'] = df['Affinity Biochemical Definition']
    bio1['relation'] = df['Affinity Biochemical Relation']
    bio1['assay_description'] = df['Affinity Biochemical Assay Type']
    bio1['source_xref'] = df['Affinity Biochemical Source Knowledge']
    
    # --- Create Cellular DataFrame ---
    bio2 = pd.DataFrame()
    bio2['inchikey'] = df[inchikey_col]
    bio2['target_key'] = df['Target ID']
    bio2['moa'] = df['Mode of Action']
    bio2['cell_line'] = np.nan
    bio2['unit'] = "nM"
    bio2['concentration'] = extracted_concentration
    bio2['concentration_unit'] = "uM"
    bio2['source_db'] = "EUbOPEN"
    bio2['source'] = np.nan
    bio2['xref_id'] = np.nan
    
    bio2['assay_type'] = "cell"
    bio2['value'] = df['Affinity On-target Cellular (nM)']
    bio2['bioactivity_type'] = df['Affinity On-target Cellular Definition']
    bio2['relation'] = df['Affinity On-target Cellular Relation']
    bio2['assay_description'] = df['Affinity On-target Cellular Assay Type']
    bio2['source_xref'] = df['Affinity on-target cellular Source Knowledge']
    
    # Combine both DataFrames
    bioactivity_df = pd.concat([bio1, bio2], ignore_index=True)
    
    # Remove rows where the bioactivity value is missing to prevent unnecessary duplications
    bioactivity_df = bioactivity_df.dropna(subset=['value'])
    bioactivity_df['relation'] = bioactivity_df['relation'].replace(0, np.nan)
    bioactivity_df = normalize_negative_control_values(bioactivity_df)
    bioactivity_df = bioactivity_df.dropna(subset=['value'])
    
    # Save to TSV
    bioactivity_df.to_csv("bioactivity.tsv", sep='\t', index=False)
    print("Successfully generated bioactivity.tsv")
    
    print("All TSV files have been generated successfully!")

if __name__ == "__main__":
    generate_tsvs()