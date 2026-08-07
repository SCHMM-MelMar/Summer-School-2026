import numpy as np
import pandas as pd
import requests
from rdkit import Chem

# ==========================================
# Helper Functions
# ==========================================

def normalize_negative_control_values(bioactivity_df):
    result = bioactivity_df.copy()
    negative_control_mask = (
        result["value"].fillna("").astype(str).str.strip().str.lower() == "negative control"
    )

    if negative_control_mask.any():
        result.loc[negative_control_mask, "moa"] = "Negative control"
        result.loc[negative_control_mask, "value"] = np.nan

    return result

def smiles_to_inchikey(smiles):
    """Generates an InChIKey from a SMILES string using RDKit."""
    if pd.isna(smiles) or not isinstance(smiles, str) or not smiles.strip():
        return np.nan
        
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol:
            return Chem.MolToInchiKey(mol)
    except Exception:
        pass 
    return np.nan

def get_chembl_id_from_inchikey(inchikey):
    """Fetches the ChEMBL ID from the ChEMBL API using an InChIKey."""
    if pd.isna(inchikey):
        return np.nan
        
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key={inchikey}"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if data.get('page_meta', {}).get('total_count', 0) > 0:
                return data['molecules'][0]['molecule_chembl_id']
    except Exception as e:
        print(f"  [API Error] for InChIKey {inchikey}: {e}")
        
    return np.nan

# ==========================================
# Main Processing Function
# ==========================================

def generate_tsvs():
    input_file = "EUbOPEN compounds_cv_approved.csv"
    print(f"Reading data from {input_file}...")
    
    # Read the CSV file once
    df = pd.read_csv(input_file, low_memory=False)
    
    # Handle slight variations in user's column naming (InChi Key vs InChiKey)
    inchikey_col = 'Compound InChi Key' if 'Compound InChi Key' in df.columns else 'Compound InChiKey'
    
    # Pre-clean InChIKeys
    df[inchikey_col] = df[inchikey_col].astype(str).str.strip()
    df[inchikey_col] = df[inchikey_col].replace(['', 'nan', 'None'], np.nan)
    
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
    
    print("Processing missing InChIKeys using RDKit...")
    compound_df['inchikey'] = compound_df.apply(
        lambda row: smiles_to_inchikey(row['smiles']) if pd.isna(row['inchikey']) else row['inchikey'], 
        axis=1
    )
    
    # Drop rows without a valid InChIKey, then deduplicate
    compound_df = compound_df.dropna(subset=['inchikey'])
    compound_df = compound_df.drop_duplicates(subset=['inchikey'])
    
    print(f"Fetching ChEMBL IDs for {len(compound_df)} unique compounds (this may take a moment)...")
    compound_df['chembl_id'] = compound_df['inchikey'].apply(get_chembl_id_from_inchikey)
    
    # Reorder columns to match the target schema
    compound_df = compound_df[['inchikey', 'smiles', 'chembl_id', 'name']]
    compound_df.to_csv("compound.tsv", sep='\t', index=False)
    print("Successfully generated compound.tsv")

    # ==========================================
    # 4. bioactivity.tsv
    # ==========================================
    # Update the raw dataframe's InChIKey column so the bioactivity table 
    # benefits from the missing InChIKeys we just generated from SMILES.
    inchikey_map = dict(zip(compound_df['smiles'], compound_df['inchikey']))
    df[inchikey_col] = df.apply(
        lambda row: inchikey_map.get(row['Compound SMILES'], row[inchikey_col]), 
        axis=1
    )

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
    bioactivity_df['relation'] = bioactivity_df['relation'].replace('0', np.nan)
    
    # Normalize 'Negative Control' text before coercing to numeric
    bioactivity_df = normalize_negative_control_values(bioactivity_df)

    # add missing InChiKeys
    print("Processing missing InChIKeys using RDKit...")
    bioactivity_df['inchikey'] = bioactivity_df.apply(
        lambda row: smiles_to_inchikey(row['smiles']) if pd.isna(row['inchikey']) else row['inchikey'], 
        axis=1
    )
    
    # Force value to numeric
    bioactivity_df['value'] = pd.to_numeric(bioactivity_df['value'], errors='coerce')
    
    bioactivity_df = bioactivity_df.dropna(subset=['value'])
    
    # Save to TSV
    bioactivity_df.to_csv("bioactivity.tsv", sep='\t', index=False)
    print("Successfully generated bioactivity.tsv")
    
    print("All TSV files have been generated successfully!")

if __name__ == "__main__":
    generate_tsvs()