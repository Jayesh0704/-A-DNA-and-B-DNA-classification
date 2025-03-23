# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import os
from itertools import combinations
import logging
import networkx as nx

# Normalization values for nucleotides
normalization_values = {
    'DA': 422.25,  # Adenine
    'DT': 347.54,  # Thymine
    'DC': 352.78,  # Cytosine
    'DG': 492.65   # Guanine
}

def load_file(file_path):
    try:
        data = pd.read_csv(file_path)
        return data
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    required_columns = ['Atom Number', 'Residue Sequence Number', 'Residue Name', 'Chain ID', 'X', 'Y', 'Z']
    missing_columns = [col for col in required_columns if col not in data.columns]
    if missing_columns:
        print(f"Missing required columns: {', '.join(missing_columns)}")
    return data

def preprocess_data(data):
    filtered_data = data[['Atom Number', 'Residue Sequence Number', 'Residue Name', 'Chain ID', 'X', 'Y', 'Z']].copy()
    residue_names = set(filtered_data['Residue Name'])
    normalised_data = set(normalization_values.keys())
    if not residue_names.issubset(normalised_data):
        print(f"Residue names not found in normalization values: {residue_names - normalised_data}")
    residue_to_nucleotide = filtered_data[['Residue Sequence Number', 'Residue Name', 'Chain ID']].drop_duplicates()
    residue_to_nucleotide = residue_to_nucleotide.set_index('Residue Sequence Number').to_dict(orient='index')

    return filtered_data, residue_to_nucleotide

def build_residue_atom_mapping(filtered_data):
    residue_atoms = filtered_data.groupby('Residue Sequence Number').apply(lambda x: x.index.tolist()).to_dict()
    return residue_atoms

def calculate_interaction_strength(residue_i, residue_j, num_contacts, residue_to_nucleotide):
    nucleotide_i = residue_to_nucleotide.get(residue_i, {}).get('Residue Name')
    nucleotide_j = residue_to_nucleotide.get(residue_j, {}).get('Residue Name')
    N_i = normalization_values.get(nucleotide_i, 1)
    N_j = normalization_values.get(nucleotide_j, 1)
    if num_contacts > 0:
        return (num_contacts * 100) / np.sqrt(N_i * N_j)
    return 0

def construct_dsg(filtered_data, residue_to_nucleotide, cutoff_distance):
    logging.info("Constructing DNA Structure Graph (DSG)...")
    residue_atoms = build_residue_atom_mapping(filtered_data)
    G = nx.Graph()
    interaction_details = []
    residues = list(residue_to_nucleotide.keys())
    total_pairs = len(residues) * (len(residues) - 1) // 2
    logging.info(f"Total unique residues: {len(residues)}")
    logging.info(f"Total residue pairs to process: {total_pairs}")
    for idx, (residue_i, residue_j) in enumerate(combinations(residues, 2), 1):
        if idx % 1000 == 0 or idx == 1:
            logging.info(f"Processing residue pair {idx}/{total_pairs}...")

        if residue_i == residue_j:
            continue

        chain_i = residue_to_nucleotide.get(residue_i, {}).get('Chain ID')
        chain_j = residue_to_nucleotide.get(residue_j, {}).get('Chain ID')
        same_chain = chain_i == chain_j

        if same_chain:
            residue_i_seq = residue_i
            residue_j_seq = residue_j
            if abs(residue_i_seq - residue_j_seq) < 2:
                continue

        atoms_i = residue_atoms.get(residue_i, [])
        atoms_j = residue_atoms.get(residue_j, [])

        if not atoms_i or not atoms_j:
            continue

        coords_i = filtered_data.loc[atoms_i, ['X', 'Y', 'Z']].values
        coords_j = filtered_data.loc[atoms_j, ['X', 'Y', 'Z']].values
        distances = cdist(coords_i, coords_j)
        num_contacts = np.sum(distances < cutoff_distance)

        if num_contacts == 0:
            continue

        interaction_strength = calculate_interaction_strength(residue_i, residue_j, num_contacts, residue_to_nucleotide)

        if interaction_strength > 0:
            G.add_node(residue_i)
            G.add_node(residue_j)
            G.add_edge(residue_i, residue_j, weight=interaction_strength)

            interaction_details.append([
                residue_i,
                residue_j,
                residue_to_nucleotide.get(residue_i, {}).get('Residue Name'),
                residue_to_nucleotide.get(residue_j, {}).get('Residue Name'),
                num_contacts,
                interaction_strength
            ])

    return G, interaction_details

def save_interaction_details(interaction_details, output_file_path):
    if not interaction_details:
        logging.warning("No interaction details to save.")
        return

    interaction_df = pd.DataFrame(
        interaction_details,
        columns=[
            'Residue_Number_1',
            'Residue_Number_2',
            'Nucleotide_1',
            'Nucleotide_2',
            'Atomic_Pair_Count',
            'Interaction_Strength'
        ]
    )
    try:
        interaction_df.to_csv(output_file_path, index=False)
        logging.info(f"Interaction details saved to {output_file_path}.")
    except Exception as e:
        logging.error(f"Failed to save interaction details to {output_file_path}: {e}")
        raise

def process_file(input_file_path, output_file_path, cutoff_distance):
    logging.info(f"Starting processing for file: {input_file_path}")
    data = load_file(input_file_path)
    filtered_data, residue_to_nucleotide = preprocess_data(data)
    G, interaction_details = construct_dsg(filtered_data, residue_to_nucleotide, cutoff_distance)
    save_interaction_details(interaction_details, output_file_path)
    logging.info(f"Completed processing for file: {input_file_path}")

def process_csv_files(input_folder, output_folder, cutoff_distance):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            input_file_path = os.path.join(input_folder, filename)
            output_file_path = os.path.join(output_folder, f"interaction_{filename}")
            process_file(input_file_path, output_file_path, cutoff_distance)
