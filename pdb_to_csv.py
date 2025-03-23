# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd


def process_pdb(pdb_content):
    atom_data = []
    chain_data = []
    current_chain = []
    chain_id = 1

    for line in pdb_content:
        if line.startswith("ATOM"):
            current_chain.append(line)
        elif line.startswith("TER") and current_chain:
            chain_data.append((current_chain, chain_id))
            current_chain = []
            chain_id += 1

    if current_chain:
        chain_data.append((current_chain, chain_id))

    for chain_atoms, cid in chain_data:
        for atom in chain_atoms:
            atom_info = {
                'Atom': 'ATOM',
                'Atom Number': atom[6:11].strip(),
                'Atom Name': atom[12:16].strip(),
                'Residue Name': atom[17:20].strip(),
                'Chain ID': str(cid),
                'Residue Sequence Number': int(atom[22:26].strip()),
                'X': float(atom[30:38].strip()),
                'Y': float(atom[38:46].strip()),
                'Z': float(atom[46:54].strip()),
                'Occupancy': float(atom[54:60].strip()),
                'Temp Factor': float(atom[60:66].strip()),
                'Atom Type': atom[76:].strip()
            }
            atom_data.append(atom_info)

    return pd.DataFrame(atom_data)

def process_pdb_files(input_folder, output_folder):
  for filename in os.listdir(input_folder):
    if filename.endswith('.pdb'):
        with open(os.path.join(input_folder, filename), 'r') as file:
            pdb_content = file.readlines()

        df_pdb = process_pdb(pdb_content)

        df_pdb = df_pdb[~df_pdb['Atom Name'].str.contains('H', na=False) & ~df_pdb['Atom Type'].str.contains('H', na=False)]

        output_file_path = os.path.join(output_folder, filename.replace('.pdb', '.csv'))
        df_pdb.to_csv(output_file_path, index=False)
        print(f"Processed and saved: {output_file_path}")
  print('All files processed.')

