

# """GENERATION OF SEQUENCE - LENGTH AND NUMBER OF SEQUENCE IS GIVEN"""

# # Production of random sequence of length (n) and no of sequence (n)
# def sequence_generation(length, number):
#   Bases = ['A','C','T','G']
#   sequence = np.random.choice(Bases, size = (number, length))
#   return sequence

# """Structure DATAFRAME with Headers"""

# df_0 = pd.DataFrame(sequence_generation(12,10))
# def structure_df(df):
#   df.index.name = 'S.no'
#   df['Sequence'] = [' '.join(row.values.astype(str)) for _, row in df.iterrows()]
#   df = df[['Sequence']]
#   return df

# df_1 = structure_df(df_0)

import numpy as np
import pandas as pd
import os

"""GENERATION OF SEQUENCE - LENGTH AND NUMBER OF SEQUENCE IS GIVEN"""

# Production of random sequence of length (n) and number of sequences (n)
def sequence_generation(length, number):
    Bases = ['A', 'C', 'T', 'G']
    sequences = [''.join(np.random.choice(Bases, length)) for _ in range(number)]
    return sequences

"""Structure DATAFRAME with Headers"""

def structure_df(sequences):
    if isinstance(sequences, pd.DataFrame):
        sequences = sequences[0].tolist()
    df = pd.DataFrame({'Sequence': sequences})
    df.index = [f'Sequence{i + 1}' for i in range(len(sequences))]
    df.index.name = 's.no'
    return df

# Generate sequences and create DataFrame
df_1 = structure_df(sequence_generation(12, 10))

"""Save into Folder"""

def save_df_to_excel(df, output_folder, filename):
    try:
        os.makedirs(output_folder, exist_ok=True)
        filepath = os.path.join(output_folder, filename)
        df.to_csv(filepath, index=True)
        print(f"DataFrame saved successfully to: {filepath}")
    except Exception as e:
        print(f"Error saving DataFrame to Excel: {e}")

