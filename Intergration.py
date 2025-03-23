import pandas as pd
from openpyxl import workbook
from random_sequence_generator import sequence_generation, structure_df, save_df_to_excel
from pdb_to_csv import process_pdb_files
from code_for_interaction import process_csv_files
from graphs import processed_graphs
from web_automation import dna_sequence_processor

# Generate a random sequence
# data_0 = pd.DataFrame(sequence_generation(12, 50000))  # Adjust arguments as needed
# data_1 = structure_df(data_0)
# file_name = 'random_sequence.csv'
# output_folder = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\input'
# save_df_to_excel(data_1, output_folder, file_name)
# input_folder_A_DNA_PDB = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\input\A_DNA'
# input_folder_B_DNA_PDB = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\input\B_DNA'
# file_path=r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\input\random_sequence.csv'

# dna_sequence_processor(file_path,input_folder_A_DNA_PDB,input_folder_B_DNA_PDB)

# Process PDB files into CSV format
input_folder_A_DNA = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\input\A_DNA'
# input_folder_B_DNA = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\input\B_DNA'


output_folder_A_DNA = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\A_DNA_CSV'
# output_folder_B_DNA = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\B_DNA_CSV'

process_pdb_files(input_folder_A_DNA, output_folder_A_DNA)
# process_pdb_files(input_folder_B_DNA, output_folder_B_DNA)



# Process interaction CSV files with a specified cutoff distance
cutoff_distance_value = 7.0  # Example: Adjust the cutoff distance as needed
input_folder_A_DNA_interaction = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\A_DNA_CSV'
# input_folder_B_DNA_interaction = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\B_DNA_CSV'
output_folder_A_DNA_interaction = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\A_DNA_Interactions'
# output_folder_B_DNA_interaction = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\B_DNA_Interactions'

process_csv_files(input_folder_A_DNA_interaction, output_folder_A_DNA_interaction, cutoff_distance_value)
# process_csv_files(input_folder_B_DNA_interaction, output_folder_B_DNA_interaction, cutoff_distance_value)

# Generate graphs for interaction data
input_folder_A_DNA_GRAPH = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\A_DNA_Interactions'
# input_folder_B_DNA_GRAPH = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\B_DNA_Interactions'
output_folder_A_DNA_GRAPH = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\Graphs\A_DNA_Graphs'
# output_folder_B_DNA_GRAPH = r'D:\SOP LOP DOP\Genomic hotspots\Automation\BIOLOP\output\Graphs\B_DNA_Graphs'

# Process graphs for A-DNA
processed_graphs(input_folder_A_DNA_GRAPH, output_folder_A_DNA_GRAPH)

# Process graphs for B-DNA
# processed_graphs(input_folder_B_DNA_GRAPH, output_folder_B_DNA_GRAPH)

print(f"Graph generation completed.")
print(f"Graphs for A-DNA are saved in {output_folder_A_DNA_GRAPH}.")
# print(f"Graphs for B-DNA are saved in {output_folder_B_DNA_GRAPH}.")