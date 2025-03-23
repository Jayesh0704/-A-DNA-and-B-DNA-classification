import os
import logging
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import plotly.graph_objects as go
import seaborn as sns

class DNAStructureGraph:
    """
    DNA Structure Graph analyzer based on the methodology from
    'Classification of DNA Structure through Structure Network Analysis'
    """
    
    def __init__(self, cutoff_distance=6.5, normalization_values=None):
        """Initialize with customizable parameters."""
        self.cutoff_distance = cutoff_distance
        self.normalization_values = normalization_values or {
            'DA': 422.25, # Adenine
            'DT': 347.54, # Thymine
            'DC': 352.78, # Cytosine
            'DG': 492.65  # Guanine
        }
        
        self.structure = None
        self.atom_data = None
        self.interaction_data = None
        self.graph = None
        self.analysis_results = None
        
        # Configure logging
        logging.basicConfig(level=logging.INFO,
                           format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
    
    def load_pdb(self, pdb_file_path):
        """Load structure from PDB file using direct file parsing."""
        self.logger.info(f"Loading PDB file: {pdb_file_path}")
        try:
            # Read the PDB file content
            with open(pdb_file_path, 'r') as file:
                pdb_content = file.readlines()
            
            # Process the PDB content similar to pdb_to_csv.py
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
                    # Extract atom information
                    atom_info = {
                        'Atom Number': int(atom[6:11].strip()),
                        'Atom': atom[12:16].strip(),
                        'Residue Name': atom[17:20].strip(),
                        'Chain ID': str(cid),
                        'Residue Sequence Number': int(atom[22:26].strip()),
                        'X': float(atom[30:38].strip()),
                        'Y': float(atom[38:46].strip()),
                        'Z': float(atom[46:54].strip()),
                        'Occupancy': float(atom[54:60].strip()),
                        'Temp Factor': float(atom[60:66].strip()),
                        'Element': atom[76:78].strip() if len(atom) > 76 else ''
                    }
                    
                    # Filter out hydrogen atoms as in the original code
                    if 'H' not in atom_info['Atom'] and 'H' not in atom_info['Element']:
                        atom_data.append(atom_info)
                        
            self.atom_data = pd.DataFrame(atom_data)
            self.logger.info(f"Extracted data for {len(atom_data)} atoms")
            return self
            
        except Exception as e:
            self.logger.error(f"Error loading PDB file: {e}")
            raise
    
    def calculate_interactions(self):
        """Calculate interactions between residues based on atomic distances."""
        if self.atom_data is None:
            raise ValueError("No atom data available. Load a PDB file first.")
            
        # Filter and prepare data - exactly as in original code
        filtered_data = self.atom_data[['Atom Number', 'Residue Sequence Number',
                                       'Residue Name', 'Chain ID', 'X', 'Y', 'Z']].copy()
        
        # Create residue to nucleotide mapping - matched to original
        residue_to_nucleotide = (
            filtered_data[['Residue Sequence Number', 'Residue Name', 'Chain ID']]
            .drop_duplicates()
            .set_index('Residue Sequence Number')
            .to_dict('index')
        )
        
        # Group atoms by residue - matched to original
        residue_atoms = filtered_data.groupby('Residue Sequence Number').apply(
            lambda x: x.index.tolist()).to_dict()
            
        # Calculate interactions
        interaction_details = []
        residues = list(residue_to_nucleotide.keys())
        
        for i, residue_i in enumerate(residues):
            for residue_j in residues[i+1:]:
                # Skip residues in the same chain that are sequentially close - CRITICAL PART FROM ORIGINAL
                chain_i = residue_to_nucleotide[residue_i]['Chain ID']
                chain_j = residue_to_nucleotide[residue_j]['Chain ID']
                
                # This is the exact condition from the original code
                if chain_i == chain_j and abs(residue_i - residue_j) < 2:
                    continue
                    
                # Calculate distances and contacts
                atoms_i = residue_atoms.get(residue_i, [])
                atoms_j = residue_atoms.get(residue_j, [])
                
                if not atoms_i or not atoms_j:
                    continue
                    
                coords_i = filtered_data.loc[atoms_i, ['X', 'Y', 'Z']].values
                coords_j = filtered_data.loc[atoms_j, ['X', 'Y', 'Z']].values
                
                # Calculate distances using cdist as in original
                distances = cdist(coords_i, coords_j)
                
                # Count number of contacts (atom pairs within cutoff distance)
                num_contacts = np.sum(distances < self.cutoff_distance)
                
                if num_contacts > 0:
                    # Calculate interaction strength - matches original formula exactly
                    nucleotide_i = residue_to_nucleotide[residue_i]['Residue Name']
                    nucleotide_j = residue_to_nucleotide[residue_j]['Residue Name']
                    N_i = self.normalization_values.get(nucleotide_i, 1)
                    N_j = self.normalization_values.get(nucleotide_j, 1)
                    interaction_strength = (num_contacts * 100) / np.sqrt(N_i * N_j)
                    
                    interaction_details.append([
                        residue_i,
                        residue_j,
                        nucleotide_i,
                        nucleotide_j,
                        num_contacts,
                        interaction_strength
                    ])
                    
        # Create interaction DataFrame with exact column names from original
        self.interaction_data = pd.DataFrame(
            interaction_details,
            columns=['Residue_Number_1', 'Residue_Number_2', 'Nucleotide_1',
                    'Nucleotide_2', 'Atomic_Pair_Count', 'Interaction_Strength']
        )
        
        self.logger.info(f"Calculated {len(interaction_details)} residue interactions")
        
        # Create a graph from interaction data
        self.graph = nx.Graph()
        for _, row in self.interaction_data.iterrows():
            self.graph.add_edge(
                row['Residue_Number_1'],
                row['Residue_Number_2'],
                weight=row['Interaction_Strength']
            )
            
        return self
    
    def save_interactions(self, csv_file_path):
        """Save interaction data to CSV file."""
        if self.interaction_data is None:
            raise ValueError("No interaction data available. Calculate interactions first.")
            
        self.interaction_data.to_csv(csv_file_path, index=False)
        self.logger.info(f"Interaction data saved at: {csv_file_path}")
        return self
    
    def load_interactions(self, csv_file_path):
        """Load interaction data from CSV file."""
        try:
            self.interaction_data = pd.read_csv(csv_file_path)
            
            # Create a graph from interaction data
            self.graph = nx.Graph()
            for _, row in self.interaction_data.iterrows():
                self.graph.add_edge(
                    row['Residue_Number_1'],
                    row['Residue_Number_2'],
                    weight=row['Interaction_Strength']
                )
                
            self.logger.info(f"Loaded interaction data from {csv_file_path}")
            return self
            
        except Exception as e:
            self.logger.error(f"Error loading interaction data: {e}")
            raise
    
    def analyze_graph(self, k_values=None, i_min_values=None):
        """Analyze graph for hubs with varying cutoffs."""
        if self.interaction_data is None or self.graph is None:
            raise ValueError("No interaction data available. Calculate or load interactions first.")
            
        # Use EXACTLY these ranges as in the original code
        if k_values is None:
            k_values = range(1, 10)  # 1 through 9
        if i_min_values is None:
            i_min_values = np.arange(0, 6.5, 0.5)  # 0 to 6 with 0.5 step
            
        self.logger.info("Analyzing graph for hub distribution")
        
        data_records = []
        for i_min in i_min_values:
            # Build filtered graph for this interaction strength
            G = nx.Graph()
            for u, v, d in self.graph.edges(data=True):
                if d['weight'] > i_min:
                    G.add_edge(u, v, weight=d['weight'])
                    
            # Get node degrees
            degrees = dict(G.degree())
            
            for k in k_values:
                # Count hubs (nodes with degree > k)
                hubs = [node for node, degree in degrees.items() if degree > k]
                H = len(hubs)
                
                data_records.append({
                    'Imin': i_min,
                    'k': k,
                    'H': H,
                    'hubs': hubs
                })
                
        self.analysis_results = pd.DataFrame(data_records)
        self.logger.info("Graph analysis completed")
        return self.analysis_results
    
    def generate_2d_hub_plots(self, output_dir, k_values=None, file_prefix=''):
        """Generate 2D plots for H vs Imin for specified k values."""
        if not hasattr(self, 'analysis_results') or self.analysis_results is None:
            self.analyze_graph()
            
        if k_values is None:
            k_values = [2, 3, 4, 5, 6]
            
        for k in k_values:
            df_k = self.analysis_results[self.analysis_results['k'] == k]
            
            plt.figure(figsize=(10, 6))
            sns.lineplot(x='Imin', y='H', data=df_k, marker='o')
            plt.xlabel('Interaction Strength (Imin)')
            plt.ylabel('Number of Hubs (H)')
            plt.title(f'Hub Distribution for k = {k}')
            plt.grid(True)
            
            # Save plot
            plot_file = os.path.join(output_dir, f"{file_prefix}H_vs_Imin_k_{k}.png")
            plt.savefig(plot_file)
            plt.close()
            
            self.logger.info(f"Generated hub plot for k={k} at {plot_file}")
    
    def plot_hub_distribution(self, k=4, output_path=None):
        """Plot hub distribution for a given k value."""
        if not hasattr(self, 'analysis_results') or self.analysis_results is None:
            self.analyze_graph()
            
        # Filter results for specific k value
        df_k = self.analysis_results[self.analysis_results['k'] == k]
        
        plt.figure(figsize=(10, 6))
        plt.plot(df_k['Imin'], df_k['H'], marker='o')
        plt.xlabel('Interaction Strength (Imin)')
        plt.ylabel('Number of Hubs (H)')
        plt.title(f'Hub Distribution for k = {k}')
        plt.grid(True)
        
        if output_path:
            plt.savefig(output_path)
            self.logger.info(f"Hub distribution plot saved at: {output_path}")
            plt.close()
        else:
            plt.show()
            
        return self
    
    def get_hub_data_table(self, k_values=None):
        """
        Get hub data in a format suitable for display in a table.
        Returns a wide-format table with Imin as rows and k values as columns.
        """
        if not hasattr(self, 'analysis_results') or self.analysis_results is None:
            self.analyze_graph()
            
        if k_values is None:
            k_values = sorted(self.analysis_results['k'].unique())
            
        # Filter for specified k values
        df_filtered = self.analysis_results[self.analysis_results['k'].isin(k_values)]
        
        # Pivot to get wide format table with Imin as rows and k values as columns
        hub_table = df_filtered.pivot(index='Imin', columns='k', values='H').reset_index()
        hub_table.columns.name = None  # Remove the columns name
        
        return hub_table
    
    def save_hub_data(self, k=None, output_path=None):
        """
        Save hub vs Imin data to CSV for a specific k value or all k values.
        This provides the tabular data requested by the user.
        """
        if not hasattr(self, 'analysis_results') or self.analysis_results is None:
            self.analyze_graph()
        
        if k is not None:
            # Filter results for specific k value
            df_k = self.analysis_results[self.analysis_results['k'] == k]
            # Extract just the columns we need for the tabular data
            hub_data = df_k[['Imin', 'H']].copy()
            # Rename the columns as requested
            hub_data.columns = ['I MIN', 'no of HUBS']
        else:
            # Get pivot table with Imin as rows and k values as columns
            hub_data = self.get_hub_data_table()
            # Rename the Imin column to I MIN
            hub_data.rename(columns={'Imin': 'I MIN'}, inplace=True)
        
        if output_path:
            hub_data.to_csv(output_path, index=False)
            self.logger.info(f"Hub distribution data saved at: {output_path}")
        
        return hub_data
    
    def visualize_graph(self, interaction_cutoff=0, output_path=None):
        """Visualize the DNA structure graph."""
        if self.interaction_data is None or self.graph is None:
            raise ValueError("No interaction data available. Calculate or load interactions first.")
            
        # Create a filtered graph based on interaction cutoff
        G = nx.Graph()
        for u, v, d in self.graph.edges(data=True):
            if d['weight'] > interaction_cutoff:
                G.add_edge(u, v, weight=d['weight'])
                
        # Create interactive graph visualization with spring layout as in original
        pos = nx.spring_layout(G, weight='weight')
        
        edge_x = []
        edge_y = []
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            
        node_x = [pos[node][0] for node in G.nodes()]
        node_y = [pos[node][1] for node in G.nodes()]
        
        # Create edge trace
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            mode='lines',
            hoverinfo='none'
        )
        
        # Create node trace - matches original implementation
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            marker=dict(
                showscale=True,
                colorscale='YlGnBu',
                size=10,
                color=[len(list(G.neighbors(node))) for node in G.nodes()],
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                )
            )
        )
        
        # Create figure with layout from original code
        fig = go.Figure(
            data=[edge_trace, node_trace],
            layout=go.Layout(
                title='DNA Structure Graph',
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(showgrid=False, zeroline=False),
                yaxis=dict(showgrid=False, zeroline=False)
            )
        )
        
        if output_path:
            fig.write_html(output_path)
            self.logger.info(f"Graph visualization saved at: {output_path}")
        else:
            fig.show()
            
        return self
