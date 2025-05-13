#!/usr/bin/env python3
"""
PerturbViz - 
Visualize perturbed genes in a network context with DESeq2 results



Created by: Gabriela Grant
Date: May 3, 2025
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
import os
import sys
import matplotlib.patheffects as path_effects

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Visualize perturbed genes in a network context')
    
    parser.add_argument('--network', '-n', required=True,
                      help='Network file with source and target columns')
    parser.add_argument('--perturbed', '-p', required=True,
                      help='File containing perturbed genes (one per line) or comma-separated list')
    parser.add_argument('--degs', '-d',
                      help='DESeq2 results file')
    parser.add_argument('--output', '-o', default='perturb_viz',
                      help='Output file prefix')
    parser.add_argument('--layout', default='spring',
                      choices=['spring', 'kamada_kawai', 'circular', 'spectral'],
                      help='Network layout algorithm')
    parser.add_argument('--node_size', type=int, default=100,
                      help='Base node size')
    parser.add_argument('--width', type=float, default=12,
                      help='Figure width in inches')
    parser.add_argument('--height', type=float, default=10,
                      help='Figure height in inches')
    parser.add_argument('--dpi', type=int, default=300,
                      help='Figure DPI')
    parser.add_argument('--label_offset', type=float, default=0.05,
                      help='Offset distance for labels')
    
    return parser.parse_args()

def load_network(network_file):
    """Load network from file with error handling."""
    try:
        print(f"Loading network from {network_file}")
        df = pd.read_csv(network_file, sep=None, engine='python')
        
        # Make sure we have source and target columns
        if 'source' not in df.columns or 'target' not in df.columns:
            # Try to use first two columns
            if len(df.columns) >= 2:
                df = df.rename(columns={df.columns[0]: 'source', df.columns[1]: 'target'})
                print("Network file didn't have 'source' and 'target' columns. Using first two columns.")
            else:
                raise ValueError("Network file must have at least two columns")
        
        # Create networkx graph
        G = nx.from_pandas_edgelist(df, 'source', 'target')
        
        print(f"Network has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        return G
    
    except Exception as e:
        print(f"Error loading network file: {str(e)}")
        sys.exit(1)

def load_gene_list(gene_input):
    """Load genes from file or parse comma-separated list."""
    if ',' in gene_input:
        # Parse comma-separated list
        genes = [g.strip() for g in gene_input.split(',')]
        print(f"Parsed {len(genes)} genes from input string")
        return genes
    
    try:
        # Read from file (one gene per line)
        with open(gene_input, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(genes)} genes from file {gene_input}")
        return genes
    except Exception as e:
        print(f"Error reading gene file: {e}")
        sys.exit(1)

def load_degs(degs_file):
    """Load DEG results from a DESeq2 output file."""
    try:
        print(f"Loading DEG results from {degs_file}")
        
        # Try different separators and formats
        try:
            # First try with comma as separator
            degs_df = pd.read_csv(degs_file)
        except:
            # Then try with tab as separator
            try:
                degs_df = pd.read_csv(degs_file, sep='\t')
            except:
                # Last try with first column as index (common DESeq2 format)
                degs_df = pd.read_csv(degs_file, sep=None, engine='python', index_col=0)
                # Add gene column from index
                degs_df = degs_df.reset_index()
        
        # Print columns for debugging
        print(f"DEG file columns: {degs_df.columns.tolist()}")
        
        # Handle typical DESeq2 output column names
        column_mapping = {
            'index': 'gene',
            'log2FoldChange': 'log2FC',
            'padj': 'padj',
            'pvalue': 'pvalue'
        }
        
        # Rename columns that exist
        for old_name, new_name in column_mapping.items():
            if old_name in degs_df.columns:
                degs_df = degs_df.rename(columns={old_name: new_name})
        
        # Check if we need to create a gene column
        if 'gene' not in degs_df.columns:
            # If no gene column, use the first column as gene identifiers
            gene_col = degs_df.columns[0]
            print(f"No 'gene' column found. Using '{gene_col}' column as gene identifiers.")
            degs_df = degs_df.rename(columns={gene_col: 'gene'})
        
        # Check if we have log2FC information
        if 'log2FC' not in degs_df.columns:
            # Look for common alternatives
            for col in degs_df.columns:
                if 'log2' in col.lower() and 'fold' in col.lower():
                    print(f"Found log2 fold change column: {col}")
                    degs_df['log2FC'] = degs_df[col]
                    break
        
        # If we still don't have log2FC, warn the user
        if 'log2FC' not in degs_df.columns:
            print("Warning: No log2 fold change column found in DEG file")
            
        print(f"Loaded {len(degs_df)} differentially expressed genes")
        return degs_df
    
    except Exception as e:
        print(f"Error loading DEG file: {str(e)}")
        print("DEG file loading failed. Creating minimal DEG dataframe.")
        return pd.DataFrame(columns=['gene', 'log2FC'])

def categorize_genes(perturbed_genes, degs_df, G):
    """Categorize genes based on perturbation and differential expression."""
    gene_categories = {}
    
    # Initialize all nodes as 'other'
    for node in G.nodes():
        gene_categories[node] = 'other'
    
    # Mark perturbed genes
    for gene in perturbed_genes:
        if gene in G.nodes():
            gene_categories[gene] = 'perturbed'
    
    # If we have DEG data, mark up- and down-regulated genes
    if degs_df is not None and 'log2FC' in degs_df.columns:
        # Create a dictionary for quick lookups
        gene_to_fc = dict(zip(degs_df['gene'], degs_df['log2FC']))
        
        # Find significantly up/down genes that are in the network
        for gene in G.nodes():
            if gene in gene_to_fc and gene not in perturbed_genes:
                fc = gene_to_fc[gene]
                # Skip if fold change is NaN
                if pd.isna(fc):
                    continue
                    
                # Classify as up or down regulated
                if fc > 1:  # Log2FC > 1 means upregulated
                    gene_categories[gene] = 'up'
                elif fc < -1:  # Log2FC < -1 means downregulated
                    gene_categories[gene] = 'down'
    
    # Count genes in each category
    category_counts = {cat: list(gene_categories.values()).count(cat) for cat in ['perturbed', 'up', 'down', 'other']}
    print(f"Gene categories: {category_counts}")
    
    return gene_categories

def create_gene_table(G, gene_categories, fold_changes, output_prefix):
    """Create a table of important genes with their network properties and expression data."""
    try:
        # Calculate network centrality metrics
        degree_centrality = nx.degree_centrality(G)
        betweenness_centrality = nx.betweenness_centrality(G)
        
        # Create data for table
        table_data = []
        
        # Get important genes
        important_genes = [gene for gene, cat in gene_categories.items() 
                        if cat in ('perturbed', 'up', 'down')]
        
        for gene in important_genes:
            row = {
                'Gene': gene,
                'Category': gene_categories[gene],
                'Degree': G.degree(gene),
                'Centrality': degree_centrality.get(gene, 0),
                'Betweenness': betweenness_centrality.get(gene, 0),
                'log2FC': fold_changes.get(gene, 'NA')
            }
            table_data.append(row)
        
        # Convert to DataFrame
        table_df = pd.DataFrame(table_data)
        
        # Sort by category and then by fold change or centrality
        if 'log2FC' in table_df.columns and not all(table_df['log2FC'] == 'NA'):
            # Convert log2FC to numeric where possible
            table_df['log2FC_num'] = pd.to_numeric(table_df['log2FC'], errors='coerce')
            table_df = table_df.sort_values(['Category', 'log2FC_num'], ascending=[True, False])
            table_df = table_df.drop(columns=['log2FC_num'])
        else:
            table_df = table_df.sort_values(['Category', 'Centrality'], ascending=[True, False])
        
        # Save to CSV
        table_file = f"{output_prefix}_important_genes.csv"
        table_df.to_csv(table_file, index=False)
        print(f"Saved important gene table to {table_file}")
    
    except Exception as e:
        print(f"Warning: Could not create gene table: {e}")

def create_labeled_network(G, pos, gene_categories, fold_changes, output_prefix, args):
    """Create a separate network visualization with clear labels for all important genes."""
    # Create figure with explicit axes (OPTION 1 FIX)
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    ax.set_title("Gene Expression Network - Detailed View with Labels", fontsize=14)
    
    # Get important nodes
    important_nodes = []
    for category in ['perturbed', 'up', 'down']:
        important_nodes.extend([node for node, cat in gene_categories.items() if cat == category])
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.2, edge_color='gray', width=0.6)
    
    # Draw nodes by category
    node_artists = {}
    for category, color in [('other', 'lightgray'), ('down', 'blue'), ('up', 'red'), ('perturbed', 'purple')]:
        nodes = [node for node, cat in gene_categories.items() if cat == category]
        if nodes:
            size = args.node_size/3 if category == 'other' else args.node_size/1.5 if category in ['up', 'down'] else args.node_size
            alpha = 0.4 if category == 'other' else 0.6 if category in ['up', 'down'] else 0.8
            nodes_drawn = nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=nodes, 
                             node_color=color, node_size=size, alpha=alpha)
            node_artists[category] = nodes_drawn
    
    # Draw offset labels with connecting lines for important nodes
    offset = args.label_offset * 1.5  # Use larger offset for the detailed view
    
    for node in important_nodes:
        if node not in pos:
            continue
            
        # Get node position
        node_x, node_y = pos[node]
        
        # Calculate offset position for label (more distributed to avoid overcrowding)
        angle = np.random.uniform(0, 2*np.pi)
        label_x = node_x + offset * np.cos(angle)
        label_y = node_y + offset * np.sin(angle)
        
        # Draw connecting line from node to label
        ax.plot([node_x, label_x], [node_y, label_y], 'k-', lw=0.5, alpha=0.6, zorder=1)
        
        # Add gene name with fold change if available
        label_text = node
        if node in fold_changes:
            fc = fold_changes[node]
            direction = "↑" if fc > 0 else "↓"
            label_text = f"{node}\n{direction}{abs(fc):.2f}"
            
        text = ax.text(label_x, label_y, label_text, fontsize=9, fontweight='bold',
                     ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8,
                                                       edgecolor='gray', boxstyle='round,pad=0.2'),
                     zorder=2)
        
        # Add path effects for better visibility
        text.set_path_effects([
            path_effects.Stroke(linewidth=2, foreground='white'),
            path_effects.Normal()
        ])
    
    # Create colorbar for perturbed genes if they have fold change data
    if 'perturbed' in node_artists and fold_changes:
        try:
            perturbed_nodes = [node for node, cat in gene_categories.items() if cat == 'perturbed' and node in fold_changes]
            if perturbed_nodes:
                # Get fold changes for perturbed genes
                perturbed_fc = [fold_changes[gene] for gene in perturbed_nodes if gene in fold_changes]
                
                if perturbed_fc:
                    # Create a diverging colormap for fold changes
                    cmap = plt.cm.RdBu_r
                    vmin = min(min(perturbed_fc), -1)
                    vmax = max(max(perturbed_fc), 1)
                    
                    # Create a ScalarMappable for the colorbar
                    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
                    sm.set_array([])
                    
                    # Add colorbar with explicit axes reference (OPTION 1 FIX)
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="2%", pad=0.5)
                    plt.colorbar(sm, cax=cax, label='log2 Fold Change')
        except Exception as e:
            print(f"Warning: Error creating detailed view colorbar: {e}")
    
    ax.axis('off')
    plt.tight_layout()
    
    # Save detailed visualization
    detailed_output_file = f"{output_prefix}_network_labeled.png"
    plt.savefig(detailed_output_file, dpi=args.dpi, bbox_inches='tight')
    print(f"Saved detailed labeled visualization to {detailed_output_file}")
    
    # Create a table of important genes and their attributes
    create_gene_table(G, gene_categories, fold_changes, output_prefix)

def visualize_network(G, gene_categories, degs_df=None, output_prefix='perturb_network', args=None):
    """
    Visualize network with highlighted gene categories and gene expression data.
    
    Args:
        G (nx.Graph): Network to visualize
        gene_categories (dict): Gene category assignments
        degs_df (pd.DataFrame): DEG results with gene and log2FC columns
        output_prefix (str): Output file prefix
        args (Namespace): Command-line arguments with visualization settings
    """
    if args is None:
        # Default settings if not provided
        class DefaultArgs:
            pass
        args = DefaultArgs()
        args.layout = 'spring'
        args.node_size = 80
        args.width = 14
        args.height = 12
        args.dpi = 300
        args.label_offset = 0.05
    
    print("\nCreating enhanced network visualization...")
    
    # Create figure with explicitly defined axes (OPTION 1 FIX)
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    
    # Generate layout
    print(f"Calculating {args.layout} layout...")
    if args.layout == 'spring':
        pos = nx.spring_layout(G, k=0.3, iterations=50, seed=42)
    elif args.layout == 'kamada_kawai':
        try:
            pos = nx.kamada_kawai_layout(G)
        except:
            print("Kamada-Kawai layout failed, falling back to spring layout")
            pos = nx.spring_layout(G, k=0.3, iterations=50, seed=42)
    elif args.layout == 'circular':
        pos = nx.circular_layout(G)
    elif args.layout == 'spectral':
        try:
            pos = nx.spectral_layout(G)
        except:
            print("Spectral layout failed, falling back to spring layout")
            pos = nx.spring_layout(G, k=0.3, iterations=50, seed=42)
    
    # Create dictionaries to store node lists and fold change data
    node_categories = {
        'perturbed': [],
        'up': [],
        'down': [],
        'other': []
    }
    
    # Store fold changes for color mapping if DEGs data is available
    fold_changes = {}
    if degs_df is not None and 'log2FC' in degs_df.columns:
        try:
            # Convert DEGs dataframe to dictionary for quick lookup
            fc_dict = dict(zip(degs_df['gene'], degs_df['log2FC']))
            
            # Get fold change for each node if available
            for node in G.nodes():
                if node in fc_dict:
                    try:
                        fc_value = fc_dict[node]
                        if not pd.isna(fc_value):  # Skip NaN values
                            fold_changes[node] = fc_value
                    except:
                        pass  # Skip if there's an issue with this node
        except Exception as e:
            print(f"Warning: Error processing fold changes: {e}")
    
    # Organize nodes by category
    for node, category in gene_categories.items():
        node_categories[category].append(node)
    
    # Draw edges first (in the background)
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, edge_color='gray', width=0.8)
    
    # Draw regular nodes (background)
    if node_categories['other']:
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=node_categories['other'],
                             node_color='lightgray', node_size=args.node_size/2, alpha=0.5)
    
    # Draw downregulated nodes
    if node_categories['down']:
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=node_categories['down'],
                             node_color='blue', node_size=args.node_size, alpha=0.8)
    
    # Draw upregulated nodes
    if node_categories['up']:
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=node_categories['up'],
                             node_color='red', node_size=args.node_size, alpha=0.8)
    
    # Create a colormap for perturbed genes if fold change data is available
    perturbed_nodes = None
    if node_categories['perturbed'] and fold_changes:
        try:
            # Get fold changes for perturbed genes
            perturbed_fc = []
            valid_perturbed = []
            for g in node_categories['perturbed']:
                if g in fold_changes:
                    perturbed_fc.append(fold_changes[g])
                    valid_perturbed.append(g)
            
            if perturbed_fc:
                # Create a diverging colormap for fold changes
                cmap = plt.cm.RdBu_r
                vmin = min(min(perturbed_fc), -1)  # At least -1 for color range
                vmax = max(max(perturbed_fc), 1)   # At least 1 for color range
                
                # Draw perturbed nodes with color based on fold change
                perturbed_nodes = nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=valid_perturbed,
                                     node_color=perturbed_fc, cmap=cmap, vmin=vmin, vmax=vmax,
                                     node_size=args.node_size*1.5, alpha=1.0, edgecolors='black', linewidths=1)
                
                # Add colorbar for fold change with explicit axes reference (OPTION 1 FIX)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="2%", pad=0.5)
                plt.colorbar(perturbed_nodes, cax=cax, label='log2 Fold Change')
            else:
                # If no fold change data for perturbed genes, draw them in purple
                nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=node_categories['perturbed'],
                                     node_color='purple', node_size=args.node_size*1.5, alpha=1.0,
                                     edgecolors='black', linewidths=1)
        except Exception as e:
            print(f"Warning: Error creating colormap: {e}")
            # Fallback to simple visualization
            nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=node_categories['perturbed'],
                                 node_color='purple', node_size=args.node_size*1.5, alpha=1.0,
                                 edgecolors='black', linewidths=1)
    else:
        # If no fold change data, just draw perturbed nodes in purple
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=node_categories['perturbed'],
                             node_color='purple', node_size=args.node_size*1.5, alpha=1.0,
                             edgecolors='black', linewidths=1)
    
    # Add labels for all important nodes using offset labels with connecting lines
    important_nodes = node_categories['perturbed'] + node_categories['up'] + node_categories['down']
    
    # For label positioning
    offset = args.label_offset  # Offset distance for labels from nodes
    
    # Determine which nodes to label
    if len(G.nodes()) <= 50:
        # For small networks, label all nodes
        nodes_to_label = list(G.nodes())
    else:
        # For larger networks, only label important nodes and top fold change genes
        nodes_to_label = important_nodes.copy()
        
        # Add additional labels for genes with highest fold changes if DEGs data is available
        if degs_df is not None and 'log2FC' in degs_df.columns:
            try:
                # Find top 10 genes with highest absolute fold changes that are in the network
                genes_in_network = set(G.nodes())
                degs_df['abs_log2FC'] = degs_df['log2FC'].abs()
                top_fc_genes = degs_df.sort_values('abs_log2FC', ascending=False)
                top_fc_genes = top_fc_genes[top_fc_genes['gene'].isin(genes_in_network)].head(10)
                
                # Add these genes to nodes_to_label if not already included
                for gene in top_fc_genes['gene']:
                    if gene not in nodes_to_label:
                        nodes_to_label.append(gene)
            except:
                print("Warning: Could not identify top fold change genes")
    
    # Draw offset labels with connecting lines to nodes
    for node in nodes_to_label:
        # Skip if node not in position dictionary
        if node not in pos:
            continue
            
        # Get node position
        node_x, node_y = pos[node]
        
        # Calculate offset position for label with a bit of randomness to avoid overlaps
        angle = np.random.uniform(0, 2*np.pi)  # Random angle
        label_x = node_x + offset * np.cos(angle)
        label_y = node_y + offset * np.sin(angle)
        
        # Draw connecting line from node to label
        ax.plot([node_x, label_x], [node_y, label_y], 'k-', lw=0.5, alpha=0.6, zorder=1)
        
        # Add gene name with white outline for better visibility
        label_text = node
        if node in fold_changes:
            fc = fold_changes[node]
            direction = "↑" if fc > 0 else "↓"
            label_text = f"{node}\n{direction}{abs(fc):.2f}"
            
        text = ax.text(label_x, label_y, label_text, fontsize=8, fontweight='bold',
                     ha='center', va='center', bbox=dict(facecolor='white', alpha=0.7, 
                                                        edgecolor='gray', boxstyle='round,pad=0.3'),
                     zorder=2)
        
        # Add path effects for better visibility
        text.set_path_effects([
            path_effects.Stroke(linewidth=2, foreground='white'),
            path_effects.Normal()
        ])
    
    # Create legend elements
    legend_elements = []
    
    if node_categories['perturbed']:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple',
                                       markersize=10, label=f'Perturbed Genes ({len(node_categories["perturbed"])})'))
    
    if node_categories['up']:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
                                       markersize=10, label=f'Upregulated Genes ({len(node_categories["up"])})'))
    
    if node_categories['down']:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
                                       markersize=10, label=f'Downregulated Genes ({len(node_categories["down"])})'))
    
    if node_categories['other']:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgray',
                                       markersize=10, label=f'Other Genes ({len(node_categories["other"])})'))
    
    # Add legend
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    # Add title with more information
    title_text = f"GEMDiff Perturbed Gene Network Visualization\n"
    
    # Add specific gene counts
    gene_counts = []
    if node_categories['perturbed']:
        gene_counts.append(f"{len(node_categories['perturbed'])} perturbed")
    if node_categories['up']:
        gene_counts.append(f"{len(node_categories['up'])} upregulated")
    if node_categories['down']:
        gene_counts.append(f"{len(node_categories['down'])} downregulated")
    if node_categories['other']:
        gene_counts.append(f"{len(node_categories['other'])} other")
    
    title_text += f"({', '.join(gene_counts)} genes, {G.number_of_edges()} interactions)"
    ax.set_title(title_text, fontsize=14)
    
    ax.axis('off')
    plt.tight_layout()
    
    # Save visualization
    output_file = f"{output_prefix}_network.png"
    plt.savefig(output_file, dpi=args.dpi, bbox_inches='tight')
    print(f"Saved network visualization to {output_file}")
    
    # Create labeled version with gene names
    try:
        create_labeled_network(G, pos, gene_categories, fold_changes, output_prefix, args)
    except Exception as e:
        print(f"Warning: Could not create labeled network: {e}")
    
    plt.close('all')
        
def main():
    """Main function to run the script."""
    args = parse_arguments()
    
    # Load network
    G = load_network(args.network)
    
    # Load perturbed genes
    perturbed_genes = load_gene_list(args.perturbed)
    
    # Load DEGs if provided
    degs_df = None
    if args.degs:
        degs_df = load_degs(args.degs)
    
    # Categorize genes
    gene_categories = categorize_genes(perturbed_genes, degs_df, G)
    
    # Create visualizations
    visualize_network(G, gene_categories, degs_df, args.output, args)
    
    print("\nVisualization complete!")

if __name__ == "__main__":
    print("PerturbViz - Network Visualization for Perturbed Genes")
    print("=======================================================")
    main()
