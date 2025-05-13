#!/usr/bin/env python3
"""
STRING Network Generator for PerturbViz

This script generates a gene interaction network file from the STRING database,
suitable for use with the PerturbViz visualization tool.

Usage:
  python generate_string_network.py --genes perturbed_genes.txt --degs deseq_results.tsv --output string_network.tsv
"""

import argparse
import pandas as pd
import requests
import io
import sys
import time

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate network file from STRING database')
    
    parser.add_argument('--genes', '-g', required=True,
                      help='File containing perturbed genes (one per line) or comma-separated list')
    parser.add_argument('--degs', '-d',
                      help='DEG file with gene and log2FC columns (optional)')
    parser.add_argument('--species', '-s', default='9606',
                      help='NCBI taxonomy ID (default: 9606 for human)')
    parser.add_argument('--score', type=int, default=400,
                      help='Minimum interaction confidence score (0-1000, default: 400)')
    parser.add_argument('--output', '-o', default='string_network.tsv',
                      help='Output network file name')
    parser.add_argument('--additional', '-a', type=int, default=50,
                      help='Number of additional interactors to include (default: 50)')
    
    return parser.parse_args()

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

def load_deg_genes(deg_file, top_n=10):
    """Load top differentially expressed genes from DEG file."""
    try:
        df = pd.read_csv(deg_file, sep='\t')
        
        # Ensure proper column names
        if 'gene' not in df.columns or 'log2FC' not in df.columns:
            if len(df.columns) >= 2:
                df.columns = ['gene', 'log2FC'] + list(df.columns[2:])
                print("Warning: Renamed columns to 'gene' and 'log2FC'")
            else:
                print("Error: DEG file must have at least two columns")
                return [], []
        
        # Sort by absolute log2FC
        df['abs_log2FC'] = df['log2FC'].abs()
        sorted_df = df.sort_values('abs_log2FC', ascending=False)
        
        # Get top upregulated genes
        up_genes = sorted_df[sorted_df['log2FC'] > 0].head(top_n)['gene'].tolist()
        
        # Get top downregulated genes
        down_genes = sorted_df[sorted_df['log2FC'] < 0].head(top_n)['gene'].tolist()
        
        print(f"Selected top {len(up_genes)} upregulated and {len(down_genes)} downregulated genes")
        return up_genes, down_genes
    
    except Exception as e:
        print(f"Error loading DEG file: {e}")
        return [], []

def get_string_network(gene_list, species='9606', score_threshold=400, additional_interactors=50):
    """
    Retrieve protein-protein interactions from STRING database.
    
    Args:
        gene_list: List of gene symbols
        species: NCBI taxonomy ID (9606 for human)
        score_threshold: Minimum interaction score (0-1000)
        additional_interactors: Number of additional interacting proteins to include
        
    Returns:
        DataFrame with source and target columns
    """
    print(f"\nQuerying STRING database for {len(gene_list)} genes...")
    print(f"Using confidence score threshold: {score_threshold}")
    
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "network"
    
    # Set parameters
    params = {
        "identifiers": "%0d".join(gene_list),  # Your proteins
        "species": species,                    # Species NCBI identifier
        "caller_identity": "PerturbViz",       # Your app name
        "add_nodes": additional_interactors,   # Additional interactors
        "network_type": "functional",         # Functional associations network
        "required_score": score_threshold      # Minimum interaction score
    }
    
    # Call STRING API
    try:
        print("Sending request to STRING database...")
        response = requests.post(f"{string_api_url}/{output_format}/{method}", data=params)
        
        # Check if the request was successful
        if response.status_code != 200:
            print(f"Error: HTTP Status Code {response.status_code}")
            print(f"Response: {response.text}")
            return None
        
        # Read the result into a DataFrame
        network_data = pd.read_csv(io.StringIO(response.text), sep='\t')
        
        if network_data.empty:
            print("No interactions found for the provided genes.")
            return None
            
        # Print information about the network
        print(f"Retrieved network with {len(network_data)} interactions")
        print(f"Network contains {len(set(network_data['preferredName_A'].tolist() + network_data['preferredName_B'].tolist()))} unique genes")
        
        # Create simplified network format
        simplified_network = network_data[['preferredName_A', 'preferredName_B']].copy()
        simplified_network.columns = ['source', 'target']
        
        return simplified_network
        
    except requests.exceptions.RequestException as e:
        print(f"Error connecting to STRING database: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def main():
    """Main function to run the script."""
    args = parse_arguments()
    
    # Load perturbed genes
    perturbed_genes = load_gene_list(args.genes)
    
    # Load DEGs if provided
    up_genes, down_genes = [], []
    if args.degs:
        up_genes, down_genes = load_deg_genes(args.degs)
    
    # Combine all gene lists, remove duplicates
    all_genes = list(set(perturbed_genes + up_genes + down_genes))
    
    print(f"Total input genes: {len(all_genes)}")
    print(f"  - Perturbed genes: {len(perturbed_genes)}")
    print(f"  - Upregulated genes: {len(up_genes)}")
    print(f"  - Downregulated genes: {len(down_genes)}")
    
    # Get network from STRING
    network_df = get_string_network(all_genes, args.species, args.score, args.additional)
    
    if network_df is None or len(network_df) == 0:
        print("\nFailed to retrieve network or no interactions found.")
        sys.exit(1)
    
    # Save network file
    try:
        network_df.to_csv(args.output, sep='\t', index=False)
        print(f"\nSuccessfully saved network file to {args.output}")
        print(f"Network has {len(network_df)} edges")
    except Exception as e:
        print(f"Error saving network file: {e}")
        sys.exit(1)
    
    print("\nNetwork generation complete! You can now use this file with PerturbViz.")
    print(f"Example: python perturbviz.py --network {args.output} --perturbed {args.genes} --degs {args.degs or 'your_degs_file.tsv'}")

if __name__ == "__main__":
    print("STRING Network Generator for PerturbViz")
    print("======================================")
    main()
