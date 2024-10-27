import os
import argparse
import yaml
import torch
import pandas as pd
from torch_geometric.data import HeteroData
from torch_geometric.transforms import ToUndirected
from tqdm import tqdm
import logging
from pathlib import Path


args = argparse.ArgumentParser()
args.add_argument('--database_files_path', type=str, required=True)
args.add_argument('--config_path', type=str, required=True)
args.add_argument('--output_path', type=str, required=True)
args = args.parse_args()


class GraphLoader:
    def __init__(self, database_files_path, config):
        self.database_files_path = database_files_path
        self.config = config

    def load_node_csvs(self, **kwargs):
        node_dict = {}
        print(f'Loading node CSVs:')
        for node_type in self.config['node_types']:
            node_file_path = f'{self.database_files_path}/{node_type}.csv'
            print(f'\nLoading {node_type} nodes from {node_file_path}')
            try:
                node_df = pd.read_csv(node_file_path, **kwargs)
            except FileNotFoundError:
                print(f"File not found: {node_file_path}")
                break
            except pd.errors.EmptyDataError:
                print(f"File is empty: {node_file_path}")
                break

            all_node_count = len(node_df)
            node_ids = node_df.iloc[:, 0].drop_duplicates()
            unique_node_count = len(node_ids)

            # first column should be the node id
            mapping = {index: i for i, index in enumerate(node_ids)}
            
            node_dict[node_type] = {}
            node_dict[node_type]['mapping'] = mapping
            node_dict[node_type]['num_nodes'] = unique_node_count
            print(f'{unique_node_count} unique nodes from {all_node_count} rows loaded')
        
        # merge and reindex GO nodes, since the edge files contain GO nodes from all three categories
        go_nodes = {**node_dict.get('Molecular_function', {}).get('mapping', {}),
                    **node_dict.get('Biological_process', {}).get('mapping', {}),
                    **node_dict.get('Cellular_component', {}).get('mapping', {})}
        go_mapping = {index: i for i, index in enumerate(go_nodes.keys())}
        node_dict['GO'] = {'mapping': go_mapping}        

        return node_dict
    
    def load_edge_csvs(self, **kwargs):
        def get_node_id_mapping(node_dict, node_type):
            """Helper function to retrieve node ID mapping with potential alternative keys."""
            if node_type not in node_dict:
                print(f"Warning: Node type '{node_type}' not found in node dictionary.")
                return {}, {}
            base_mapping = node_dict[node_type]['mapping']
            try:
                alt_mapping = {k.split(':')[1]: v for k, v in base_mapping.items() if ':' in k}
            except Exception as e:
                raise Exception(f"Error creating alternative mapping for {node_type}: {str(e)}")
            return base_mapping, alt_mapping

        def map_edge_ids(edge_df, column, node_mapping, alt_mapping, type_normalize=False):
            """Maps edge IDs from the CSV to node IDs with primary and alternative mappings."""
            ids = edge_df[column].astype(str).str.lower() if type_normalize else edge_df[column].astype(str)
            mapped_ids = [node_mapping.get(x) or alt_mapping.get(x) or None for x in ids]
            unmapped_ids = {x for x in ids if x not in node_mapping and x not in alt_mapping}

            return [int(x) if x is not None else None for x in mapped_ids], unmapped_ids
        
        edge_dict = {}
        unmatched_node_dict = {}

        for edge in self.config['edge_types']:
            if len(edge) != 5:
                print(f"Warning: Invalid edge configuration: {edge}")
                break
            
            edge_type, src_type, dst_type, src_column, dst_column = edge
            src_node_dict, alt_src_node_dict = get_node_id_mapping(node_dict, src_type)
            dst_node_dict, alt_dst_node_dict = get_node_id_mapping(node_dict, dst_type)

            edge_file_path = f'{self.database_files_path}/{edge_type}.csv'
            print(f'\nLoading {edge_type} edges from {edge_file_path}')
            try:
                edge_df = pd.read_csv(edge_file_path, **kwargs)
            except FileNotFoundError:
                print(f"Error: File not found: {edge_file_path}")
                break
            except pd.errors.EmptyDataError:
                print(f"Warning: File is empty: {edge_file_path}")
                break

            # Map source and destination IDs
            src, unmapped_src_ids = map_edge_ids(edge_df, src_column, src_node_dict, alt_src_node_dict, src_type in ['Disease', 'Phenotype'])
            if src_type in unmatched_node_dict:
                unmatched_node_dict[src_type].update(unmapped_src_ids)
            else:
                unmatched_node_dict[src_type] = unmapped_src_ids

            dst, unmapped_dst_ids = map_edge_ids(edge_df, dst_column, dst_node_dict, alt_dst_node_dict, dst_type in ['Disease', 'Phenotype'])
            if dst_type in unmatched_node_dict:
                unmatched_node_dict[dst_type].update(unmapped_dst_ids)
            else:
                unmatched_node_dict[dst_type] = unmapped_dst_ids

            # Filter out edges with missing node mappings or non-integer values
            src_none = sum(1 for x in src if x is None)
            dst_none = sum(1 for x in dst if x is None)
            src_dst_pairs = list(filter(lambda x: x[0] is not None and x[1] is not None, zip(src, dst)))
            none_count = len(src) - len(src_dst_pairs)
            print(f'First 5 {edge_type} pairs: {src_dst_pairs[:5]}')

            if src_dst_pairs:
                src, dst = zip(*src_dst_pairs)
                edge_dict[edge_type] = (list(src), list(dst))
                print(f'{len(src)} {edge_type} edges loaded from {len(edge_df)} rows')
                print(f'{none_count} edges could not be loaded due to unmatched IDs: {src_none} unmatched {src_type} IDs, {dst_none} unmatched {dst_type} IDs')

            else:
                edge_dict[edge_type] = ([], [])
                print(f'No valid {edge_type} edges loaded from {len(edge_df)} rows')
                break
        return edge_dict, unmatched_node_dict


    def load_data(self, node_dict, edge_dict, node_types, edge_types):
        data = HeteroData()
        for node_type in node_types:
            data[node_type].mapping = node_dict[node_type]['mapping']
            data[node_type].num_nodes = node_dict[node_type]['num_nodes']
        for edge_type in edge_types:
            src, dst = edge_dict[edge_type[0]]
            data[(edge_type[1], edge_type[0], edge_type[2])].edge_index = torch.tensor([src, dst], dtype=torch.long)
            
        data = ToUndirected(merge=False)(data)
        return data
    
    
if __name__ == '__main__':
    with open(args.config_path, 'r') as f:
        config = yaml.safe_load(f)    
    graph_loader = GraphLoader(args.database_files_path, config)
    node_dict = graph_loader.load_node_csvs(low_memory=False)
    edge_dict, unmatched_node_dict = graph_loader.load_edge_csvs(low_memory=False)

    for node_type, unmapped_ids in unmatched_node_dict.items():
        unmatched_ids_df = pd.DataFrame({'node_id': list(unmapped_ids)})
        os.makedirs(f'{args.output_path}/unmatched_node_lists', exist_ok=True)
        unmatched_ids_df.to_csv(f'{args.output_path}/unmatched_node_lists/unmatched_{node_type}.csv', index=False)
    print(f'Unmatched node lists saved to {args.output_path}/unmatched_node_lists')

    data = graph_loader.load_data(node_dict, edge_dict, config['node_types'], config['edge_types'])
    output_file_path = f'{args.output_path}/crossbar_heterodata.pt'
    print(f'Loaded data into HeteroData object:')
    print(data.edge_index_dict)

    torch.save(data, output_file_path)
    print(f'HeteroData object saved to {output_file_path}')