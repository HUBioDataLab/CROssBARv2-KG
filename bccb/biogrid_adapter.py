import os
import sys
import pandas as pd
import numpy as np
import time
import collections
from typing import Union

from pathlib import Path
from time import time


from pypath.inputs import biogrid, uniprot
from pypath.share import curl, settings

from tqdm import tqdm # progress bar

from biocypher._logger import logger
from pypath.resources import urls
from contextlib import ExitStack

from bioregistry import normalize_curie

from enum import Enum

global adapter_name
adapter_name = "biogrid" # might be useful for future


class BiogridEdgeFields(Enum):
    SOURCE = "source"
    PUBMED_IDS = "pmid"
    EXPERIMENTAL_SYSTEM = "experimental_system"
    

class BioGRID:
    def __init__(self, output_dir = None, export_csvs = False, split_output = False, cache=False, debug=False, retries=6,
                organism=9606, biogrid_fields: Union[None, list[BiogridEdgeFields]] = None, add_prefix = True, test_mode = False, aggregate_pubmed_ids: bool = True,
                aggregate_methods: bool = True):
        """
        Downloads and processes BioGRID data

            Args:
                export_csvs: Flag for whether or not create csvs of outputs of databases
                split_csvs: whether or not to split output csv files to multiple parts
                cache: if True, it uses the cached version of the data, otherwise
                forces download.
                debug: if True, turns on debug mode in pypath.
                retries: number of retries in case of download error.
                organism: taxonomy id number of selected organism, if it is None, downloads all organism data.
                biogrid_fields: biogrid fields to be used in the graph.
                add_prefix: if True, add prefix to uniprot ids
                test_mode: if True, take small sample from data for testing
                aggregate_pubmed_ids: if True, collects all pubmed ids that belongs to same protein pair
                aggregate_methods = if True, collects all experiemental methods that belongs to same protein pair
                
        """
        
        self.export_csvs = export_csvs
        self.split_output = split_output
        # Swiss-Prot accession list is downloaded lazily in download_biogrid_data()
        # so that merely constructing the adapter (e.g. for configuration) does not
        # trigger a large network download.
        self.swissprots = None
        self.cache = cache
        self.debug = debug
        self.retries = retries        
        self.organism = organism
        self.biogrid_fields = biogrid_fields
        self.add_prefix = add_prefix
        self.test_mode = test_mode
        
        self.aggregate_dict = {BiogridEdgeFields.PUBMED_IDS.value:aggregate_pubmed_ids,
                              BiogridEdgeFields.EXPERIMENTAL_SYSTEM.value:aggregate_methods}

        
        if export_csvs:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            self.output_dir = output_dir

    
    def export_dataframe(self, dataframe, data_label):
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # TODO: activate this block after adding n_rows_in_file setting to config file
        # if self.split_output:
        #     n_chunks = round(len(dataframe) / n_rows_in_file)
            # output_path =  os.path.join(self.output_dir, data_label)
            # for id, chunk in  enumerate(np.array_split(dataframe, n_chunks)):
            #     chunk.to_csv(os.path.join(output_path, f"crossbar_ppi_data_{data_label}_{id+1}.csv"), index=False)
        # else:
        output_path = os.path.join(self.output_dir, f"crossbar_ppi_data_{data_label}.csv")
        dataframe.to_csv(output_path, index=False)

        return output_path

                         
    def download_biogrid_data(self):
        """
        Wrapper function to download BioGRID data using pypath; used to access
        settings.
            
        To do: Make arguments of biogrid.biogrid_all_interactions selectable for user. 
        """
        
        logger.info(f"This is the link of BioGRID data we downloaded:{urls.urls['biogrid']['all']}. Please check if it is up to date")    
        logger.debug("Started downloading BioGRID data")
        t0 = time()

        with ExitStack() as stack:                         
            stack.enter_context(settings.context(retries=self.retries))
            
            if self.debug:                
                stack.enter_context(curl.debug_on())
            if not self.cache:
                stack.enter_context(curl.cache_off())

            # download biogrid data
            self.biogrid_ints = biogrid.biogrid_all_interactions(self.organism, 9999999999, False)
                        
            # download these fields for mapping from gene symbol to uniprot id
            self.uniprot_to_gene = uniprot.uniprot_data("gene_names", "*", True)
            self.uniprot_to_tax = uniprot.uniprot_data("organism_id", "*", True)

            # swiss-prot accession list (used as a filter in biogrid_process);
            # downloaded here (lazily) rather than in __init__.
            self.swissprots = list(uniprot._all_uniprots("*", True))


        if self.test_mode:
            # Bound the number of raw records so processing stays fast, but keep enough
            # rows to survive gene->uniprot mapping and the Swiss-Prot filter. The final
            # sample is trimmed to a small representative set at the end of
            # biogrid_process() (after mapping/filtering), so it is never empty.
            self.biogrid_ints = self.biogrid_ints[:10000]

        t1 = time()
        logger.info(f'BioGRID data is downloaded in {round((t1-t0) / 60, 2)} mins')
                         

    def biogrid_process(self, rename_selected_fields: Union[None, list[str]] = None) -> None:
        """
        Processor function for BioGRID data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. In addition, it
        maps entries to uniprot ids using gene name and tax id information in the BioGRID data. Also, it filters protein pairs found in swissprot.
        
         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        
        selected_fields = self.set_edge_fields()
            
        default_field_names = {"source":"source", "pmid":"pubmed_ids", "experimental_system":"method"}
        
        self.biogrid_field_new_names = {}
        
        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception("Length of selected_fields variable should be equal to length of rename_selected_fields variable")
            
            for field_old_name, field_new_name in list(zip(selected_fields, rename_selected_fields)):
                self.biogrid_field_new_names[field_old_name] = field_new_name
            
            self.biogrid_field_new_names["uniprot_a"] = "uniprot_a"
            self.biogrid_field_new_names["uniprot_b"] = "uniprot_b"
        else:
            for field_old_name in selected_fields:
                self.biogrid_field_new_names[field_old_name] = default_field_names[field_old_name]
            
            self.biogrid_field_new_names["uniprot_a"] = "uniprot_a"
            self.biogrid_field_new_names["uniprot_b"] = "uniprot_b"
        
        
        logger.debug("Started processing BioGRID data")
        t1 = time()
                         
        # create dataframe          
        biogrid_df = pd.DataFrame.from_records(self.biogrid_ints, columns=self.biogrid_ints[0]._fields)

        # biogrid id (gene symbols) to uniprot id mapping
        biogrid_df['partner_a'] = biogrid_df['partner_a'].str.upper()
        biogrid_df['partner_b'] = biogrid_df['partner_b'].str.upper()
                         
        gene_to_uniprot = collections.defaultdict(list)
        for k,v in self.uniprot_to_gene.items():
            for gene in v.split():
                gene_to_uniprot[gene.upper()].append(k)

        prot_a_uniprots = []
        for prot, tax in zip(biogrid_df['partner_a'], biogrid_df['tax_a']):
            uniprot_id_a = (
                ";".join([_id for _id in gene_to_uniprot[prot] if tax == self.uniprot_to_tax[_id]])
                    if prot in gene_to_uniprot else
                None)
            prot_a_uniprots.append(uniprot_id_a)

        prot_b_uniprots = []
        for prot, tax in zip(biogrid_df['partner_b'], biogrid_df['tax_b']):
            uniprot_id_b = (
                ";".join([_id for _id in gene_to_uniprot[prot] if tax == self.uniprot_to_tax[_id]])
                    if prot in gene_to_uniprot else
                None)
            prot_b_uniprots.append(uniprot_id_b)

        biogrid_df["uniprot_a"] = prot_a_uniprots
        biogrid_df["uniprot_b"] = prot_b_uniprots
        
        biogrid_df.fillna(value=np.nan, inplace=True)
        
        # add source database info
        biogrid_df["source"] = "BioGRID"
        # filter selected fields
        biogrid_df = biogrid_df[list(self.biogrid_field_new_names.keys())]
        # rename columns
        biogrid_df.rename(columns=self.biogrid_field_new_names, inplace=True)
        
        # drop rows that have semicolon (";")
        biogrid_df.drop(biogrid_df[(biogrid_df["uniprot_a"].str.contains(";")) | (biogrid_df["uniprot_b"].str.contains(";"))].index, axis=0, inplace=True)
        biogrid_df.reset_index(drop=True, inplace=True)
        
        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        biogrid_df = biogrid_df[(biogrid_df["uniprot_a"].isin(self.swissprots)) & (biogrid_df["uniprot_b"].isin(self.swissprots))]
        biogrid_df.reset_index(drop=True, inplace=True)
        
        # drop duplicates if same a x b pair exists multiple times 
        # keep the first pair and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same experimental system type with b x a pair, drop b x a pair
        biogrid_df_unique = biogrid_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)        
        
        
        def aggregate_fields(element):
            element = "|".join([str(e) for e in set(element.dropna())])
            if not element:
                return np.nan
            else:
                return element
            
        # Build the aggregation map for every selected field. Fields whose aggregate
        # flag is True are collapsed with aggregate_fields ("|"-joined); all others keep
        # the first value. This is built unconditionally: previously agg_dict was only
        # created inside `if any(self.aggregate_dict.values())`, which raised a NameError
        # when both aggregate_pubmed_ids and aggregate_methods were False.
        agg_field_list = [k for k, v in self.aggregate_dict.items() if v]

        agg_dict = {}
        for k, v in self.biogrid_field_new_names.items():
            if k in agg_field_list:
                agg_dict[v] = aggregate_fields
            else:
                agg_dict[v] = "first"

        biogrid_df_unique = biogrid_df_unique.groupby(["uniprot_a", "uniprot_b"], sort=False, as_index=False).aggregate(agg_dict)
        #biogrid_df_unique["pubmed_id"].replace("", np.nan, inplace=True)
        
        if "experimental_system" in self.biogrid_field_new_names.keys():            
            biogrid_df_unique = biogrid_df_unique[~biogrid_df_unique[["uniprot_a", "uniprot_b", self.biogrid_field_new_names["experimental_system"]]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)
        else:
            biogrid_df_unique = biogrid_df_unique[~biogrid_df_unique[["uniprot_a", "uniprot_b"]].apply(frozenset, axis=1).duplicated()].reset_index(drop=True)

        # In test mode, trim to a small representative sample AFTER mapping, the
        # Swiss-Prot filter and deduplication, so the sample is never empty (unlike
        # trimming the raw records before mapping in download_biogrid_data()).
        if self.test_mode:
            biogrid_df_unique = biogrid_df_unique.head(100).reset_index(drop=True)

        if self.export_csvs:
            biogrid_output_path = self.export_dataframe(biogrid_df_unique, "biogrid")
            logger.info(f'Final BioGRID data is written: {biogrid_output_path}')

        self.final_biogrid_ints = biogrid_df_unique
        
        t2 = time()
        logger.info(f'BioGRID data is processed in {round((t2-t1) / 60, 2)} mins')
            
    def set_edge_fields(self) -> list:
        """
        Sets biogrid edge fields
        Returns:
            selected field list
        """
        if self.biogrid_fields is None:
            return [field.value for field in BiogridEdgeFields]
        else:
            return [field.value for field in self.biogrid_fields]
        
    def add_prefix_to_id(self, prefix="uniprot", identifier=None, sep=":") -> str:
        """
        Adds prefix to uniprot id
        """
        if self.add_prefix and identifier:
            return normalize_curie( prefix + sep + str(identifier))
        
        return identifier
        
        
    def get_biogrid_edges(self) -> list:
        """
        Get PPI edges from biogrid data
        """
        
        # create edge list
        edge_list = []
        
        for index, row in tqdm(self.final_biogrid_ints.iterrows(), total=self.final_biogrid_ints.shape[0]):
            _dict = row.to_dict()
            
            _source = self.add_prefix_to_id(identifier = _dict["uniprot_a"])
            _target = self.add_prefix_to_id(identifier = _dict["uniprot_b"])
            
            del _dict["uniprot_a"], _dict["uniprot_b"]
            
            _props = dict()
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        # if column has multiple entries create list
                        _props[str(k).replace(" ","_").lower()] = v.split("|")
                    else:
                        _props[str(k).replace(" ","_").lower()] = v
           

            edge_list.append((None, _source, _target, "Interacts_With", _props))
            
            
        return edge_list
