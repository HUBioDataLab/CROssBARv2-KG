from time import time
import collections
from typing import Optional, Union, Literal
from collections.abc import Generator
from enum import Enum, EnumMeta, auto
from functools import lru_cache
import pandas as pd
import numpy as np

import os
import requests
import h5py

from tqdm import tqdm  # progress bar
from pypath.share import curl, settings
from pypath.utils import mapping
from pypath.inputs import uniprot
from biocypher._logger import logger
from contextlib import ExitStack
from bioregistry import normalize_curie

from pydantic import BaseModel, DirectoryPath, FilePath, HttpUrl, validate_call

logger.debug(f"Loading module {__name__}.")


class UniprotEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class UniprotNodeType(Enum, metaclass=UniprotEnumMeta):
    """
    Node types of the UniProt API represented in this adapter.
    """

    PROTEIN = auto()
    GENE = auto()
    ORGANISM = auto()
    CELLULAR_COMPARTMENT = auto()


class UniprotNodeField(Enum, metaclass=UniprotEnumMeta):
    """
    Fields of nodes the UniProt API represented in this adapter. Overview of
    uniprot fields: https://www.uniprot.org/help/return_fields
    """

    # core attributes
    LENGTH = "length"
    SUBCELLULAR_LOCATION = "subcellular_location"
    MASS = "mass"
    ORGANISM = "organism_name"
    ORGANISM_ID = "organism_id"
    PROTEIN_NAMES = "protein_name"
    EC = "ec"
    PROTEIN_GENE_NAMES = "gene_names"
    PRIMARY_GENE_NAME = "gene_primary"
    SEQUENCE = "sequence"

    # xref attributes
    ENSEMBL_TRANSCRIPT_IDS = "xref_ensembl"
    PROTEOME = "xref_proteomes"
    ENTREZ_GENE_IDS = "xref_geneid"
    KEGG_IDS = "xref_kegg"

    # not from uniprot REST
    # we provide these by mapping ENSTs via pypath
    ENSEMBL_GENE_IDS = "ensembl_gene_ids"

    # not from uniprot REST
    # we provide these by downloading the ProtT5 embeddings from uniprot
    PROTT5_EMBEDDING = "prott5_embedding"

    # not from uniprot REST
    # we provide these by getting embeddings from ESM2 650M model
    ESM2_EMBEDDING = "esm2_embedding"

    # not from uniprot REST
    # we provide these by getting embeddings from Nucletide tranformer 2 model
    NT_EMBEDDING = "nt_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None
    @classmethod
    def get_protein_properties(cls):
        return [
            cls.LENGTH.value,
            cls.MASS.value,
            cls.PROTEIN_NAMES.value,
            cls.PROTEOME.value,
            cls.EC.value,
            cls.ORGANISM.value,
            cls.ORGANISM_ID.value,
            cls.SEQUENCE.value,
            cls.PROTT5_EMBEDDING.value,
            cls.ESM2_EMBEDDING.value,
        ]
    @classmethod
    def get_gene_properties(cls) -> list:
        return [
            cls.PROTEIN_GENE_NAMES.value,
            cls.ENTREZ_GENE_IDS.value,
            cls.KEGG_IDS.value,
            cls.ENSEMBL_TRANSCRIPT_IDS.value,
            cls.ENSEMBL_GENE_IDS.value,
            cls.PRIMARY_GENE_NAME.value,
            cls.NT_EMBEDDING.value,
        ]
    
    @classmethod
    def get_organism_properties(cls) -> list:
        return [cls.ORGANISM.value]
    
    @classmethod
    def get_split_fields(cls) -> list:
        return [
            cls.PROTEOME.value,
            cls.PROTEIN_GENE_NAMES.value,
            cls.EC.value,
            cls.ENTREZ_GENE_IDS.value,
            cls.ENSEMBL_TRANSCRIPT_IDS.value,
            cls.KEGG_IDS.value,
        ]

class UniprotEdgeType(Enum, metaclass=UniprotEnumMeta):
    """
    Edge types of the UniProt API represented in this adapter.
    """

    PROTEIN_TO_ORGANISM = auto()
    GENE_TO_PROTEIN = auto()


class UniprotIDField(Enum, metaclass=UniprotEnumMeta):
    """
    Fields of edges of the UniProt API represented in this adapter. Used to
    assign source and target identifiers in `get_edges()`.
    """

    # default
    PROTEIN_UNIPROT_ACCESSION = auto()
    GENE_ENTREZ_ID = auto()
    ORGANISM_NCBI_TAXONOMY_ID = auto()

    # optional
    GENE_ENSEMBL_GENE_ID = auto()


class UniProtModel(BaseModel):
    organism: Literal["*"] | int | None = "*"
    rev: bool = True
    node_types: Optional[Union[list[UniprotNodeType], None]] = None
    node_fields: Optional[Union[list[UniprotNodeField], None]] = None
    edge_types: Optional[Union[list[UniprotEdgeType], None]] = None
    id_fields: Optional[Union[list[UniprotIDField], None]] = None
    add_prefix: bool = True
    test_mode: bool = False


class Uniprot:
    """
    Class that downloads uniprot data using pypath and reformats it to be ready
    for import into a BioCypher database.

    Args:
        organism: organism code in NCBI taxid format, e.g. "9606" for human.
        rev: if True, it downloads swissprot (i.e., reviewed) entries only.
        node_types: `UniprotNodeType` fields that will be included in graph, if it is None, select all fields.
        node_fields: `UniprotNodeField` fields that will be included in graph, if it is None, select all fields.
        edge_types: `UniprotEdgeType` fields that will be included in graph, if it is None, select all fields.
        id_fields: `UniprotIDField` field that will be included in graph as node identifier, if it is None, selects first 3 fields.
        add_prefix: if True, add prefix to database identifiers.
        test_mode: if True, limits amount of output data.
    """

    def __init__(
        self,
        organism: Optional[Literal["*"] | int | None] = "*",
        rev: Optional[bool] = True,
        node_types: Optional[Union[list[UniprotNodeType], None]] = None,
        node_fields: Optional[Union[list[UniprotNodeField], None]] = None,
        edge_types: Optional[Union[list[UniprotEdgeType], None]] = None,
        id_fields: Optional[Union[list[UniprotIDField], None]] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
    ):
        model = UniProtModel(
            organism=organism,
            rev=rev,
            node_types=node_types,
            node_fields=node_fields,
            edge_types=edge_types,
            id_fields=id_fields,
            add_prefix=add_prefix,
            test_mode=test_mode,
        ).model_dump()

        # params
        self.organism = model["organism"]
        self.rev = model["rev"]
        self.add_prefix = model["add_prefix"]
        self.test_mode = model["test_mode"]

        # provenance
        self.data_source = "uniprot"
        self.data_version = "2024_03"
        self.data_licence = "CC BY 4.0"

        self._configure_fields()

        self._set_node_and_edge_fields(
            node_types=model["node_types"],
            node_fields=model["node_fields"],
            edge_types=model["edge_types"],
        )

        self.set_id_fields(id_fields=model["id_fields"])

        # loading of ligands and receptors sets
        self.ligands = self._read_ligands_set()
        self.receptors = self._read_receptors_set()

        # loading of subcellular locations set
        self.locations = set()

        # gene property name mappings that will be used gene node properties in KG
        self.gene_property_name_mappings = {"gene_primary":"gene_symbol",
                                            "xref_ensembl":"ensembl_transcript_ids",
                                            "xref_kegg":"kegg_ids",
                                            }
        # protein property name mappings that will be used protein node properties in KG
        self.protein_property_name_mappings = {"protein_name":"protein_names"}

    def _read_ligands_set(self) -> set:
        # check if ligands file exists
        if not os.path.isfile("data/ligands_curated.csv"):
            return set()

        ligand_file = pd.read_csv("data/ligands_curated.csv", header=None)
        return set(ligand_file[0])

    def _read_receptors_set(self) -> set:
        # check if receptors file exists
        if not os.path.isfile("data/receptors_curated.csv"):
            return set()

        receptor_file = pd.read_csv("data/receptors_curated.csv", header=None)
        return set(receptor_file[0])

    @validate_call
    def download_uniprot_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
        prott5_embedding_output_path: FilePath | None = None,
        esm2_embedding_path: FilePath | None = None,
        nucleotide_transformer_embedding_path: FilePath | None = None,
    ):
        """
        Wrapper function to download uniprot data using pypath; used to access
        settings.

        Args:
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """

        # stack pypath context managers
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())

            self._download_uniprot_data(
                prott5_embedding_output_path=prott5_embedding_output_path,
                esm2_embedding_path=esm2_embedding_path,
                nucleotide_transformer_embedding_path=nucleotide_transformer_embedding_path,
            )

            # preprocess data
            self._preprocess_uniprot_data()

    @validate_call
    def _download_uniprot_data(
        self, 
        prott5_embedding_output_path: FilePath | None = None,
        esm2_embedding_path: FilePath | None = None,
        nucleotide_transformer_embedding_path: FilePath | None = None,
    ):
        """
        Download uniprot data from uniprot.org through pypath.

        Here is an overview of uniprot return fields:
        https://www.uniprot.org/help/return_fields

        TODO make use of multi-field query
        """

        logger.info("Downloading uniprot data...")

        t0 = time()

        # download all swissprot ids
        self.uniprot_ids = set(uniprot._all_uniprots(self.organism, self.rev))

        # limit to 100 for testing
        if self.test_mode:
            self.uniprot_ids = set(list(self.uniprot_ids)[:100])

        # download attribute dicts
        self.data = {}
        for query_key in tqdm(self.node_fields, desc="Downloading uniprot fields"):
            if query_key in [
                UniprotNodeField.ENSEMBL_GENE_IDS.value,
                UniprotNodeField.PROTT5_EMBEDDING.value,
                UniprotNodeField.ESM2_EMBEDDING.value,
                UniprotNodeField.NT_EMBEDDING.value,
            ]:
                continue

            elif query_key == UniprotNodeField.SUBCELLULAR_LOCATION.value:
                self.data[query_key] = uniprot.uniprot_locations(
                    self.organism, self.rev
                )
            else:
                self.data[query_key] = uniprot.uniprot_data(
                    query_key, self.organism, self.rev
                )

            logger.debug(f"{query_key} field is downloaded")

        # add ensembl gene ids
        self.data[UniprotNodeField.ENSEMBL_GENE_IDS.value] = {}

        if UniprotNodeField.PROTT5_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.PROTT5_EMBEDDING.value] = {}
            self.download_prott5_embeddings(
                prott5_embedding_output_path=prott5_embedding_output_path
            )
        
        if UniprotNodeField.ESM2_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.ESM2_EMBEDDING.value] = {}
            self.retrieve_esm2_embeddings(esm2_embedding_path)

        if UniprotNodeField.NT_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.NT_EMBEDDING.value] = {}
            self.retrieve_nucleotide_transformer_embeddings(
                nucleotide_transformer_embedding_path
            )

        t1 = time()
        msg = f"Acquired UniProt data in {round((t1-t0) / 60, 2)} mins."
        logger.info(msg)

    @validate_call
    def download_prott5_embeddings(
        self, 
        prott5_embedding_output_path: FilePath | None = None
    ):
        """
        Downloads ProtT5 embedding from uniprot website
        If the files exists in a defined file path, then
        directly read it.

        Args:
            prott5_embedding_output_path (FilePath, optional): Defaults to None.
        """
        url: HttpUrl = (
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/uniprot_sprot/per-protein.h5"
        )
        if prott5_embedding_output_path:
            full_path = prott5_embedding_output_path
        else:
            full_path = os.path.join(os.getcwd(), "prott5_protein_embeddings.h5")

        if not os.path.isfile(full_path):
            logger.info("Downloading ProtT5 embeddings...")

            with requests.get(url, stream=True) as response:
                with open(full_path, "wb") as f:
                    for chunk in response.iter_content(512 * 1024):
                        if chunk:
                            f.write(chunk)
        else:
            logger.info("ProtT5 Embedding file is exists. Reading from the file..")

        df_list = []
        with h5py.File(full_path, "r") as file:

            for uniprot_id, embedding in file.items():
                
                embedding = np.array(embedding).astype(np.float16)
                count = np.count_nonzero(~np.isnan(embedding))
                if (
                    self.organism not in ("*", None)
                    and uniprot_id in self.uniprot_ids
                    and count == 1024
                ):
                    # self.data[UniprotNodeField.PROTT5_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))
                elif count == 1024:
                    # self.data[UniprotNodeField.PROTT5_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))

        
        self.prott5_embedding_df = pd.DataFrame(df_list, columns=['uniprot_id', 'embedding'])

        del df_list
                    
    @validate_call
    def retrieve_esm2_embeddings(self, 
                                 esm2_embedding_path: FilePath | None = None) -> None:
        
        logger.info("Retrieving ESM2 embeddings...")

        df_list = []
        with h5py.File(esm2_embedding_path, "r") as file:

            for uniprot_id, embedding in file.items():
                embedding = np.array(embedding).astype(np.float16)
                count = np.count_nonzero(~np.isnan(embedding))
                if (
                    self.organism not in ("*", None)
                    and uniprot_id in self.uniprot_ids
                    and count == 1280
                ):
                    # self.data[UniprotNodeField.ESM2_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))
                elif count == 1280:
                    # self.data[UniprotNodeField.ESM2_EMBEDDING.value][
                    #     uniprot_id
                    # ] = embedding
                    df_list.append((uniprot_id, embedding))

        self.esm2_embedding_df = pd.DataFrame(df_list, columns=['uniprot_id', 'embedding'])

        del df_list

    def retrieve_nucleotide_transformer_embeddings(self,
                                                   nucleotide_transformer_embedding_path: FilePath | None = None) -> None:
        
        logger.info("Retrieving Nucleotide Transformer embeddings...")

        self.entrez_id_to_nucleotide_transformer_embedding = {}
        with h5py.File(nucleotide_transformer_embedding_path, "r") as file:
            for entrez_id, embedding in file.items():
                embedding = np.array(embedding).astype(np.float16)
                self.entrez_id_to_nucleotide_transformer_embedding[entrez_id] = embedding
                
    def _preprocess_uniprot_data(self):
        """
        Preprocess uniprot data to make it ready for import. First, three types
        of processing are applied:
        - nothing is done (for ensembl gene ids, which come from pypath)
        - simple string replacement
        - replace separators in integers and convert to int
        - field splitting

        Then, special treatment is applied to some fields:
        - ensg ids are extracted from the ensembl transcript ids
        - protein names and virus hosts have dedicated normalisation functions
        """

        logger.info("Preprocessing UniProt data.")

        for arg in tqdm(self.node_fields, desc="Processing uniprot fields"):

            # do not process ensembl gene ids (we will get them from pypath)
            # and prott5 embeddings
            if arg in [
                UniprotNodeField.ENSEMBL_GENE_IDS.value,
                UniprotNodeField.PROTT5_EMBEDDING.value,
                UniprotNodeField.ESM2_EMBEDDING.value,
                UniprotNodeField.NT_EMBEDDING.value,
            ]:
                pass

            elif arg in [
                UniprotNodeField.LENGTH.value,
                UniprotNodeField.MASS.value,
                UniprotNodeField.ORGANISM_ID.value,
            ]:
                for protein, attribute_value in self.data.get(arg).items():
                    self.data[arg][protein] = int(
                        str(attribute_value).replace(",", "")
                    )

            elif arg not in self.split_fields:
                if arg != UniprotNodeField.SUBCELLULAR_LOCATION.value:
                    for protein, attribute_value in self.data.get(arg).items():

                        self.data[arg][protein] = (
                            attribute_value.replace("|", ",")
                            .replace("'", "^")
                            .strip()
                        )

            else:

                for protein, attribute_value in self.data.get(arg).items():
                    # Field splitting
                    self.data[arg][protein] = self._split_fields(
                        arg, attribute_value
                    )

            # Special treatment
            # ENST and ENSG ids
            if arg == UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS.value:

                for protein, attribute_value in self.data.get(arg).items():
                    

                    attribute_value, ensg_ids = self._find_ensg_from_enst(
                        attribute_value
                    )
                    
                    # update enst in data dict
                    self.data[UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS.value][
                        protein
                    ] = attribute_value

                    if ensg_ids:
                        # add ensgs to data dict
                        self.data[UniprotNodeField.ENSEMBL_GENE_IDS.value][
                            protein
                        ] = ensg_ids

            # Protein names
            elif arg == UniprotNodeField.PROTEIN_NAMES.value:

                for protein, attribute_value in self.data.get(arg).items():

                    self.data[arg][protein] = self._split_protein_names_field(
                        attribute_value
                    )

            elif arg == UniprotNodeField.SUBCELLULAR_LOCATION.value:
                for protein, attribute_value in self.data.get(arg).items():
                    individual_protein_locations = []
                    for element in attribute_value:
                        loc = (
                            str(element.location)
                            .replace("'", "")
                            .replace("[", "")
                            .replace("]", "")
                            .strip()
                        )
                        individual_protein_locations.append(loc)
                        self.locations.add(loc)

                    self.data[arg][protein] = individual_protein_locations

    @validate_call
    def _get_ligand_or_receptor(self, uniprot_id: str):
        """
        Tell if UniProt protein node is a L, R or nothing.
        """

        uniprot_id = uniprot_id[8:]

        if uniprot_id in self.ligands:
            return "ligand"
        return "receptor" if uniprot_id in self.receptors else "protein"

    @validate_call
    def get_nodes(
        self, 
        ligand_or_receptor: bool = False,
        protein_label: str = "protein",
        gene_label: str = "gene",
        organism_label: str = "organism"
    ) -> Generator[tuple[str, str, dict]]:
        """
        Yield nodes (protein, gene, organism) from UniProt data.
        """

        # raise error if ligand_or_receptor is True but self.ligands or
        # self.receptors are empty
        if ligand_or_receptor and (not self.ligands or not self.receptors):
            raise ValueError(
                "No ligands or receptors found in the 'data' directory. "
                "Please set ligand_or_receptor to False or add the files."
            )

        logger.info(
            "Preparing UniProt nodes of the types "
            f"{[type.name for type in self.node_types]}."
        )
        
        for protein_id, all_props in self._reformat_and_filter_proteins():
            
            if UniprotNodeType.PROTEIN in self.node_types:
                protein_props = self._get_protein_properties(all_props, protein_id.split(":")[1])

                # append protein node to output
                if ligand_or_receptor:
                    ligand_or_receptor = self._get_ligand_or_receptor(protein_id)
                    yield (protein_id, ligand_or_receptor, protein_props)
                else:
                    yield (protein_id, protein_label, protein_props)

            # append gene node to output if desired
            if UniprotNodeType.GENE in self.node_types:

                gene_list = self._get_gene(all_props)

                for gene_id, gene_props in gene_list:

                    if gene_id:
                        yield (gene_id, gene_label, gene_props)

            # append organism node to output if desired
            if UniprotNodeType.ORGANISM in self.node_types:

                organism_id, organism_props = self._get_organism(all_props)

                if organism_id:
                    yield (
                        organism_id,
                        organism_label,
                        organism_props,
                    )

    @validate_call
    def get_edges(self, gene_to_protein_label: str = "Gene_encodes_protein",
                  protein_to_organism_label: str = "Protein_belongs_to_organism") -> Generator[tuple[None, str, str, str, dict]]:
        """
        Get nodes and edges from UniProt data.
        """

        logger.info(
            "Preparing UniProt edges of the types "
            f"{[type.name for type in self.edge_types]}."
        )

        # create lists of edges
        edge_list = []

        # generic properties for all edges for now
        properties = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

        for protein in tqdm(self.uniprot_ids, desc="Getting edges"):

            protein_id = self.add_prefix_to_id("uniprot", protein)

            if UniprotEdgeType.GENE_TO_PROTEIN in self.edge_types:

                type_dict = {
                    UniprotNodeField.ENTREZ_GENE_IDS.value: "ncbigene",
                    UniprotNodeField.ENSEMBL_GENE_IDS.value: "ensembl",
                }

                # find preferred identifier for gene
                if UniprotIDField.GENE_ENTREZ_ID in self.id_fields:

                    id_type = UniprotNodeField.ENTREZ_GENE_IDS.value

                elif UniprotIDField.GENE_ENSEMBL_GENE_ID in self.id_fields:

                    id_type = UniprotNodeField.ENSEMBL_GENE_IDS.value

                if genes := self.data.get(id_type).get(protein):
                    genes = self._ensure_iterable(genes)

                    for gene in genes:

                        if not gene:
                            continue

                        gene_id = self.add_prefix_to_id(
                            type_dict[id_type],
                            gene,
                        )
                        edge_list.append(
                            (
                                None,
                                gene_id,
                                protein_id,
                                gene_to_protein_label,
                                properties,
                            )
                        )

            if UniprotEdgeType.PROTEIN_TO_ORGANISM in self.edge_types:

                # TODO all of this processing in separate function
                # is it even still necessary?

                organism_id = (
                    str(
                        self.data.get(UniprotNodeField.ORGANISM_ID.value).get(
                            protein
                        )
                    )
                    if self.data.get(UniprotNodeField.ORGANISM_ID.value).get(
                        protein
                    )
                    else None
                )

                if organism_id:

                    organism_id = self.add_prefix_to_id(
                        "ncbitaxon", organism_id
                    )
                    edge_list.append(
                        (
                            None,
                            protein_id,
                            organism_id,
                            protein_to_organism_label,
                            properties,
                        )
                    )

        if edge_list:

            return edge_list

    def _reformat_and_filter_proteins(self):
        """
        For each uniprot id, select desired fields and reformat to give a tuple
        containing id and properties. Yield a tuple for each protein.
        """

        for protein in tqdm(self.uniprot_ids):

            protein_id = self.add_prefix_to_id("uniprot", protein)

            _props = {arg: self.data.get(arg).get(protein) for arg in self.node_fields}
            yield protein_id, _props

    @validate_call
    def _get_gene(self, all_props: dict) -> list:
        """
        Get gene node representation from UniProt data per protein. Since one
        protein can have multiple genes, return a list of tuples.
        """

        # if genes and database(GeneID) fields exist, define gene_properties
        if not (
            UniprotNodeField.PROTEIN_GENE_NAMES.value in all_props.keys()
            and UniprotNodeField.ENTREZ_GENE_IDS.value in all_props.keys()
        ):
            return []

        # Find preferred identifier for gene and check if it exists
        if UniprotIDField.GENE_ENTREZ_ID in self.id_fields:

            id_type = UniprotNodeField.ENTREZ_GENE_IDS.value

        elif UniprotIDField.GENE_ENSEMBL_GENE_ID in self.id_fields:

            id_type = UniprotNodeField.ENSEMBL_GENE_IDS.value

        gene_raw = all_props.pop(id_type)

        if not gene_raw:
            return []

        type_dict = {
            UniprotNodeField.ENTREZ_GENE_IDS.value: "ncbigene",
            UniprotNodeField.ENSEMBL_GENE_IDS.value: "ensembl",
        }

        genes = self._ensure_iterable(gene_raw)
        gene_props = {}

        for k in all_props.keys():

            if k not in self.gene_properties:
                continue

            # select parenthesis content in field names and make lowercase
            gene_props[
                (
                    self.gene_property_name_mappings[k]
                    if self.gene_property_name_mappings.get(k)
                    else k.replace(" ", "_").replace("-", "_").lower()
                )
            ] = all_props[k]

            if k == UniprotNodeField.NT_EMBEDDING.value and self.entrez_id_to_nucleotide_transformer_embedding.get(genes[0]) is not None:

                gene_props[k] = [str(emb) for emb in self.entrez_id_to_nucleotide_transformer_embedding[genes[0]]]

        # source, licence, and version fields
        gene_props["source"] = self.data_source
        gene_props["licence"] = self.data_licence
        gene_props["version"] = self.data_version

        gene_list = []

        for gene in genes:

            gene_id = self.add_prefix_to_id(
                type_dict[id_type],
                gene,
            )

            gene_list.append((gene_id, gene_props))

        return gene_list

    @validate_call
    def _get_organism(self, all_props: dict):

        organism_props = {}

        organism_id = self.add_prefix_to_id(
            "ncbitaxon",
            str(all_props.pop(UniprotNodeField.ORGANISM_ID.value)),
        )

        for k in all_props.keys():

            if k in self.organism_properties:
                organism_props[k] = all_props[k]

        # source, licence, and version fields
        organism_props["source"] = self.data_source
        organism_props["licence"] = self.data_licence
        organism_props["version"] = self.data_version

        return organism_id, organism_props

    @validate_call
    def _get_protein_properties(self, all_props: dict, protein_id: str) -> dict:
        protein_props = {}
        
        for k in all_props.keys():

            # define protein_properties
            if k not in self.protein_properties:
                continue

            elif k == UniprotNodeField.PROTEIN_NAMES.value:
                protein_props["primary_protein_name"] = (
                    self._ensure_iterable(all_props[k])[0]
                    if all_props[k]
                    else None
                )

                protein_props[
                    (
                        self.protein_property_name_mappings[k]
                        if self.protein_property_name_mappings.get(k)
                        else k.replace(" ", "_").replace("-", "_")
                    )
                ] = all_props[k]

            elif k == UniprotNodeField.PROTT5_EMBEDDING.value:
                res = self.prott5_embedding_df[self.prott5_embedding_df["uniprot_id"] == protein_id]["embedding"]
                if not res.empty:
                # embedding = np.array(self.prott5_embedding[k]).astype(np.float16)
                    protein_props[k.replace(" ", "_").replace("-", "_")] = [str(emb) for emb in res.values[0]]

            elif k == UniprotNodeField.ESM2_EMBEDDING.value:
                res = self.esm2_embedding_df[self.esm2_embedding_df["uniprot_id"] == protein_id]["embedding"]
                if not res.empty:
                # embedding = np.array(self.esm2_embedding[k]).astype(np.float16)
                    protein_props[k.replace(" ", "_").replace("-", "_")] = [str(emb) for emb in res.values[0]]
            
            else:
                # replace hyphens and spaces with underscore
                protein_props[
                    (
                        self.protein_property_name_mappings[k]
                        if self.protein_property_name_mappings.get(k)
                        else k.replace(" ", "_").replace("-", "_")
                    )
                ] = all_props[k]

        # source, licence, and version fields
        protein_props["source"] = self.data_source
        protein_props["licence"] = self.data_licence
        protein_props["version"] = self.data_version

        return protein_props

    def _split_fields(self, field_key, field_value):
        """
        Split fields with multiple entries in uniprot
        Args:
            field_key: field name
            field_value: entry of the field
        """
        if not field_value:
            return None
        # replace sensitive elements for admin-import
        field_value = (
            field_value.replace("|", ",").replace("'", "^").strip()
        )

        # define fields that will not be splitted by semicolon
        split_dict = {
            UniprotNodeField.PROTEOME.value: ",",
            UniprotNodeField.PROTEIN_GENE_NAMES.value: " ",
        }

        # if field in split_dict split accordingly
        if field_key in split_dict:
            field_value = field_value.split(split_dict[field_key])
            # if field has just one element in the list make it string
            if len(field_value) == 1:
                field_value = field_value[0]

        else:
            field_value = field_value.strip().strip(";").split(";")

                # split colons (":") in kegg field
            if field_key == UniprotNodeField.KEGG_IDS.value:
                _list = [e.split(":")[1].strip() for e in field_value]
                field_value = _list

            # take first element in database(GeneID) field
            if field_key == UniprotNodeField.ENTREZ_GENE_IDS.value:
                field_value = field_value[0]

            # if field has just one element in the list make it string
            if isinstance(field_value, list) and len(field_value) == 1:
                field_value = field_value[0]

        return field_value

    def _split_protein_names_field(self, field_value):
        """
        Split protein names field in uniprot
        Args:
            field_value: entry of the protein names field
        Example:
            "Acetate kinase (EC 2.7.2.1) (Acetokinase)" -> ["Acetate kinase", "Acetokinase"]
        """
        field_value = field_value.replace("|", ",").replace(
            "'", "^"
        )  # replace sensitive elements

        if "[Cleaved" in field_value:
            # discarding part after the "[Cleaved"
            clip_index = field_value.index("[Cleaved")
            protein_names = (
                field_value[:clip_index].replace("(Fragment)", "").strip()
            )

            # handling multiple protein names
            if "(EC" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []

                for name in splitted:
                    if not name.strip().startswith("EC") and not name.strip().startswith("Fragm"):
                        protein_names.append(name.rstrip(")").strip())

            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = [
                    name.rstrip(")").strip()
                    for name in splitted
                    if not name.strip().startswith("Fragm")
                ]
        elif "[Includes" in field_value:
            # discarding part after the "[Includes"
            clip_index = field_value.index("[Includes")
            protein_names = (
                field_value[:clip_index].replace("(Fragment)", "").strip()
            )
            # handling multiple protein names
            if "(EC" in protein_names[0]:

                splitted = protein_names[0].split(" (")
                protein_names = [
                    name.rstrip(")").strip()
                    for name in splitted
                    if not name.strip().startswith("EC")
                    and not name.strip().startswith("Fragm")
                ]
            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = [
                    name.rstrip(")").strip()
                    for name in splitted
                    if not name.strip().startswith("Fragm")
                ]
        elif "(EC" in field_value.replace("(Fragment)", ""):
            splitted = field_value.split(" (")
            protein_names = [
                name.rstrip(")").strip()
                for name in splitted
                if not name.strip().startswith("EC")
                and not name.strip().startswith("Fragm")
            ]
        elif " (" in field_value.replace("(Fragment)", ""):
            splitted = field_value.split(" (")
            protein_names = [
                name.rstrip(")").strip()
                for name in splitted
                if not name.strip().startswith("Fragm")
            ]
        else:
            protein_names = field_value.replace("(Fragment)", "").strip()

        return protein_names

    def _find_ensg_from_enst(self, enst_list):
        """
        take ensembl transcript ids, return ensembl gene ids by using pypath mapping tool

        Args:
            field_value: ensembl transcript list

        """

        enst_list = self._ensure_iterable(enst_list)

        enst_list = [enst.split(" [")[0] for enst in enst_list]

        ensg_ids = set()
        for enst_id in enst_list:
            ensg_id = list(
                mapping.map_name(
                    enst_id.split(".")[0], "enst_biomart", "ensg_biomart"
                )
            )
            ensg_id = ensg_id[0] if ensg_id else None
            if ensg_id:
                ensg_ids.add(ensg_id)

        ensg_ids = list(ensg_ids)

        if len(ensg_ids) == 1:
            ensg_ids = ensg_ids[0]

        if len(enst_list) == 1:
            enst_list = enst_list[0]

        return enst_list, ensg_ids

    @lru_cache
    def _normalise_curie_cached(
        self, prefix: str, identifier: str, sep: str = ":"
    ) -> Optional[str]:
        """
        Wrapper to call and cache `normalize_curie()` from Bioregistry.
        """

        if not self.normalise_curies:
            return identifier

        return normalize_curie(f"{prefix}{sep}{identifier}", sep=sep)

    @lru_cache
    @validate_call
    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)

        return identifier

    def _configure_fields(self):
        # fields that need splitting
        self.split_fields = UniprotNodeField.get_split_fields()

        # properties of nodes
        self.protein_properties = UniprotNodeField.get_protein_properties()

        self.gene_properties = UniprotNodeField.get_gene_properties()

        self.organism_properties = UniprotNodeField.get_organism_properties()

    def _set_node_and_edge_fields(
        self,
        node_types,
        node_fields,
        edge_types,
    ):

        # ensure computation of ENSGs
        if (
            UniprotNodeField.ENSEMBL_GENE_IDS in node_fields
            and UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS not in node_fields
        ):
            node_fields.append(UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS)

        # check which node types and fields to include
        self.node_types = node_types or list(UniprotNodeType)

        self.node_fields = [field.value for field in node_fields] if node_fields else [field.value for field in UniprotNodeField]

        # check which edge types and fields to include
        self.edge_types = edge_types or list(UniprotEdgeType)

    def set_id_fields(self, id_fields):
        self.id_fields = id_fields or list(UniprotIDField)[:3]

    def _ensure_iterable(self, value):
        return [value] if isinstance(value, str) else value

    @validate_call
    def export_data_to_csv(
        self,
        node_data: Generator[tuple[str, str, dict]] = None,
        edge_data: Generator[tuple[None, str, str, str, dict]] = None,
        path: DirectoryPath | None = None,
    ) -> None:
        """
        Save node and edge data to csv
            node_data: output of `get_nodes()` function
            edge_data: output of `get_edges()` function
            path: Directory to save the output csv file
        """
        if node_data:
            logger.debug("Saving uniprot node data as csv")
            node_types_dict = collections.defaultdict(list)
            for _id, _type, props in node_data:
                _dict = {"id": _id} | props
                node_types_dict[_type].append(_dict)

            for _type, values in node_types_dict.items():
                df = pd.DataFrame.from_records(values)
                if path:
                    full_path = os.path.join(path, f"{_type.capitalize()}.csv")
                else:
                    full_path = os.path.join(
                        os.getcwd(), f"{_type.capitalize()}.csv"
                    )

                df.to_csv(full_path, index=False)
                logger.info(
                    f"{_type.capitalize()} data is written: {full_path}"
                )

        if edge_data:
            logger.debug("Saving uniprot edge data as csv")
            edge_types_dict = collections.defaultdict(list)
            for _, source_id, target_id, _type, props in edge_data:
                _dict = {"source_id": source_id, "target_id": target_id} | props
                edge_types_dict[_type].append(_dict)

            for _type, values in edge_types_dict.items():
                df = pd.DataFrame.from_records(values)
                if path:
                    full_path = os.path.join(path, f"{_type.capitalize()}.csv")
                else:
                    full_path = os.path.join(
                        os.getcwd(), f"{_type.capitalize()}.csv"
                    )

                df.to_csv(full_path, index=False)
                logger.info(
                    f"{_type.capitalize()} data is written: {full_path}"
                )
