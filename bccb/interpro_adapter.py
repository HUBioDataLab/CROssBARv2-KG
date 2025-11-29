from __future__ import annotations

import os
import h5py
import pandas as pd
import numpy as np

from pypath.share import curl, settings
from pypath.inputs import interpro

from contextlib import ExitStack
from bioregistry import normalize_curie
from tqdm import tqdm

from time import time
from biocypher._logger import logger

from typing import Literal, Union, Optional
from pydantic import BaseModel, DirectoryPath, FilePath, validate_call

from enum import Enum, EnumMeta

logger.debug(f"Loading module {__name__}.")


class InterProEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class InterProNodeField(Enum, metaclass=InterProEnumMeta):
    """
    Domain node fields in InterPro
    """

    # primary attributes
    PROTEIN_COUNT = "protein_count"
    NAME = "name"
    TYPE = "type"
    PARENT_LIST = "parent_list"
    CHILD_LIST = "child_list"

    # member list attributes
    PFAM = "PFAM"

    # external attributes
    EC = "EC"

    # structural attributes
    PDB = "PDB"

    # embedding
    DOM2VEC_EMBEDDING = "dom2vec_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def get_primary_attributes(cls):
        """
        Returns primary InterPro attributes
        """
        return [
            cls.PROTEIN_COUNT.value,
            cls.NAME.value,
            cls.TYPE.value,
            cls.PARENT_LIST.value,
            cls.CHILD_LIST.value,
        ]

    @classmethod
    def get_member_list_attributes(cls):
        """
        Returns external InterPro attributes
        """
        return [cls.PFAM.value]

    @classmethod
    def get_external_attributes(cls):
        """
        Returns external InterPro attributes
        """
        return [cls.EC.value]

    @classmethod
    def get_structural_attributes(cls):
        """
        Returns structural InterPro attributes
        """
        return [cls.PDB.value]


class InterProEdgeField(Enum, metaclass=InterProEnumMeta):
    """
    Domain edge fields in InterPro
    """

    START = "start"
    END = "end"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class InterProModel(BaseModel):
    page_size: int = 150
    organism: int | Literal["*"] | None = None
    add_prefix: bool = True
    node_fields: Union[list[InterProNodeField], None] = None
    edge_fields: Union[list[InterProEdgeField], None] = None
    test_mode: bool = False


class InterPro:
    """
    Class that downloads InterPro data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(
        self,
        page_size: Optional[int] = 150,
        organism: Optional[int | Literal["*"] | None] = None,
        add_prefix: Optional[bool] = True,
        node_fields: Optional[Union[list[InterProNodeField], None]] = None,
        edge_fields: Optional[Union[list[InterProEdgeField], None]] = None,
        test_mode: Optional[bool] = False,
    ):
        """
        Args:
            retries: number of retries in case of download error.
            page_size: page size of downloaded annotation data
            organism: rganism code in NCBI taxid format, e.g. "9606" for human. If it is None or "*", downloads all organism data.
            add_prefix: if True, add prefix to database identifiers
            node_fields: `InterProNodeField` fields to be used in the graph, if it is None, select all fields.
            edge_fields: `InterProEdgeField` fields to be used in the graph, if it is None, select all fields.
            test_mode: limits amount of data for testing
        """

        model = InterProModel(
            page_size=page_size,
            organism=organism,
            add_prefix=add_prefix,
            node_fields=node_fields,
            edge_fields=edge_fields,
            test_mode=test_mode,
        ).model_dump()

        self.page_size = model["page_size"]
        self.organism = (
            None if model["organism"] in ("*", None) else model["organism"]
        )
        self.add_prefix = model["add_prefix"]

        # set node and edge fields
        self.set_node_and_edge_fields(
            node_fields=model["node_fields"], edge_fields=model["edge_fields"]
        )

        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100
            
    @validate_call
    def download_interpro_data(self, cache: bool = False,
                               debug: bool = False,
                               retries: int = 3,
                               dom2vec_embedding_path: FilePath | None = None) -> None:
        """
        Wrapper function to download InterPro data using pypath; used to access
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

            self.download_domain_node_data(dom2vec_embedding_path=dom2vec_embedding_path)
            self.download_domain_edge_data()

    @validate_call
    def download_domain_node_data(self, dom2vec_embedding_path: FilePath | None = None) -> None:
        """
        Downloads domain node data from Interpro
        """

        logger.debug("Started downloading InterPro domain data")
        t0 = time()

        # To do: Filter according to tax id??
        self.interpro_entries = (
            interpro.interpro_entries()
        )  # returns a list of namedtuples
        self.interpro_structural_xrefs = interpro.interpro_xrefs(
            db_type="structural"
        )
        self.interpro_external_xrefs = interpro.interpro_xrefs(
            db_type="external"
        )

        if InterProNodeField.DOM2VEC_EMBEDDING.value in self.node_fields:
            self.retrieve_dom2vec_embeddings(dom2vec_embedding_path=dom2vec_embedding_path)

        t1 = time()
        logger.info(
            f"InterPro domain data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def download_domain_edge_data(self) -> None:
        """
        Downloads Uniprot annotation data from Interpro
        """

        logger.debug("Started downloading InterPro annotation data")
        t0 = time()

        if self.organism:
            # WARNING: decrease page_size parameter if there is a curl error about timeout in the pypath_log
            self.interpro_annotations = interpro.interpro_annotations(
                page_size=self.page_size,
                reviewed=True,
                tax_id=self.organism,
            )
        else:
            self.interpro_annotations = interpro.interpro_annotations(
                page_size=self.page_size, reviewed=True, tax_id=""
            )

        t1 = time()
        logger.info(
            f"InterPro annotation data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def retrieve_dom2vec_embeddings(self, 
                                    dom2vec_embedding_path: FilePath | None = None) -> None:

        logger.info("Retrieving dom2vec domain embeddings.")

        self.interpro_id_to_dom2vec_embedding = {}
        with h5py.File(dom2vec_embedding_path, "r") as f:
            for interpro_id, embedding in f.items():
                self.interpro_id_to_dom2vec_embedding[interpro_id] = np.array(embedding).astype(np.float16)
                    
    @validate_call
    def get_interpro_nodes(self, node_label: str = "domain") -> list[tuple]:
        """
        Prepares InterPro domain nodes for BioCypher
        Args:
            node_label : label of interpro nodes
        """

        if not hasattr(self, "interpro_entries"):
            self.download_domain_node_data()

        # create list of nodes
        node_list = []

        # define primary and external attributes
        primary_attributes = InterProNodeField.get_primary_attributes()
        member_list_attributes = InterProNodeField.get_member_list_attributes()
        external_attributes = InterProNodeField.get_external_attributes()
        structural_attributes = InterProNodeField.get_structural_attributes()

        logger.debug("Creating domain nodes")
        t0 = time()

        # set counter for early stopping
        counter = 0

        for entry in tqdm(self.interpro_entries):
            props = {}
            interpro_props = entry._asdict()

            domain_id = self.add_prefix_to_id("interpro", entry.interpro_id)

            # get primary InterPro attributes
            for element in primary_attributes:

                if element in self.node_fields and interpro_props.get(element):
                    if element == InterProNodeField.PROTEIN_COUNT.value:
                        props[element.replace(" ", "_").lower()] = int(
                            interpro_props.get(element)
                        )
                    elif element == InterProNodeField.NAME.value:
                        props[element.replace(" ", "_").lower()] = interpro_props.get(element).replace("'","^").replace("|",",") if interpro_props.get(element) else None
                    else:
                        props[element.replace(" ", "_").lower()] = (
                            self.check_length(interpro_props.get(element))
                        )

            # get member list InterPro attributes
            for element in member_list_attributes:
                if element in self.node_fields and interpro_props.get(
                    "member_list"
                ).get(element):
                    props[element.replace(" ", "_").lower()] = (
                        self.check_length(
                            interpro_props.get("member_list").get(element)
                        )
                    )

            # get external InterPro attributes
            for element in external_attributes:
                if (
                    element in self.node_fields
                    and self.interpro_external_xrefs.get(entry.interpro_id).get(
                        element
                    )
                ):
                    props[element.replace(" ", "_").lower()] = (
                        self.check_length(
                            self.interpro_external_xrefs.get(
                                entry.interpro_id
                            ).get(element)
                        )
                    )

            # get structural InterPro attributes
            for element in structural_attributes:
                if (
                    element in self.node_fields
                    and self.interpro_structural_xrefs.get(
                        entry.interpro_id
                    ).get(element)
                ):
                    props[element.replace(" ", "_").lower()] = (
                        self.check_length(
                            self.interpro_structural_xrefs.get(
                                entry.interpro_id
                            ).get(element)
                        )
                    )

            # get dom2vec embedding
            if InterProNodeField.DOM2VEC_EMBEDDING.value in self.node_fields and self.interpro_id_to_dom2vec_embedding.get(entry.interpro_id) is not None:
                props[InterProNodeField.DOM2VEC_EMBEDDING.value] = [str(emb) for emb in self.interpro_id_to_dom2vec_embedding[entry.interpro_id]]

            # add node to list
            node_list.append((domain_id, node_label, props))

            counter += 1

            if self.early_stopping and counter == self.early_stopping:
                break
            node_list.append((domain_id, node_label, props))

            counter += 1

            if self.early_stopping and counter == self.early_stopping:
                break

        t1 = time()
        logger.info(f"InterPro nodes created in {round((t1-t0) / 60, 2)} mins")

        return node_list

    @validate_call
    def get_interpro_edges(
        self, edge_label: str = "protein_has_domain"
    ) -> list[tuple]:
        """
        Prepares Protein-Domain edges for BioCypher
        Args:
            edge_label: label of protein-domain edge
        """

        if not hasattr(self, "interpro_annotations"):
            self.download_domain_edge_data()

        # create list of edges
        edge_list = []

        logger.debug("Creating protein-domain edges")
        t0 = time()

        # set counter for early stopping
        counter = 0

        # DOMAIN-PROTEIN EDGES
        for k, v in tqdm(self.interpro_annotations.items()):
            # k -> uniprot id
            for annotation in v:

                interpro_props = annotation._asdict()
                props = {}

                for field in self.edge_fields:
                    if interpro_props.get(field, None):
                        props[field.replace(" ", "_").lower()] = (
                            self.check_length(interpro_props[field])
                        )

                interpro_id = self.add_prefix_to_id(
                    "interpro", annotation.interpro_id
                )
                uniprot_id = self.add_prefix_to_id("uniprot", k)

                edge_list.append(
                    (None, uniprot_id, interpro_id, edge_label, props)
                )

                counter += 1

            if self.early_stopping and counter >= self.early_stopping:
                break

        t1 = time()
        logger.info(f"InterPro edges created in {round((t1-t0) / 60, 2)} mins")

        return edge_list

    @validate_call
    def check_length(self, element: str | list | int | float) -> str | list | int | float:
        """
        If the type of given entry is a list and has just one element returns this one element
        """
        if isinstance(element, list) and len(element) == 1:
            return element[0]
        else:
            return element

    @validate_call
    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to ids
        """
        if self.add_prefix:
            return normalize_curie(prefix + sep + identifier)

        return identifier

    def set_node_and_edge_fields(self, node_fields, edge_fields) -> None:
        """
        Sets Interpro node and edge fields
        """

        if node_fields:
            self.node_fields = [field.value for field in node_fields]
        else:
            self.node_fields = [field.value for field in InterProNodeField]

        if edge_fields:
            self.edge_fields = [field.value for field in edge_fields]
        else:
            self.edge_fields = [field.value for field in InterProEdgeField]

    def export_as_csv(self, path: DirectoryPath | None = None):
        if path:
            node_full_path = os.path.join(path, "Domain.csv")
            edge_full_path = os.path.join(path, "Protein_has_domain.csv")
        else:
            node_full_path = os.path.join(os.getcwd(), "Domain.csv")
            edge_full_path = os.path.join(os.getcwd(), "Protein_has_domain.csv")

        # write nodes
        nodes = self.get_interpro_nodes()
        node_df_list = []
        for n in nodes:
            props = {"id": n[0]}
            props |= n[2]
            node_df_list.append(props)

        nodes_df = pd.DataFrame.from_records(node_df_list)
        nodes_df.to_csv(node_full_path, index=False)
        logger.info(f"Domain node data is written: {node_full_path}")

        # write edges
        edges = self.get_interpro_edges()
        edges_df_list = []
        for e in edges:
            props = {"source_id": e[1], "target_id": e[2]}
            props |= e[4]
            edges_df_list.append(props)

        edges_df = pd.DataFrame.from_records(edges_df_list)
        edges_df.to_csv(edge_full_path, index=False)
        logger.info(f"Domain edge data is written: {edge_full_path}")
