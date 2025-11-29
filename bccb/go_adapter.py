from pypath.share import curl, settings
from pypath.inputs import interpro, uniprot
from pypath.inputs import go as go_input
from pypath.utils import go as go_util
from contextlib import ExitStack
from bioregistry import normalize_curie

import collections
import os
import h5py
import requests

import pandas as pd
import numpy as np

from typing import Literal, Union, Optional
from pydantic import BaseModel, FilePath, HttpUrl, validate_call

from time import time

from tqdm import tqdm

from biocypher._logger import logger

from enum import Enum, EnumMeta, auto

logger.debug(f"Loading module {__name__}.")


class GOEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class GONodeField(Enum, metaclass=GOEnumMeta):
    NAME = "name"
    ANC2VEC_EMBBEDDING = "anc2vec_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class GOEdgeField(Enum, metaclass=GOEnumMeta):
    REFERENCE = "reference"
    EVIDENCE_CODE = "evidence_code"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class GONodeType(Enum, metaclass=GOEnumMeta):
    PROTEIN = auto()
    CELLULAR_COMPONENT = auto()
    BIOLOGICAL_PROCESS = auto()
    MOLECULAR_FUNCTION = auto()
    DOMAIN = auto()


class GOEdgeType(Enum, metaclass=GOEnumMeta):
    # Protein-Go
    PROTEIN_TO_CELLULAR_COMPONENT = auto()
    PROTEIN_TO_BIOLOGICAL_PROCESS = auto()
    PROTEIN_TO_MOLECULAR_FUNCTION = auto()

    # Domain-GO
    DOMAIN_TO_CELLULAR_COMPONENT = auto()
    DOMAIN_TO_BIOLOGICAL_PROCESS = auto()
    DOMAIN_TO_MOLECULAR_FUNCTION = auto()

    # Go-Go
    CELLULAR_COMPONENT_TO_CELLULAR_COMPONENT = auto()
    BIOLOGICAL_PROCESS_TO_BIOLOGICAL_PROCESS = auto()
    MOLECULAR_FUNCTION_TO_MOLECULAR_FUNCTION = auto()
    BIOLOGICAL_PROCESS_TO_MOLECULAR_FUNCTION = auto()


class ProteinToCellularComponentEdgeLabel(Enum, metaclass=GOEnumMeta):
    LOCATED_IN = "located_in"
    IS_ACTIVE_IN = "is_active_in"
    PART_OF = "part_of"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.PROTEIN_TO_CELLULAR_COMPONENT


class ProteinToBiologicalProcessEdgeLabel(Enum, metaclass=GOEnumMeta):
    INVOLVED_IN = "involved_in"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.PROTEIN_TO_BIOLOGICAL_PROCESS


class ProteinToMolecularFunctionEdgeLabel(Enum, metaclass=GOEnumMeta):
    ENABLES = "enables"
    CONTRIBUTES_TO = "contributes_to"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.PROTEIN_TO_MOLECULAR_FUNCTION


class DomainToCellularComponentEdgeLabel(Enum, metaclass=GOEnumMeta):
    LOCATED_IN = "located_in"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.DOMAIN_TO_CELLULAR_COMPONENT


class DomainToBiologicalProcessEdgeLabel(Enum, metaclass=GOEnumMeta):
    INVOLVED_IN = "involved_in"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.DOMAIN_TO_BIOLOGICAL_PROCESS


class DomainToMolecularFunctionEdgeLabel(Enum, metaclass=GOEnumMeta):
    ENABLES = "enables"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.DOMAIN_TO_MOLECULAR_FUNCTION


class MolecularFunctionToMolecularFunctionEdgeLabel(Enum, metaclass=GOEnumMeta):
    IS_A = "is_a"
    POSITIVELY_REGULATES = "positively_regulates"
    NEGATIVELY_REGULATES = "negatively_regulates"
    PART_OF = "part_of"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.MOLECULAR_FUNCTION_TO_MOLECULAR_FUNCTION


class BiologicalProcessToBiologicalProcessEdgeLabel(Enum, metaclass=GOEnumMeta):
    IS_A = "is_a"
    POSITIVELY_REGULATES = "positively_regulates"
    NEGATIVELY_REGULATES = "negatively_regulates"
    PART_OF = "part_of"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.BIOLOGICAL_PROCESS_TO_BIOLOGICAL_PROCESS


class CellularComponentToCellularComponentEdgeLabel(Enum, metaclass=GOEnumMeta):
    IS_A = "is_a"
    PART_OF = "part_of"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.CELLULAR_COMPONENT_TO_CELLULAR_COMPONENT


class BiologicalProcessToMolecularFunctionEdgeLabel(Enum, metaclass=GOEnumMeta):
    POSITIVELY_REGULATES = "positively_regulates"
    NEGATIVELY_REGULATES = "negatively_regulates"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def neccessary_edge_type(cls):
        return GOEdgeType.BIOLOGICAL_PROCESS_TO_MOLECULAR_FUNCTION


class GOModel(BaseModel):
    organism: Union[int, Literal["*"], None] = None
    node_types: Union[list[GONodeType], None] = None
    go_node_fields: Union[list[GONodeField], None] = None
    edge_types: Union[list[GOEdgeType], None] = None
    go_edge_fields: Union[list[GOEdgeField], None] = None
    edge_labels: Union[
        list[ProteinToCellularComponentEdgeLabel],
        list[ProteinToBiologicalProcessEdgeLabel],
        list[ProteinToMolecularFunctionEdgeLabel],
        list[DomainToCellularComponentEdgeLabel],
        list[DomainToBiologicalProcessEdgeLabel],
        list[DomainToMolecularFunctionEdgeLabel],
        list[MolecularFunctionToMolecularFunctionEdgeLabel],
        list[BiologicalProcessToBiologicalProcessEdgeLabel],
        list[CellularComponentToCellularComponentEdgeLabel],
        list[BiologicalProcessToMolecularFunctionEdgeLabel],
        None,
    ] = None
    add_prefix: bool = True
    test_mode: bool = False
    remove_selected_annotations: list[str] = ["IEA"]


class GO:
    """
    Class that downloads Gene Ontology data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(
        self,
        organism: Optional[int | Literal["*"] | None] = None,
        node_types: Optional[Union[list[GONodeType], None]] = None,
        go_node_fields: Optional[Union[list[GONodeField], None]] = None,
        edge_types: Optional[Union[list[GOEdgeType], None]] = None,
        go_edge_fields: Optional[Union[list[GOEdgeField], None]] = None,
        edge_labels: Union[
            list[ProteinToCellularComponentEdgeLabel],
            list[ProteinToBiologicalProcessEdgeLabel],
            list[ProteinToMolecularFunctionEdgeLabel],
            list[DomainToCellularComponentEdgeLabel],
            list[DomainToBiologicalProcessEdgeLabel],
            list[DomainToMolecularFunctionEdgeLabel],
            list[MolecularFunctionToMolecularFunctionEdgeLabel],
            list[BiologicalProcessToBiologicalProcessEdgeLabel],
            list[CellularComponentToCellularComponentEdgeLabel],
            list[BiologicalProcessToMolecularFunctionEdgeLabel],
            None,
        ] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
        remove_selected_annotations: Optional[list[str]] = ["IEA"],
    ):
        """
        Args:
            organism: ncbi tax id or known name of organism of interest
            node_types: node types that will be included in graph, if defined it must be element(s) from GONodeType enum class (not the value of the element)
            go_node_fields: node properties that will be included in graph, if defined it must be element(s) from GONodeField enum class (not the value of the element)
            edge_types: edge types that will be included in graph, if defined it must be element(s) from GOEdgeType enum class (not the value of the element)
            go_edge_fields: protein-go edge properties that will be included in graph, if defined it must be element(s) from GOEdgeField enum class (not the value of the element)
            edge_labels: edge labels that will be included in graph, if defined it must be element(s) from ...EdgeLabel enum classes
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of data for testing
            remove_selected_annotations: removes selected annotations from protein-go edges, by default it removes electronic annotations
        """
        model = GOModel(
            organism=organism,
            node_types=node_types,
            go_node_fields=go_node_fields,
            edge_types=edge_types,
            go_edge_fields=go_edge_fields,
            edge_labels=edge_labels,
            add_prefix=add_prefix,
            test_mode=test_mode,
            remove_selected_annotations=remove_selected_annotations,
        ).model_dump()

        self.organism = (
            "*" if model["organism"] in ("*", None) else model["organism"]
        )
        self.add_prefix = model["add_prefix"]
        self.remove_selected_annotations = model["remove_selected_annotations"]

        # for checking source and target node types of selected edge types
        self.check_node_types_of_edges = {
            GOEdgeType.PROTEIN_TO_CELLULAR_COMPONENT: [
                GONodeType.PROTEIN,
                GONodeType.CELLULAR_COMPONENT,
            ],
            GOEdgeType.PROTEIN_TO_BIOLOGICAL_PROCESS: [
                GONodeType.PROTEIN,
                GONodeType.BIOLOGICAL_PROCESS,
            ],
            GOEdgeType.PROTEIN_TO_MOLECULAR_FUNCTION: [
                GONodeType.PROTEIN,
                GONodeType.MOLECULAR_FUNCTION,
            ],
            GOEdgeType.DOMAIN_TO_CELLULAR_COMPONENT: [
                GONodeType.DOMAIN,
                GONodeType.CELLULAR_COMPONENT,
            ],
            GOEdgeType.DOMAIN_TO_BIOLOGICAL_PROCESS: [
                GONodeType.DOMAIN,
                GONodeType.BIOLOGICAL_PROCESS,
            ],
            GOEdgeType.DOMAIN_TO_MOLECULAR_FUNCTION: [
                GONodeType.DOMAIN,
                GONodeType.MOLECULAR_FUNCTION,
            ],
            GOEdgeType.BIOLOGICAL_PROCESS_TO_BIOLOGICAL_PROCESS: [
                GONodeType.BIOLOGICAL_PROCESS
            ],
            GOEdgeType.CELLULAR_COMPONENT_TO_CELLULAR_COMPONENT: [
                GONodeType.CELLULAR_COMPONENT
            ],
            GOEdgeType.MOLECULAR_FUNCTION_TO_MOLECULAR_FUNCTION: [
                GONodeType.MOLECULAR_FUNCTION
            ],
            GOEdgeType.BIOLOGICAL_PROCESS_TO_MOLECULAR_FUNCTION: [
                GONodeType.BIOLOGICAL_PROCESS,
                GONodeType.MOLECULAR_FUNCTION,
            ],
        }
        # for checking edge labels
        self.domain_to_go_edge_types = [
            GOEdgeType.DOMAIN_TO_CELLULAR_COMPONENT,
            GOEdgeType.DOMAIN_TO_BIOLOGICAL_PROCESS,
            GOEdgeType.DOMAIN_TO_MOLECULAR_FUNCTION,
        ]
        self.protein_to_go_edge_types = [
            GOEdgeType.PROTEIN_TO_CELLULAR_COMPONENT,
            GOEdgeType.PROTEIN_TO_BIOLOGICAL_PROCESS,
            GOEdgeType.PROTEIN_TO_MOLECULAR_FUNCTION,
        ]
        self.go_to_go_edge_types = [
            GOEdgeType.CELLULAR_COMPONENT_TO_CELLULAR_COMPONENT,
            GOEdgeType.BIOLOGICAL_PROCESS_TO_BIOLOGICAL_PROCESS,
            GOEdgeType.MOLECULAR_FUNCTION_TO_MOLECULAR_FUNCTION,
            GOEdgeType.BIOLOGICAL_PROCESS_TO_MOLECULAR_FUNCTION,
        ]

        self.check_edge_type_duals = {
            GOEdgeType.PROTEIN_TO_CELLULAR_COMPONENT: "protein-cellular component",
            GOEdgeType.PROTEIN_TO_BIOLOGICAL_PROCESS: "protein-biological process",
            GOEdgeType.PROTEIN_TO_MOLECULAR_FUNCTION: "protein-molecular function",
            GOEdgeType.CELLULAR_COMPONENT_TO_CELLULAR_COMPONENT: "cellular component-cellular component",
            GOEdgeType.BIOLOGICAL_PROCESS_TO_BIOLOGICAL_PROCESS: "biological process-biological process",
            GOEdgeType.MOLECULAR_FUNCTION_TO_MOLECULAR_FUNCTION: "molecular function-molecular function",
            GOEdgeType.BIOLOGICAL_PROCESS_TO_MOLECULAR_FUNCTION: "biological process-molecular function",
            GOEdgeType.DOMAIN_TO_CELLULAR_COMPONENT: "domain-cellular component",
            GOEdgeType.DOMAIN_TO_BIOLOGICAL_PROCESS: "domain-biological process",
            GOEdgeType.DOMAIN_TO_MOLECULAR_FUNCTION: "domain-molecular function",
        }

        # will be used for edge filtering
        self.edge_filterer = set()

        # set node and edge types
        self.set_node_and_edge_types(
            node_types=model["node_types"], edge_types=model["edge_types"]
        )

        # set edge labels
        self.set_edge_labels(edge_labels=model["edge_labels"])

        # set node and edge properties
        self.set_node_and_edge_properties(
            go_node_fields=model["go_node_fields"],
            go_edge_fields=model["go_edge_fields"],
        )

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_go_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 6,
        all_go_annotations_url: HttpUrl = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gaf.gz",
        all_annotations_output_path: FilePath | None = None,
        anc2vec_embedding_path: FilePath | None = None
    ):
        """
        Wrapper function to download Gene Ontology data using pypath; used to access
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

            logger.debug("Started downloading Gene Ontology entry data")
            t0 = time()

            self.go_ontology = go_util.GeneOntology()

            t1 = time()
            logger.info(
                f"Gene Ontology entry data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )
            
            if GONodeField.ANC2VEC_EMBBEDDING.value in self.go_node_fields:
                self.retrieve_anc2vec_embedding(anc2vec_embedding_path)

            if any(
                [
                    True if et in self.edge_types else False
                    for et in self.protein_to_go_edge_types
                ]
            ):
                t0 = time()
                self.swissprots = set(
                    uniprot._all_uniprots(organism="*", swissprot=True)
                )

                if self.organism in ("*", None):
                    logger.debug(
                        "Started downloading Gene Ontology annotation data for all organisms"
                    )

                    if all_annotations_output_path:
                        full_path = all_annotations_output_path
                    else:
                        full_path = os.path.join(
                            os.getcwd(), "goa_uniprot_gcrp.gaf.gz"
                        )

                    if not os.path.isfile(full_path):
                        with requests.get(
                            all_go_annotations_url, stream=True
                        ) as response:
                            with open(full_path, "wb") as f:
                                for chunk in response.iter_content(1024):
                                    if chunk:
                                        f.write(chunk)

                    colnames = [
                        "DB",
                        "DB_Object_ID",
                        "DB_Object_Symbol",
                        "Qualifier",
                        "GO_ID",
                        "DB:Reference",
                        "Evidence Code",
                        "With (or) From",
                        "Aspect",
                        "DB_Object_Name",
                        "DB_Object_Synonym",
                        "DB_Object_Type",
                        "Taxon and Interacting taxon",
                        "Date",
                        "Assigned_By",
                        "Annotation_Extension",
                        "Gene_Product_Form_ID",
                    ]
                    
                    go_annots = pd.read_csv(
                        full_path,
                        skiprows=10,
                        sep="\t",
                        names=colnames,
                        usecols=[
                            "DB_Object_ID",
                            "Qualifier",
                            "GO_ID",
                            "DB:Reference",
                            "Evidence Code",
                        ],
                        header=None,
                        low_memory=False,
                        chunksize=10_000_000,
                    )

                    self.go_annots_df = pd.DataFrame()
                    for chunk in go_annots:
                        chunk.columns = [
                            "entry",
                            "qualifier",
                            "go_id",
                            "reference",
                            "evidence_code",
                        ]
                        filtered = chunk[chunk["entry"].isin(self.swissprots)]
                        self.go_annots_df = pd.concat(
                            [self.go_annots_df, filtered], ignore_index=True
                        )

                    if self.remove_selected_annotations:
                        self.go_annots_df = self.go_annots_df[~self.go_annots_df["evidence_code"].isin(self.remove_selected_annotations)]
                    
                    self.go_annots_df.drop_duplicates(
                        subset=["entry", "go_id"], ignore_index=True, inplace=True
                        )
                else:
                    logger.debug(
                        f"Started downloading Gene Ontology annotation data for tax id {self.organism}"
                    )

                    self.go_annots = go_input.go_annotations_all(
                        organism=self.organism,
                        fields=[
                            "qualifier",
                            "go_id",
                            "reference",
                            "evidence_code",
                        ],
                    )  # returns dict of uniprot ids as keys and go term annotations as values

                t1 = time()

                logger.info(
                    f"Gene Ontology annotation data is downloaded in {round((t1-t0) / 60, 2)} mins"
                )

            if any(
                [
                    True if et in self.edge_types else False
                    for et in self.domain_to_go_edge_types
                ]
            ):
                logger.debug("Started downloading Interpro2go data")
                t0 = time()

                self.interpro2go = interpro.interpro_xrefs(
                    db_type="go"
                )  # returns dict of interpro ids as keys and go term annotations as values
                t1 = time()
                logger.info(
                    f"Interpro2go data is downloaded in {round((t1-t0) / 60, 2)} mins"
                )
    def retrieve_anc2vec_embedding(self, anc2vec_embedding_path: FilePath | None = None):
        logger.info("Retrieving Anc2vec go term embeddings")

        self.go_term_to_anc2vec_embedding = {}
        with h5py.File(anc2vec_embedding_path, "r") as f:
            for go_term, embedding in tqdm(f.items(), total=len(f.keys())):
                self.go_term_to_anc2vec_embedding[go_term] = np.array(embedding).astype(np.float16)

    def set_node_and_edge_types(
        self, node_types: list, edge_types: list
    ) -> None:
        """
        Prepare node and edge types

        Warning: If you don't provide required node types of an edge type, it will give an error.
        """
        self.node_types = node_types or list(GONodeType)

        if edge_types:
            self.edge_types = edge_types

            for edge_type in edge_types:
                self.create_edge_filterer(edge_type)
                check_node_types = self.check_node_types_of_edges[edge_type]

                for nt in check_node_types:
                    if nt not in self.node_types:
                        logger.error(
                            f"{nt} must be included in node_types list"
                        )
                        raise ValueError(
                            f"{nt} must be included in node_types list"
                        )

        else:
            self.edge_types = list(GOEdgeType)

            for _type in GOEdgeType:
                self.create_edge_filterer(_type)

    def set_edge_labels(self, edge_labels: list) -> None:
        """
        Prepare edge labels

        Warning: If you don't provide required edge type of an edge label, it will give an error
        """

        # create predefined list for edge label enum classes
        protein_enum_classes = [
            ProteinToCellularComponentEdgeLabel,
            ProteinToBiologicalProcessEdgeLabel,
            ProteinToMolecularFunctionEdgeLabel,
        ]
        go_enum_classes = [
            BiologicalProcessToBiologicalProcessEdgeLabel,
            BiologicalProcessToMolecularFunctionEdgeLabel,
            MolecularFunctionToMolecularFunctionEdgeLabel,
            CellularComponentToCellularComponentEdgeLabel,
        ]
        domain_enum_classes = [
            DomainToCellularComponentEdgeLabel,
            DomainToBiologicalProcessEdgeLabel,
            DomainToMolecularFunctionEdgeLabel,
        ]

        if edge_labels:
            # create sets for edge label selection
            self.protein_to_go_edge_labels = set()
            self.go_to_go_edge_labels = set()
            self.domain_to_go_edge_labels = set()

            for edge_label in edge_labels:
                # define Protein-Go edge labels
                for enum_class in protein_enum_classes:
                    if edge_label in enum_class:
                        self.protein_to_go_edge_labels.add(edge_label.value)

                        if (
                            enum_class.neccessary_edge_type()
                            not in self.edge_types
                        ):
                            logger.error(
                                f"{enum_class.neccessary_edge_type()} must be included in edge_types list"
                            )
                            raise ValueError(
                                f"{enum_class.neccessary_edge_type()} must be included in edge_types list"
                            )

                # define Go-Go edge labels
                for enum_class in go_enum_classes:
                    if edge_label in enum_class:
                        self.go_to_go_edge_labels.add(edge_label.value)

                        if (
                            enum_class.neccessary_edge_type()
                            not in self.edge_types
                        ):
                            logger.error(
                                f"{enum_class.neccessary_edge_type()} must be included in edge_types list"
                            )
                            raise ValueError(
                                f"{enum_class.neccessary_edge_type()} must be included in edge_types list"
                            )

                # define Domain-Go edge labels
                for enum_class in domain_enum_classes:
                    if edge_label in enum_class:
                        self.domain_to_go_edge_labels.add(edge_label.value)

                        if (
                            enum_class.neccessary_edge_type()
                            not in self.edge_types
                        ):
                            logger.error(
                                f"{enum_class.neccessary_edge_type()} must be included in edge_types list"
                            )
                            raise ValueError(
                                f"{enum_class.neccessary_edge_type()} must be included in edge_types list"
                            )

        else:
            # create lists for edge label selection
            self.protein_to_go_edge_labels = []
            self.go_to_go_edge_labels = []
            self.domain_to_go_edge_labels = []

            # define Protein-Go edge labels
            for enum_class in protein_enum_classes:
                self.protein_to_go_edge_labels.extend(
                    [label.value for label in enum_class]
                )

            # define Go-Go edge labels
            for enum_class in go_enum_classes:
                self.go_to_go_edge_labels.extend(
                    [label.value for label in enum_class]
                )

            # define Domain-Go edge labels
            for enum_class in domain_enum_classes:
                self.domain_to_go_edge_labels.extend(
                    [label.value for label in enum_class]
                )

            self.protein_to_go_edge_labels = set(self.protein_to_go_edge_labels)
            self.go_to_go_edge_labels = set(self.go_to_go_edge_labels)
            self.domain_to_go_edge_labels = set(self.domain_to_go_edge_labels)

    def set_node_and_edge_properties(
        self, go_node_fields: list, go_edge_fields: list
    ) -> None:
        """
        Prepare node and edge properties
        """
        if go_node_fields:
            self.go_node_fields = [field.value for field in go_node_fields]
        else:
            self.go_node_fields = [field.value for field in GONodeField]

        if go_edge_fields:
            self.go_edge_fields = [field.value for field in go_edge_fields]
        else:
            self.go_edge_fields = [field.value for field in GOEdgeField]

    def create_aspect_to_node_label_dict(self) -> None:
        """
        Creates a dictionary for node label annotation
        """
        self.aspect_to_node_label_dict = {}

        if GONodeType.CELLULAR_COMPONENT in self.node_types:
            self.aspect_to_node_label_dict["C"] = "cellular component"
        if GONodeType.BIOLOGICAL_PROCESS in self.node_types:
            self.aspect_to_node_label_dict["P"] = "biological process"
        if GONodeType.MOLECULAR_FUNCTION in self.node_types:
            self.aspect_to_node_label_dict["F"] = "molecular function"

    def create_edge_filterer(self, edge_type) -> None:
        """
        Creates a dictionary for edge type filtering
        """
        if self.check_edge_type_duals.get(edge_type, None):
            self.edge_filterer.add(self.check_edge_type_duals[edge_type])

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

    def get_go_nodes(self) -> list[tuple]:
        """
        Prepare nodes ready to import into BioCypher
        """

        if not hasattr(self, "go_ontology"):
            self.download_go_data(cache=True)

        logger.info("Preparing GO nodes.")

        # create a list for nodes
        node_list = []

        # create node label dict
        self.create_aspect_to_node_label_dict()

        counter = 0
        for go_term in tqdm(
            self.go_ontology.name.keys()
        ):  # keys in self.go_ontology.name is current go ids in the ontology

            if label := self.aspect_to_node_label_dict.get(
                self.go_ontology.aspect[go_term], None
            ):
                go_id = self.add_prefix_to_id("go", go_term)

                node_props = {}
                if (
                    GONodeField.NAME.value in self.go_node_fields
                    or GONodeField.NAME in self.go_node_fields
                ):
                    node_props[GONodeField.NAME.value] = (
                        self.go_ontology.name[go_term]
                        .replace("'", "^")
                        .replace("|", "")
                    )
                
                if GONodeField.ANC2VEC_EMBBEDDING.value in self.go_node_fields and self.go_term_to_anc2vec_embedding.get(go_term) is not None:
                    node_props[GONodeField.ANC2VEC_EMBBEDDING.value] = [str(emb) for emb in self.go_term_to_anc2vec_embedding[go_term]]
                    

                node_list.append((go_id, label, node_props))

                counter += 1

            if self.early_stopping and counter >= self.early_stopping:
                break

        return node_list

    def get_go_edges(self) -> list[tuple]:
        """
        Prepare edges ready to import into BioCypher
        """
        if not hasattr(self, "go_annots_df") and not hasattr(self, "go_annots"):
            self.download_go_data(cache=True)

        # in case someone wants get only edges, run this function again
        self.create_aspect_to_node_label_dict()

        # create a list for edges
        edge_list = []

        # PROTEIN-GO EDGES
        if any(
            [
                True if et in self.edge_types else False
                for et in self.protein_to_go_edge_types
            ]
        ):
            logger.info("Preparing Protein-GO edges.")

            self.protein_to_go_edges = []

            if self.organism in ("*", None):
                for index, row in tqdm(
                    self.go_annots_df.iterrows(),
                    total=self.go_annots_df.shape[0],
                ):
                    if (
                        row["go_id"] in self.go_ontology.aspect.keys()
                        and str(row["qualifier"])
                        in self.protein_to_go_edge_labels
                        and (
                            self.aspect_to_node_label_dict.get(
                                self.go_ontology.aspect[row["go_id"]], None
                            )
                            and "protein-"
                            + self.aspect_to_node_label_dict[
                                self.go_ontology.aspect[row["go_id"]]
                            ]
                            in self.edge_filterer
                        )
                    ):
                        protein_id = self.add_prefix_to_id(
                            "uniprot", row["entry"]
                        )
                        go_id = self.add_prefix_to_id("go", row["go_id"])
                        edge_label = "_".join(
                            [
                                "protein",
                                str(row["qualifier"]).replace(" ", "_"),
                                self.aspect_to_node_label_dict[
                                    self.go_ontology.aspect[row["go_id"]]
                                ].replace(" ", "_"),
                            ]
                        )

                        props = {}

                        if (
                            GOEdgeField.REFERENCE in self.go_edge_fields
                            or GOEdgeField.REFERENCE.value
                            in self.go_edge_fields
                        ):
                            props[GOEdgeField.REFERENCE.value] = row[
                                "reference"
                            ]

                        if (
                            GOEdgeField.EVIDENCE_CODE in self.go_edge_fields
                            or GOEdgeField.EVIDENCE_CODE.value
                            in self.go_edge_fields
                        ):
                            props[GOEdgeField.EVIDENCE_CODE.value] = row[
                                "evidence_code"
                            ]

                        self.protein_to_go_edges.append(
                            (None, protein_id, go_id, edge_label, props)
                        )  # TODO delete this row after checking data
                        edge_list.append(
                            (None, protein_id, go_id, edge_label, props)
                        )

                    if self.early_stopping and index >= self.early_stopping:
                        break
            else:
                counter = 0
                for k, v in tqdm(self.go_annots.items()):
                    if k in self.swissprots:
                        protein_id = self.add_prefix_to_id("uniprot", k)
                        for annotation in list(v):
                            # subtract annotations and qualifiers that are not in self.protein_to_go_edge_labels and the ones that not in go ontology
                            if (
                                annotation.go_id
                                in self.go_ontology.aspect.keys()
                                and annotation.evidence_code
                                not in self.remove_selected_annotations
                                and str(annotation.qualifier)
                                in self.protein_to_go_edge_labels
                                and (
                                    self.aspect_to_node_label_dict.get(
                                        self.go_ontology.aspect[
                                            annotation.go_id
                                        ],
                                        None,
                                    )
                                    and "protein-"
                                    + self.aspect_to_node_label_dict[
                                        self.go_ontology.aspect[
                                            annotation.go_id
                                        ]
                                    ]
                                    in self.edge_filterer
                                )
                            ):
                                go_id = self.add_prefix_to_id(
                                    "go", annotation.go_id
                                )
                                edge_label = "_".join(
                                    [
                                        "protein",
                                        str(annotation.qualifier).replace(
                                            " ", "_"
                                        ),
                                        self.aspect_to_node_label_dict[
                                            self.go_ontology.aspect[
                                                annotation.go_id
                                            ]
                                        ].replace(" ", "_"),
                                    ]
                                )

                                props = {}
                                if (
                                    GOEdgeField.REFERENCE in self.go_edge_fields
                                    or GOEdgeField.REFERENCE.value
                                    in self.go_edge_fields
                                ):
                                    props[GOEdgeField.REFERENCE.value] = (
                                        annotation.reference
                                    )

                                if (
                                    GOEdgeField.EVIDENCE_CODE
                                    in self.go_edge_fields
                                    or GOEdgeField.EVIDENCE_CODE.value
                                    in self.go_edge_fields
                                ):
                                    props[GOEdgeField.EVIDENCE_CODE.value] = (
                                        annotation.evidence_code
                                    )

                                self.protein_to_go_edges.append(
                                    (None, protein_id, go_id, edge_label, props)
                                )  # TODO delete this row after checking data
                                edge_list.append(
                                    (None, protein_id, go_id, edge_label, props)
                                )

                                counter += 1

                    if self.early_stopping and counter >= self.early_stopping:
                        break

        # GO-GO EDGES
        if any(
            [
                True if et in self.edge_types else False
                for et in self.go_to_go_edge_types
            ]
        ):
            logger.info("Preparing GO-GO edges.")

            self.go_to_go_edges = []

            counter = 0
            for k, v in tqdm(self.go_ontology.ancestors.items()):
                source_go_id = self.add_prefix_to_id("go", k)

                for ancestor in list(v):
                    if (
                        str(ancestor[1]) in self.go_to_go_edge_labels
                        and str(k) in self.go_ontology.aspect.keys()
                        and ancestor[0] in self.go_ontology.aspect.keys()
                        and (
                            self.aspect_to_node_label_dict.get(
                                self.go_ontology.aspect[k], None
                            )
                            and self.aspect_to_node_label_dict.get(
                                self.go_ontology.aspect[ancestor[0]], None
                            )
                            and (
                                self.aspect_to_node_label_dict[
                                    self.go_ontology.aspect[k]
                                ]
                                + "-"
                                + self.aspect_to_node_label_dict[
                                    self.go_ontology.aspect[ancestor[0]]
                                ]
                            )
                            in self.edge_filterer
                        )
                    ):
                        target_go_id = self.add_prefix_to_id("go", ancestor[0])
                        edge_label = "_".join(
                            [
                                self.aspect_to_node_label_dict[
                                    self.go_ontology.aspect[k]
                                ].replace(" ", "_"),
                                ancestor[1],
                                self.aspect_to_node_label_dict[
                                    self.go_ontology.aspect[ancestor[0]]
                                ].replace(" ", "_"),
                            ]
                        )
                        self.go_to_go_edges.append(
                            (None, source_go_id, target_go_id, edge_label, {})
                        )  # TODO delete this row after checking data and keep only self.edge_list.append() line
                        edge_list.append(
                            (None, source_go_id, target_go_id, edge_label, {})
                        )

                        counter += 1

                if self.early_stopping and counter >= self.early_stopping:
                    break

        # DOMAIN-GO EDGES
        if any(
            [
                True if et in self.edge_types else False
                for et in self.domain_to_go_edge_types
            ]
        ):
            logger.info("Preparing Domain-GO edges.")

            domain_function_label_dict = {
                "P": "involved_in",
                "F": "enables",
                "C": "located_in",
            }

            self.domain_to_go_edges = []

            counter = 0
            for k, v in tqdm(self.interpro2go.items()):
                if v:
                    for go_term in v:
                        if go_term in self.go_ontology.aspect.keys():
                            aspect = self.go_ontology.aspect.get(go_term)
                            if (
                                domain_function_label_dict.get(aspect)
                                in self.domain_to_go_edge_labels
                                and self.aspect_to_node_label_dict.get(
                                    self.go_ontology.aspect[go_term], None
                                )
                                and "domain-"
                                + self.aspect_to_node_label_dict[
                                    self.go_ontology.aspect[go_term]
                                ]
                                in self.edge_filterer
                            ):

                                edge_label = "_".join(
                                    [
                                        "protein_domain",
                                        domain_function_label_dict.get(aspect),
                                        self.aspect_to_node_label_dict[
                                            self.go_ontology.aspect[go_term]
                                        ].replace(" ", "_"),
                                    ]
                                )
                                interpro_id = self.add_prefix_to_id(
                                    "interpro", k
                                )
                                go_id = self.add_prefix_to_id("go", go_term)
                                self.domain_to_go_edges.append(
                                    (None, interpro_id, go_id, edge_label, {})
                                )  # TODO delete this row after checking data and keep only self.edge_list.append() line
                                edge_list.append(
                                    (None, interpro_id, go_id, edge_label, {})
                                )

                                counter += 1

                if self.early_stopping and counter >= self.early_stopping:
                    break

        return edge_list

    def export_as_csv(self, path: str | None = None) -> None:
        # Write nodes
        nodes = self.get_go_nodes()

        nodes_df_dict = collections.defaultdict(list)
        for n in nodes:
            row = {"id": n[0]} | n[2]
            nodes_df_dict[n[1]].append(row)

        for label, data in nodes_df_dict.items():
            df = pd.DataFrame.from_records(data)

            if path:
                full_path = os.path.join(
                    path, f"{label.replace(' ','_').capitalize()}.csv"
                )
            else:
                full_path = f"{label.replace(' ','_').capitalize()}.csv"

            df.to_csv(full_path, index=False)
            logger.info(
                f"{label.replace(' ','_').capitalize()} data is written: {full_path}"
            )

        # Write edges
        if (
            not hasattr(self, "protein_to_go_edges")
            or not hasattr(self, "go_to_go_edges")
            or not hasattr(self, "domain_to_go_edges")
        ):
            _ = self.get_go_edges()

        edge_data_dict = {
            "protein_to_go": self.protein_to_go_edges,
            "go_to_go": self.go_to_go_edges,
            "domain_to_go": self.domain_to_go_edges,
        }
        edges_df_dict = collections.defaultdict(list)

        for label, edge_type in edge_data_dict.items():
            for e in edge_type:
                row = {"source": e[1], "target": e[2], "label": e[3]} | e[4]
                edges_df_dict[label].append(row)

        for label, data in edges_df_dict.items():
            df = pd.DataFrame.from_records(data)

            if path:
                full_path = os.path.join(path, f"{label.capitalize()}.csv")
            else:
                full_path = f"{label.capitalize()}.csv"

            df.to_csv(full_path, index=False)
            logger.info(f"{label.capitalize()} data is written: {full_path}")
