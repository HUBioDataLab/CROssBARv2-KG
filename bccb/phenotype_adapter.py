from __future__ import annotations

from pypath.share import curl, settings

from pypath.inputs import hpo, ontology

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time

import os
import h5py

from biocypher._logger import logger
from enum import Enum, EnumMeta, auto
from pydantic import BaseModel, DirectoryPath, FilePath, validate_call

from typing import Union

import pandas as pd
import numpy as np

logger.debug(f"Loading module {__name__}.")


class PhenotypeEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class PhenotypeNodeField(Enum, metaclass=PhenotypeEnumMeta):
    NAME = "name"
    SYNONYMS = "synonyms"

    # embedding
    CADA_EMBEDDING = "cada_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class PhenotypeDiseaseEdgeField(Enum, metaclass=PhenotypeEnumMeta):
    PUBMED_IDS = "pubmed_ids"
    EVIDENCE = "evidence"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class PhenotypeEdgeType(Enum, metaclass=PhenotypeEnumMeta):
    PROTEIN_TO_PHENOTYPE = auto()
    PHENOTYPE_HIERARCHICAL_EDGES = auto()
    PHENOTYPE_TO_DISEASE = auto()


class HPOModel(BaseModel):
    phenotype_node_fields: Union[list[PhenotypeNodeField], None] = None
    phenotype_disease_edge_fields: Union[
        list[PhenotypeDiseaseEdgeField], None
    ] = None
    edge_types: Union[list[PhenotypeEdgeType], None] = None
    remove_selected_annotations: list = ["IEA"]
    add_prefix: bool = True
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None


class HPO:
    def __init__(
        self,
        phenotype_node_fields: Union[list[PhenotypeNodeField], None] = None,
        phenotype_disease_edge_fields: Union[list[PhenotypeDiseaseEdgeField], None] = None,
        edge_types: Union[list[PhenotypeEdgeType], None] = None,
        remove_selected_annotations: list = ["IEA"],
        add_prefix: bool = True,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
    ):
        """
        Args:
            phenotype_node_fields: Phenotype node fields that will be included in graph, if defined it must be values of elements from PhenotypeNodeField enum class (not the names)
            phenotype_disease_edge_fields: Phenotype-Disease edge fields that will included in grah, if defined it must be values of elements from PhenotypeDiseaseEdgeField enum class (not the names)
            edge_types: list of edge types that will be included in graph, if defined it must be elements (not values of elements) from PhenotypeEdgeType enum class
            remove_selected_annotations: removes selected annotations from phenotype-disease edges, by default it removes electronic annotations
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of output data
            export_csv: if True, export data as csv
            output_dir: Location of csv export if `export_csv` is True, if not defined and `export_csv` is True, it will be current directory
        """

        model = HPOModel(
            phenotype_node_fields=phenotype_node_fields,
            phenotype_disease_edge_fields=phenotype_disease_edge_fields,
            edge_types=edge_types,
            remove_selected_annotations=remove_selected_annotations,
            add_prefix=add_prefix,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
        ).model_dump()

        self.add_prefix = model["add_prefix"]
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]
        self.remove_selected_annotations = model["remove_selected_annotations"]

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set node fields
        self.set_node_fields(
            phenotype_node_fields=model["phenotype_node_fields"]
        )

        # set edge fields
        self.set_edge_fields(
            phenotype_disease_edge_fields=model["phenotype_disease_edge_fields"]
        )

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_hpo_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
        cada_embedding_path: FilePath | None = None
    ) -> None:
        """
        Wrapper function to download hpo data from various databases using pypath.
        Args
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """

        with ExitStack() as stack:
            stack.enter_context(settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())

            logger.debug("Started downloading HPO data")
            t0 = time()

            if PhenotypeEdgeType.PROTEIN_TO_PHENOTYPE in self.edge_types:
                self.protein_hpo_annotations = hpo.hpo_annotations()

            self.hpo_ontology = hpo.hpo_ontology()

            self.hpo_terms = hpo.hpo_terms()

            if PhenotypeEdgeType.PHENOTYPE_TO_DISEASE in self.edge_types:
                self.hpo_phenotype_disease = hpo.hpo_diseases()

            if PhenotypeNodeField.CADA_EMBEDDING.value in self.phenotype_node_fields:
                self.retrieve_cada_embeddings(cada_embedding_path=cada_embedding_path)

            t1 = time()
            logger.info(
                f"HPO data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )
    @validate_call
    def retrieve_cada_embeddings(self, cada_embedding_path: FilePath | None = None) -> None:
        logger.info("Retrieving CADA phenotype embeddings.")

        self.hpo_id_to_cada_embedding = {}
        with h5py.File(cada_embedding_path, "r") as f:
            for hpo_id, embedding in f.items():
                self.hpo_id_to_cada_embedding[hpo_id] = np.array(embedding).astype(np.float16)

    def process_phenotype_disease(self) -> pd.DataFrame:

        if not hasattr(self, "hpo_phenotype_disease"):
            self.download_hpo_data()
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mondo_mappings()

        logger.debug(
            "Started processing HPO phenotype-disease interaction data"
        )
        t0 = time()

        df_list = []
        for hpo_id, diseases in self.hpo_phenotype_disease.items():
            if hpo_id == "hpo_id":
                continue

            for disease in diseases:
                if (
                    disease.evidence not in self.remove_selected_annotations
                    and disease.omim.split(":")[0] == "OMIM"
                    and self.mondo_mappings.get(disease.omim.split(":")[1])
                ):
                    if disease.pmid:
                        if ";" in disease.pmid:
                            pmid = "|".join(
                                [
                                    i.replace("PMID:", "")
                                    for i in disease.pmid.split(";")
                                ]
                            )
                        else:
                            pmid = disease.pmid.replace("PMID:", "")
                    else:
                        pmid = None

                    df_list.append(
                        (
                            hpo_id,
                            self.mondo_mappings.get(disease.omim.split(":")[1]),
                            pmid,
                            disease.evidence,
                        )
                    )

        df = pd.DataFrame(
            df_list, columns=["hpo_id", "disease_id", "pubmed_ids", "evidence"]
        )
        df.fillna(value=np.nan, inplace=True)

        df = df.groupby(
            ["hpo_id", "disease_id"], sort=False, as_index=False
        ).aggregate(
            {
                "hpo_id": "first",
                "disease_id": "first",
                "pubmed_ids": self.merge_source_column,
                "evidence": "first",
            }
        )
        df.replace("", np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"HPO phenotype-disease interaction data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        # write phenotype-disease edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Phenotype_to_disease.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Phenotype_to_disease.csv"
                )

            df.to_csv(full_path, index=False)
            logger.info(f"Phenotype-disease data is written: {full_path}")

        return df

    @validate_call
    def get_nodes(self, label: str = "phenotype") -> list[tuple]:

        if not hasattr(self, "hpo_terms"):
            self.download_hpo_data()

        logger.debug("Preparing phenotype nodes")

        node_list = []

        for index, (term, name) in tqdm(enumerate(self.hpo_terms.items())):
            hpo_id = self.add_prefix_to_id(prefix="hp", identifier=term)

            props = {}
            if PhenotypeNodeField.NAME.value in self.phenotype_node_fields:
                props[PhenotypeNodeField.NAME.value] = name.replace(
                    "|", ","
                ).replace("'", "^")

            if (
                PhenotypeNodeField.SYNONYMS.value in self.phenotype_node_fields
                and self.hpo_ontology["synonyms"].get(term)
            ):
                if len(self.hpo_ontology["synonyms"].get(term)) == 1:
                    props[PhenotypeNodeField.SYNONYMS.value] = (
                        list(self.hpo_ontology["synonyms"].get(term))[0]
                        .replace("|", ",")
                        .replace("'", "^")
                    )
                else:
                    props[PhenotypeNodeField.SYNONYMS.value] = [
                        t.replace("|", ",").replace("'", "^")
                        for t in self.hpo_ontology["synonyms"].get(term)
                    ]

            # get CADA phenotype embeddings
            if PhenotypeNodeField.CADA_EMBEDDING.value in self.phenotype_node_fields and self.hpo_id_to_cada_embedding.get(term):
                props[PhenotypeNodeField.CADA_EMBEDDING.value] = [str(emb) for emb in self.hpo_id_to_cada_embedding[term]]

            node_list.append((hpo_id, label, props))

            if self.early_stopping and index >= self.early_stopping:
                break

        # write node data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Phenotype.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Phenotype.csv")

            phenotype_df_list = [
                {"hpo_id": hpo_id} | props for hpo_id, _, props in node_list
            ]

            df = pd.DataFrame.from_records(phenotype_df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Phenotype node data is written: {full_path}")

        return node_list

    @validate_call
    def get_edges(self,
                  protein_phenotype_label: str = "protein_is_associated_with_phenotype",
                  phenotype_hierarchical_label: str = "phenotype_is_a_phenotype",
                  phenotype_disease_label: str = "phenotype_is_associated_with_disease") -> list[tuple]:

        logger.info("Preparing all edge types")

        edge_list = []

        if PhenotypeEdgeType.PROTEIN_TO_PHENOTYPE in self.edge_types:
            edge_list.extend(self.get_protein_phenotype_edges(protein_phenotype_label))

        if PhenotypeEdgeType.PHENOTYPE_HIERARCHICAL_EDGES in self.edge_types:
            edge_list.extend(self.get_phenotype_hierarchical_edges(phenotype_hierarchical_label))

        if PhenotypeEdgeType.PHENOTYPE_TO_DISEASE in self.edge_types:
            edge_list.extend(self.get_phenotype_disease_edges(phenotype_disease_label))

        return edge_list

    @validate_call
    def get_protein_phenotype_edges(
        self, label: str = "protein_is_associated_with_phenotype"
    ) -> list[tuple]:
        if not hasattr(self, "protein_hpo_annotations"):
            self.download_hpo_data()

        logger.debug("Preparing protein-phenotype edges")

        edge_list = set()

        for index, (uniprot_id, annotations) in tqdm(
            enumerate(self.protein_hpo_annotations.items())
        ):
            protein_id = self.add_prefix_to_id(
                prefix="uniprot", identifier=uniprot_id
            )
            for annot in annotations:
                hpo_id = self.add_prefix_to_id(
                    prefix="hp", identifier=annot.hpo_id
                )
                edge_list.add((None, protein_id, hpo_id, label))

            if self.early_stopping and index >= self.early_stopping:
                break

        # write protein-phenotype edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Protein_to_phenotype.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Protein_to_phenotype.csv"
                )

            df_list = [
                {"protein_id": protein_id, "hpo_id": hpo_id}
                for _, protein_id, hpo_id, _ in edge_list
            ]

            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Protein-phenotype edge data is written: {full_path}")

        return [i + ({},) for i in edge_list]

    @validate_call
    def get_phenotype_hierarchical_edges(
        self, label: str = "phenotype_is_a_phenotype"
    ) -> list[tuple]:

        if not hasattr(self, "hpo_ontology"):
            self.download_hpo_data()

        logger.debug("Preparing phenotype hierarchical edges")

        edge_list = []

        for index, (child, parents) in tqdm(
            enumerate(self.hpo_ontology["parents"].items())
        ):
            child_id = self.add_prefix_to_id(prefix="hp", identifier=child)
            for parent in parents:
                parent_id = self.add_prefix_to_id(
                    prefix="hp", identifier=parent
                )
                edge_list.append((None, child_id, parent_id, label, {}))

            if self.early_stopping and index >= self.early_stopping:
                break

        # write phenotype hierarchical edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Phenotype_hierarchical_edges.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Phenotype_hierarchical_edges.csv"
                )

            df_list = [
                {"child_id": child_id, "parent_id": parent_id}
                for _, child_id, parent_id, _, _ in edge_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(
                f"Phenotype hierarchical edge data is written: {full_path}"
            )

        return edge_list

    @validate_call
    def get_phenotype_disease_edges(
        self, label: str = "phenotype_is_associated_with_disease"
    ) -> list[tuple]:

        phenotype_disease_df = self.process_phenotype_disease()

        logger.debug("Preparing phenotype-disease edges")

        edge_list = []

        for index, row in tqdm(
            phenotype_disease_df.iterrows(), total=phenotype_disease_df.shape[0]
        ):
            _dict = row.to_dict()

            hpo_id = self.add_prefix_to_id(
                prefix="hp", identifier=_dict["hpo_id"]
            )
            disease_id = self.add_prefix_to_id(
                prefix="MONDO", identifier=_dict["disease_id"]
            )

            del _dict["disease_id"], _dict["hpo_id"]

            props = {}
            for k, v in _dict.items():
                if k in self.phenotype_disease_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v

            edge_list.append((None, hpo_id, disease_id, label, props))

            if self.early_stopping and index >= self.early_stopping:
                break

        return edge_list

    def prepare_mondo_mappings(self):

        logger.debug("Preparing mondo mappings to OMIM")

        mondo = ontology.ontology(
            ontology="mondo", fields=["is_obsolete", "obo_xref"]
        )

        self.mondo_mappings = {}

        mapping_db_list = ["OMIM"]

        for term in mondo:
            if (
                not term.is_obsolete
                and term.obo_id
                and "MONDO" in term.obo_id
                and term.obo_xref
            ):
                for xref in term.obo_xref:
                    if xref.get("database") in mapping_db_list:
                        self.mondo_mappings[xref["id"]] = term.obo_id

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

    def merge_source_column(self, element, joiner="|"):
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                _list.extend(iter(e.split(joiner)))
            else:
                _list.append(e)

        return joiner.join(list(dict.fromkeys(_list).keys()))

    def set_node_fields(self, phenotype_node_fields):
        if phenotype_node_fields:
            self.phenotype_node_fields = [
                field.value for field in phenotype_node_fields
            ]
        else:
            self.phenotype_node_fields = [
                field.value for field in PhenotypeNodeField
            ]

    def set_edge_fields(self, phenotype_disease_edge_fields):
        if phenotype_disease_edge_fields:
            self.phenotype_disease_edge_fields = [
                field.value for field in phenotype_disease_edge_fields
            ]
        else:
            self.phenotype_disease_edge_fields = [
                field.value for field in PhenotypeDiseaseEdgeField
            ]

    def set_edge_types(self, edge_types):
        self.edge_types = edge_types or list(PhenotypeEdgeType)
