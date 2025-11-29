from __future__ import annotations
from pypath.share import curl, settings

from pypath.inputs import reactome, uniprot, ctdbase, compath, unichem, drugbank
from pypath.inputs import ontology
from . import kegg_local

from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
from biocypher._logger import logger

import collections
import os
import h5py
import gzip

from pydantic import BaseModel, DirectoryPath, FilePath, EmailStr, validate_call
from typing import Union

from enum import Enum, EnumMeta, auto

import pandas as pd
import numpy as np


class PathwayEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class PathwayNodeField(Enum, metaclass=PathwayEnumMeta):
    NAME = "name"
    ORGANISM = "organism"
    BIOKEEN_EMBEDDING = "biokeen_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class ProteinPathwayEdgeField(Enum, metaclass=PathwayEnumMeta):
    SOURCE = "source"
    EVIDENCE_CODE = "evidence_code"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None
    
class DiseasePathwayEdgeField(Enum, metaclass=PathwayEnumMeta):
    SOURCE = "source"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None
    
class DrugPathwayEdgeField(Enum, metaclass=PathwayEnumMeta):
    SOURCE = "source"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None



class PathwayEdgeType(Enum, metaclass=PathwayEnumMeta):
    PROTEIN_TO_PATHWAY = auto()
    REACTOME_HIERARCHICAL_RELATIONS = auto()
    DRUG_TO_PATHWAY = auto()
    DISEASE_TO_PATHWAY = auto()
    PATHWAY_TO_PATHWAY = auto()
    PATHWAY_ORTHOLOGY = auto()

class PathwayNodeType(Enum, metaclass=PathwayEnumMeta):
    PATHWAY = auto()


logger.debug(f"Loading module {__name__}.")


class PathwayModel(BaseModel):
    drugbank_user: EmailStr
    drugbank_passwd: str
    pathway_node_fields: Union[list[PathwayNodeField], None] = None
    protein_pathway_edge_fields: Union[list[ProteinPathwayEdgeField], None] = None
    disease_pathway_edge_fields: Union[list[DiseasePathwayEdgeField], None] = None
    drug_pathway_edge_fields: Union[list[DrugPathwayEdgeField], None] = None
    node_types: Union[list[PathwayNodeType], None] = None
    edge_types: Union[list[PathwayEdgeType], None] = None
    remove_selected_annotations: list[str] = ["IEA"]
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True
    kegg_organism: list[str] | str | None = None


# ADD evidence_code to schema
class Pathway:
    def __init__(
        self,
        drugbank_user: EmailStr,
        drugbank_passwd: str,
        pathway_node_fields: Union[list[PathwayNodeField], None] = None,
        protein_pathway_edge_fields: Union[list[ProteinPathwayEdgeField], None] = None,
        disease_pathway_edge_fields: Union[list[DiseasePathwayEdgeField], None] = None,
        drug_pathway_edge_fields: Union[list[DrugPathwayEdgeField], None] = None,
        node_types: Union[list[PathwayNodeType], None] = None,
        edge_types: Union[list[PathwayEdgeType], None] = None,
        remove_selected_annotations: list[str] = ["IEA"],
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
        add_prefix: bool = True,
        kegg_organism: list[str] | str | None = None,
    ):
        """
        Args:
            drugbank_user: drugbank username
            drugbank_passwd: drugbank password
            pathway_node_fields: Pathway node fields that will be included in graph, if defined it must be values of elements from PathwayNodeField enum class (not the names)
            protein_pathway_edge_fields: Protein-pathway edge fields that will included in grah, if defined it must be values of elements from ProteinPathwayEdgeField enum class (not the names)
            edge_types: list of edge types that will be included in graph, if defined it must be elements (not values of elements) from PathwayEdgeType enum class
            remove_selected_annotations: removes selected annotations from phenotype-disease edges, by default it removes electronic annotations
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of output data
            export_csv: if True, export data as csv
            output_dir: Location of csv export if `export_csv` is True, if not defined and `export_csv` is True, it will be current directory
            kegg_organism: which kegg organisms' pathway data that will be included in graph, by default it include all the kegg organisms
            if defined, it should kegg organism prefixes.
        """

        model = PathwayModel(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            pathway_node_fields=pathway_node_fields,
            protein_pathway_edge_fields=protein_pathway_edge_fields,
            disease_pathway_edge_fields=disease_pathway_edge_fields,
            drug_pathway_edge_fields=drug_pathway_edge_fields,
            node_types=node_types,
            edge_types=edge_types,
            remove_selected_annotations=remove_selected_annotations,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
            add_prefix=add_prefix,
            kegg_organism=kegg_organism,
        ).model_dump()

        self.drugbank_user = model["drugbank_user"]
        self.drugbank_passwd = model["drugbank_passwd"]
        self.add_prefix = model["add_prefix"]
        self.remove_selected_annotations = model["remove_selected_annotations"]
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]

        # set kegg organisms list
        if not model["kegg_organism"]:
            self.kegg_organism = list(kegg_local._Organism()._data.keys())
        else:
            self.kegg_organism = self.ensure_iterable(model["kegg_organism"])

        # set node fields
        self.set_node_fields(pathway_node_fields=model["pathway_node_fields"])

        # set edge fields
        self.set_edge_fields(
            protein_pathway_edge_fields=model["protein_pathway_edge_fields"],
            disease_pathway_edge_fields=model["disease_pathway_edge_fields"],
            drug_pathway_edge_fields=model["drug_pathway_edge_fields"],
        )

        # set edge types
        self.set_node_and_edge_types(node_types=model["node_types"],
                                     edge_types=model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_pathway_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
        biokeen_embedding_path: FilePath | None = None,
    ) -> None:
        """
        Wrapper function to download pathway data from various databases using pypath.
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

            self.download_reactome_data()

            self.download_kegg_data()

            self.download_ctd_data()

            self.download_compath_data()

            if PathwayNodeField.BIOKEEN_EMBEDDING.value in self.pathway_node_fields:
                self.retrieve_biokeen_embeddings(biokeen_embedding_path=biokeen_embedding_path)

    def download_reactome_data(self) -> None:

        logger.debug("Started downloading Reactome data")
        t0 = time()

        if PathwayNodeType.PATHWAY in self.node_types:
            self.reactome_pathways = list(reactome.reactome_pathways())

        if PathwayEdgeType.REACTOME_HIERARCHICAL_RELATIONS in self.edge_types:
            self.reactome_hierarchial_relations = (
                reactome.reactome_pathway_relations()
            )

        if PathwayEdgeType.PROTEIN_TO_PATHWAY in self.edge_types:
            self.reactome_uniprot_pathway = reactome.reactome_uniprots()

        if PathwayEdgeType.DRUG_TO_PATHWAY in self.edge_types:
            self.reactome_chebi_pathway = reactome.reactome_chebis()
            self.chebi_to_drugbank = {
                list(v)[0]: k
                for k, v in unichem.unichem_mapping("drugbank", "chebi").items()
            }

        t1 = time()
        logger.info(
            f"Reactome data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def download_kegg_data(self) -> None:

        logger.debug("Started downloading KEGG data")
        t0 = time()

        if PathwayNodeType.PATHWAY in self.node_types:

            self.kegg_pathway_abbv_organism_name_dict = {
                k: v[1] for k, v in kegg_local._Organism()._data.items()
            }

            self.kegg_pathways = []
            for org in tqdm(self.kegg_organism):
                try:
                    self.kegg_pathways.extend(
                        kegg_local._kegg_list("pathway", org=org)
                    )
                except (IndexError, UnicodeDecodeError, gzip.BadGzipFile) as e:
                    logger.debug(
                        f"Error occured in {org} organism in pathway data downloading with an {e} error"
                    )

        if PathwayEdgeType.PROTEIN_TO_PATHWAY in self.edge_types:
            self.kegg_gene_to_pathway = {}
            for org in tqdm(self.kegg_organism):
                try:
                    self.kegg_gene_to_pathway = (
                        self.kegg_gene_to_pathway
                        | kegg_local.gene_to_pathway(org=org)
                    )
                except (IndexError, UnicodeDecodeError, gzip.BadGzipFile) as e:
                    logger.debug(
                        f"Error occured in {org} organism  in gene-pathway data downloading with an {e} error"
                    )

            self.kegg_to_uniprot = {
                v.strip(";").split(";")[0]: k
                for k, v in uniprot.uniprot_data(
                    "xref_kegg", 9606, True
                ).items()
            }

        if PathwayEdgeType.DRUG_TO_PATHWAY in self.edge_types:
            self.kegg_drug_to_pathway = kegg_local.drug_to_pathway()
            drugbank_data = drugbank.DrugbankFull(
                user=self.drugbank_user, passwd=self.drugbank_passwd
            )
            drugbank_drugs_external_ids = (
                drugbank_data.drugbank_external_ids_full()
            )
            self.kegg_drug_to_drugbank = {
                v.get("KEGG Drug"): k
                for k, v in drugbank_drugs_external_ids.items()
                if v.get("KEGG Drug")
            }

        if PathwayEdgeType.DISEASE_TO_PATHWAY in self.edge_types:

            self.kegg_disease_to_pathway = kegg_local.disease_to_pathway()

            self.kegg_diseases_mappings = {}
            for dis in self.kegg_disease_to_pathway.keys():
                result = kegg_local.get_diseases(dis)
                self.kegg_diseases_mappings[dis] = result[0].db_links

        t1 = time()
        logger.info(f"KEGG data is downloaded in {round((t1-t0) / 60, 2)} mins")

    def download_ctd_data(self) -> None:

        logger.debug("Started downloading CTD data")
        t0 = time()

        if PathwayEdgeType.DISEASE_TO_PATHWAY in self.edge_types:
            self.ctd_disease_pathway = ctdbase.ctdbase_relations(
                "disease_pathway"
            )

        t1 = time()
        logger.info(f"CTD data is downloaded in {round((t1-t0) / 60, 2)} mins")

    def download_compath_data(self) -> None:

        logger.debug("Started downloading Compath data")
        t0 = time()

        if PathwayEdgeType.PATHWAY_TO_PATHWAY in self.edge_types:
            self.compath_pathway_pathway = compath.compath_mappings(
                source_db=None, target_db=None
            )

        t1 = time()
        logger.info(
            f"Compath data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )
    @validate_call
    def retrieve_biokeen_embeddings(self, biokeen_embedding_path: FilePath) -> None:
        logger.info("Retrieving BioKenn pathway embeddings.")

        self.pathway_id_to_biokeen_embedding = {}
        with h5py.File(biokeen_embedding_path, "r") as f:
            for patway_id, embedding in f.items():
                self.pathway_id_to_biokeen_embedding[patway_id] = np.array(embedding).astype(np.float16)
    def process_reactome_protein_pathway(self) -> pd.DataFrame:

        if not hasattr(self, "reactome_uniprot_pathway"):
            self.download_reactome_data()

        logger.debug("Started processing Reactome protein-pathway data")
        t0 = time()

        df_list = [
            (pp.uniprot_id, pp.pathway_id, pp.evidence_code)
            for pp in self.reactome_uniprot_pathway
            if pp.evidence_code not in self.remove_selected_annotations
        ]
        df = pd.DataFrame(
            df_list, columns=["uniprot_id", "pathway_id", "evidence_code"]
        )

        df.drop_duplicates(ignore_index=True, inplace=True)

        df["source"] = "Reactome"

        t1 = time()
        logger.info(
            f"Reactome protein-pathway data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_kegg_protein_pathway(self) -> pd.DataFrame:

        if not hasattr(self, "kegg_gene_to_pathway"):
            self.download_kegg_data()

        logger.debug("Started processing KEGG protein-pathway data")
        t0 = time()

        df_list = []
        for gene, pathways in self.kegg_gene_to_pathway.items():
            if not gene.startswith("org") and self.kegg_to_uniprot.get(gene):
                organism_prefix = gene.split(":")[0].strip()
                for pathway in pathways.PathwayEntries:
                    df_list.append(
                        (
                            self.kegg_to_uniprot[gene],
                            pathway.pathway_id.replace("map", organism_prefix),
                        )
                    )

        df = pd.DataFrame(df_list, columns=["uniprot_id", "pathway_id"])

        df.drop_duplicates(ignore_index=True, inplace=True)

        df["source"] = "KEGG"

        t1 = time()
        logger.info(
            f"KEGG protein-pathway data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_reactome_drug_pathway(self) -> pd.DataFrame:

        if not hasattr(self, "reactome_chebi_pathway"):
            self.download_reactome_data()

        logger.debug("Started processing Reactome drug-pathway data")
        t0 = time()

        df_list = [
            (self.chebi_to_drugbank[cp.chebi_id], cp.pathway_id)
            for cp in self.reactome_chebi_pathway
            if cp.evidence_code not in self.remove_selected_annotations
            and self.chebi_to_drugbank.get(cp.chebi_id)
        ]
        df = pd.DataFrame(df_list, columns=["drug_id", "pathway_id"])

        df.drop_duplicates(ignore_index=True, inplace=True)

        df["source"] = "Reactome"

        t1 = time()
        logger.info(
            f"Reactome drug-pathway data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_kegg_drug_pathway(self) -> pd.DataFrame:

        if not hasattr(self, "kegg_drug_to_pathway"):
            self.download_kegg_data()

        logger.debug("Started processing KEGG drug-pathway data")
        t0 = time()

        df_list = []
        for drug, pathways in self.kegg_drug_to_pathway.items():
            if self.kegg_drug_to_drugbank.get(drug):
                for pathway in pathways.PathwayEntries:
                    df_list.append(
                        (
                            self.kegg_drug_to_drugbank[drug],
                            pathway.pathway_id.replace("map", "hsa"),
                        )
                    )

        df = pd.DataFrame(df_list, columns=["drug_id", "pathway_id"])

        df.drop_duplicates(ignore_index=True, inplace=True)

        df["source"] = "KEGG"

        t1 = time()
        logger.info(
            f"KEGG drug-pathway data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_kegg_disease_pathway(self) -> pd.DataFrame:

        if not hasattr(self, "kegg_disease_to_pathway"):
            self.download_kegg_data()
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mondo_mappings()

        logger.debug("Started processing KEGG disease-pathway data")
        t0 = time()

        kegg_dbs_to_mondo_dbs = {
            "MeSH": "MESH",
            "OMIM": "OMIM",
            "ICD-10": "ICD10CM",
        }

        df_list = []
        for disease, pathways in self.kegg_disease_to_pathway.items():
            found = False
            disease_id = None
            for db in kegg_dbs_to_mondo_dbs.keys():
                if found:
                    break
                if self.kegg_diseases_mappings[disease].get(db):
                    for ref in self.ensure_iterable(
                        self.kegg_diseases_mappings[disease][db]
                    ):
                        if found:
                            break

                        if self.mondo_mappings[kegg_dbs_to_mondo_dbs[db]].get(
                            ref
                        ):
                            disease_id = self.mondo_mappings[
                                kegg_dbs_to_mondo_dbs[db]
                            ].get(ref)

            if disease_id:
                for pathway in pathways.PathwayEntries:
                    df_list.append(
                        (
                            disease_id,
                            pathway.pathway_id.replace("map", "hsa"),
                        )
                    )

        df = pd.DataFrame(df_list, columns=["disease_id", "pathway_id"])

        df.drop_duplicates(ignore_index=True, inplace=True)

        df["source"] = "KEGG"

        t1 = time()
        logger.info(
            f"KEGG disease-pathway data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_ctd_disease_pathway(self) -> pd.DataFrame:

        if not hasattr(self, "ctd_disease_pathway"):
            self.download_ctd_data()
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mondo_mappings()

        kegg_pathways_checker_list = [
            i[0] for i in kegg_local._kegg_list("pathway", org="hsa")
        ]

        logger.debug("Started processing CTD disease-pathway data")
        t0 = time()

        df_list = []
        for dp in self.ctd_disease_pathway:
            disease_db = dp.DiseaseID.split(":")[0]
            pathway_db = dp.PathwayID.split(":")[0]
            if self.mondo_mappings[disease_db].get(dp.DiseaseID.split(":")[1]):
                if pathway_db == "KEGG":
                    if (
                        dp.PathwayID.split(":")[1].replace("_", "")
                        in kegg_pathways_checker_list
                    ):
                        df_list.append(
                            (
                                self.mondo_mappings[disease_db].get(
                                    dp.DiseaseID.split(":")[1]
                                ),
                                dp.PathwayID.split(":")[1].replace("_", ""),
                            )
                        )
                else:
                    df_list.append(
                        (
                            self.mondo_mappings[disease_db].get(
                                dp.DiseaseID.split(":")[1]
                            ),
                            dp.PathwayID.split(":")[1],
                        )
                    )

        df = pd.DataFrame(df_list, columns=["disease_id", "pathway_id"])

        df.drop_duplicates(ignore_index=True, inplace=True)

        df["source"] = "CTD"

        t1 = time()
        logger.info(
            f"CTD disease-pathway data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def merge_protein_pathway_data(self) -> pd.DataFrame:

        kegg_df = self.process_kegg_protein_pathway()

        reactome_df = self.process_reactome_protein_pathway()

        logger.debug("Started merging protein-pathway data")
        t0 = time()

        merged_df = pd.concat([kegg_df, reactome_df], ignore_index=True)

        t1 = time()
        logger.info(
            f"Protein-pathway edge data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        # write protein-pathway edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Protein_to_pathway.csv"
                )
            else:
                full_path = os.path.join(os.getcwd(), "Protein_to_pathway.csv")

            merged_df.to_csv(full_path, index=False)
            logger.info(
                f"Protein-pathway edge data data is written: {full_path}"
            )

        return merged_df

    def merge_drug_pathway_data(self) -> pd.DataFrame:

        kegg_df = self.process_kegg_drug_pathway()

        reactome_df = self.process_reactome_drug_pathway()

        logger.debug("Started merging drug-pathway data")
        t0 = time()

        merged_df = pd.concat([kegg_df, reactome_df], ignore_index=True)

        t1 = time()
        logger.info(
            f"Drug-pathway edge data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        # write drug-pathway edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Drug_to_pathway.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Drug_to_pathway.csv")

            merged_df.to_csv(full_path, index=False)
            logger.info(f"Drug-pathway edge data data is written: {full_path}")

        return merged_df

    def merge_disease_pathway_data(self) -> pd.DataFrame:

        kegg_df = self.process_kegg_disease_pathway()

        ctd_df = self.process_ctd_disease_pathway()

        logger.debug("Started merging disease-pathway edge data")
        t0 = time()

        merged_df = pd.merge(
            kegg_df, ctd_df, how="outer", on=["disease_id", "pathway_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        t1 = time()
        logger.info(
            f"Disease-pathway edge data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        # write disease-pathway edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Disease_to_pathway.csv"
                )
            else:
                full_path = os.path.join(os.getcwd(), "Disease_to_pathway.csv")

            merged_df.to_csv(full_path, index=False)
            logger.info(
                f"Disease-pathway edge data data is written: {full_path}"
            )

        return merged_df

    @validate_call
    def get_nodes(self, label: str = "pathway") -> list[tuple]:

        if not hasattr(self, "reactome_pathways"):
            self.download_reactome_data()
        if not hasattr(self, "kegg_pathways"):
            self.download_kegg_data()

        logger.info("Started writing pathway nodes")

        node_list = []

        for index, p in tqdm(enumerate(self.reactome_pathways)):
            pathway_id = self.add_prefix_to_id(
                prefix="reactome", identifier=p.pathway_id
            )

            props = {}
            if PathwayNodeField.NAME.value in self.pathway_node_fields:
                props[PathwayNodeField.NAME.value] = p.pathway_name.replace(
                    "'", "^"
                )

            if PathwayNodeField.ORGANISM.value in self.pathway_node_fields:
                props[PathwayNodeField.ORGANISM.value] = p.organism.replace("'","^") if p.organism else None

            if PathwayNodeField.BIOKEEN_EMBEDDING.value in self.pathway_node_fields and self.pathway_id_to_biokeen_embedding.get(p.pathway_id) is not None:
                props[PathwayNodeField.BIOKEEN_EMBEDDING.value] = [str(emb) for emb in self.pathway_id_to_biokeen_embedding[p.pathway_id]]

            node_list.append((pathway_id, label, props))

            if self.early_stopping and index >= self.early_stopping:
                break

        for index, p in tqdm(enumerate(self.kegg_pathways)):
            pathway_id = self.add_prefix_to_id(
                prefix="kegg.pathway", identifier=p[0]
            )

            props = {}
            if PathwayNodeField.NAME.value in self.pathway_node_fields:
                props[PathwayNodeField.NAME.value] = (
                    p[1].split("-")[0].strip().replace("'", "^")
                )

            if PathwayNodeField.ORGANISM.value in self.pathway_node_fields:
                props[PathwayNodeField.ORGANISM.value] = (
                    self.kegg_pathway_abbv_organism_name_dict.get(p[0][:3]).replace("'", "^") if self.kegg_pathway_abbv_organism_name_dict.get(p[0][:3]) else None
                )
            
            if PathwayNodeField.BIOKEEN_EMBEDDING.value in self.pathway_node_fields and self.pathway_id_to_biokeen_embedding.get(p[0]) is not None:
                props[PathwayNodeField.BIOKEEN_EMBEDDING.value] = [str(emb) for emb in self.pathway_id_to_biokeen_embedding[p[0]]]

            node_list.append((pathway_id, label, props))

            if self.early_stopping and index >= self.early_stopping:
                break

        # write pathway node data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Pathway.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Pathway.csv")

            df_list = [
                {"pathway_id": pathway_id} | props
                for pathway_id, _, props in node_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Pathway node data is written: {full_path}")

        return node_list

    @validate_call
    def get_edges(self,
                  protein_pathway_label: str = "protein_take_part_in_pathway",
                  drug_pathway_label: str = "drug_has_target_in_pathway",
                  disease_pathway_label: str = "disease_modulates_pathway",
                  reactome_hierarchy_label: str = "pathway_participates_pathway",
                  pathway_orthology_label: str = "pathway_is_ortholog_to_pathway") -> list[tuple]:

        logger.info("Started writing all pathway edges")

        edge_list = []

        if PathwayEdgeType.PROTEIN_TO_PATHWAY in self.edge_types:
            edge_list.extend(self.get_protein_pathway_edges(protein_pathway_label))

        if PathwayEdgeType.DRUG_TO_PATHWAY in self.edge_types:
            edge_list.extend(self.get_drug_pathway_edges(drug_pathway_label))

        if PathwayEdgeType.DISEASE_TO_PATHWAY in self.edge_types:
            edge_list.extend(self.get_disease_pathway_edges(disease_pathway_label))

        if PathwayEdgeType.PATHWAY_TO_PATHWAY in self.edge_types:
            edge_list.extend(self.get_pathway_pathway_edges())

        if PathwayEdgeType.REACTOME_HIERARCHICAL_RELATIONS in self.edge_types:
            edge_list.extend(self.get_reactome_hierarchical_edges(reactome_hierarchy_label))

        if PathwayEdgeType.PATHWAY_ORTHOLOGY in self.edge_types:
            edge_list.extend(self.get_pathway_pathway_orthology_edges(pathway_orthology_label))

        return edge_list

    @validate_call
    def get_protein_pathway_edges(
        self, label: str = "protein_take_part_in_pathway"
    ) -> list[tuple]:

        protein_pathway_edges_df = self.merge_protein_pathway_data()

        logger.info("Started writing protein-pathway edges")

        edge_list = []
        for index, row in tqdm(
            protein_pathway_edges_df.iterrows(),
            total=protein_pathway_edges_df.shape[0],
        ):
            _dict = row.to_dict()

            if _dict["source"] == "Reactome":
                pathway_id = self.add_prefix_to_id(
                    prefix="reactome", identifier=_dict["pathway_id"]
                )
            else:
                pathway_id = self.add_prefix_to_id(
                    prefix="kegg.pathway", identifier=_dict["pathway_id"]
                )

            uniprot_id = self.add_prefix_to_id(
                prefix="uniprot", identifier=_dict["uniprot_id"]
            )

            del _dict["uniprot_id"], _dict["pathway_id"]

            props = {}
            for k, v in _dict.items():
                if k in self.protein_pathway_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v

            edge_list.append((None, uniprot_id, pathway_id, label, props))

            if self.early_stopping and index >= self.early_stopping:
                break

        return edge_list

    @validate_call
    def get_drug_pathway_edges(
        self, label: str = "drug_has_target_in_pathway"
    ) -> list[tuple]:

        drug_pathway_edges_df = self.merge_drug_pathway_data()

        logger.info("Started writing drug-pathway edges")

        edge_list = []
        for index, row in tqdm(
            drug_pathway_edges_df.iterrows(),
            total=drug_pathway_edges_df.shape[0],
        ):
            _dict = row.to_dict()

            if _dict["source"] == "Reactome":
                pathway_id = self.add_prefix_to_id(
                    prefix="reactome", identifier=_dict["pathway_id"]
                )
            else:
                pathway_id = self.add_prefix_to_id(
                    prefix="kegg.pathway", identifier=_dict["pathway_id"]
                )

            drug_id = self.add_prefix_to_id(
                prefix="drugbank", identifier=_dict["drug_id"]
            )

            del _dict["drug_id"], _dict["pathway_id"]

            props = {}
            for k, v in _dict.items():
                if k in self.drug_pathway_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v

            edge_list.append((None, drug_id, pathway_id, label, props))

            if self.early_stopping and index >= self.early_stopping:
                break

        return edge_list

    @validate_call
    def get_disease_pathway_edges(
        self, label: str = "disease_modulates_pathway"
    ) -> list[tuple]:

        disease_pathway_edges_df = self.merge_disease_pathway_data()

        logger.info("Started writing disease-pathway edges")

        edge_list = []
        for index, row in tqdm(
            disease_pathway_edges_df.iterrows(),
            total=disease_pathway_edges_df.shape[0],
        ):
            _dict = row.to_dict()

            if _dict["pathway_id"].startswith("R-"):
                pathway_id = self.add_prefix_to_id(
                    prefix="reactome", identifier=_dict["pathway_id"]
                )
            else:
                pathway_id = self.add_prefix_to_id(
                    prefix="kegg.pathway", identifier=_dict["pathway_id"]
                )

            disease_id = self.add_prefix_to_id(
                prefix="MONDO", identifier=_dict["disease_id"]
            )

            del _dict["disease_id"], _dict["pathway_id"]

            props = {}
            for k, v in _dict.items():
                if k in self.disease_pathway_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v

            edge_list.append((None, disease_id, pathway_id, label, props))

            if self.early_stopping and index >= self.early_stopping:
                break

        return edge_list

    def get_pathway_pathway_edges(self) -> list[tuple]:
        if not hasattr(self, "compath_pathway_pathway"):
            self.download_compath_data()

        logger.info("Started writing pathway-pathway edges")

        edge_list = []
        for index, pp in tqdm(enumerate(self.compath_pathway_pathway)):
            if pp.source_db in ["kegg", "reactome"] and pp.target_db in [
                "kegg",
                "reactome",
            ]:
                if pp.relation == "isPartOf":
                    label = "pathway_is_part_of_pathway"
                elif pp.relation == "equivalentTo":
                    label = "pathway_is_equivalent_to_pathway"

                if pp.pathway_id_1.startswith("R-"):
                    pathway_id1 = self.add_prefix_to_id(
                        prefix="reactome", identifier=pp.pathway_id_1
                    )
                else:
                    pathway_id1 = self.add_prefix_to_id(
                        prefix="kegg.pathway", identifier=pp.pathway_id_1
                    )

                if pp.pathway_id_2.startswith("R-"):
                    pathway_id2 = self.add_prefix_to_id(
                        prefix="reactome", identifier=pp.pathway_id_2
                    )
                else:
                    pathway_id2 = self.add_prefix_to_id(
                        prefix="kegg.pathway", identifier=pp.pathway_id_2
                    )

                edge_list.append((None, pathway_id1, pathway_id2, label, {}))

                if self.early_stopping and index >= self.early_stopping:
                    break

        # write pathway-pathway edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Pathway_to_pathway.csv"
                )
            else:
                full_path = os.path.join(os.getcwd(), "Pathway_to_pathway.csv")

            df_list = [
                {
                    "pathway_id1": pathway_id1,
                    "pathway_id2": pathway_id2,
                    "label": label,
                }
                for _, pathway_id1, pathway_id2, label, _ in edge_list
            ]

            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Pathway-pathway edge data is written: {full_path}")

        return edge_list

    @validate_call
    def get_reactome_hierarchical_edges(
        self, label: str = "pathway_participates_pathway"
    ) -> list[tuple]:

        if not hasattr(self, "reactome_hierarchial_relations"):
            self.download_reactome_data()

        logger.info("Started writing reactome hierarchial edges")

        edge_list = []
        for index, pp in tqdm(enumerate(self.reactome_hierarchial_relations)):
            parent_id = self.add_prefix_to_id(
                prefix="reactome", identifier=pp.parent
            )
            child_id = self.add_prefix_to_id(
                prefix="reactome", identifier=pp.child
            )

            edge_list.append((None, child_id, parent_id, label, {}))

            if self.early_stopping and index >= self.early_stopping:
                break

        # write Reactome hierarchical edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Reactome_hierarchical_edges.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Reactome_hierarchical_edges.csv"
                )

            df_list = [
                {"child_id": child_id, "parent_id": parent_id, "label": label}
                for _, child_id, parent_id, label, _ in edge_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(
                f"Reactome hierarchical edge data is written: {full_path}"
            )

        return edge_list

    @validate_call
    def get_pathway_pathway_orthology_edges(
        self, label: str = "pathway_is_ortholog_to_pathway"
    ) -> list[tuple]:

        if not hasattr(self, "kegg_pathways"):
            self.download_kegg_data()

        if not hasattr(self, "reactome_pathways"):
            self.download_reactome_data()

        logger.info("Started writing pathway orthology edges")

        edge_list = []
        index = 0
        for p1 in tqdm([p for p in self.kegg_pathways if p[0][:3] == "hsa"]):
            p1_prefix_removed = p1[0][3:]

            for p2 in self.kegg_pathways:
                if p1 == p2:
                    continue

                p2_prefix_removed = p2[0][3:]

                if p1_prefix_removed == p2_prefix_removed:
                    pathway1_id = self.add_prefix_to_id(
                        prefix="kegg.pathway", identifier=p1[0]
                    )
                    pathway2_id = self.add_prefix_to_id(
                        prefix="kegg.pathway", identifier=p2[0]
                    )
                    edge_list.append(
                        (None, pathway1_id, pathway2_id, label, {})
                    )

                    index += 1

            if self.early_stopping and index >= self.early_stopping:
                break

        index = 0
        for p1 in tqdm(
            [p for p in self.reactome_pathways if p.organism == "Homo sapiens"]
        ):
            p1_id_last_element = p1.pathway_id.split("-")[-1]

            for p2 in self.reactome_pathways:
                p2_id_last_element = p2.pathway_id.split("-")[-1]

                if p1.pathway_id == p2.pathway_id:
                    continue

                if p1_id_last_element == p2_id_last_element:
                    pathway1_id = self.add_prefix_to_id(
                        prefix="reactome", identifier=p1.pathway_id
                    )
                    pathway2_id = self.add_prefix_to_id(
                        prefix="reactome", identifier=p2.pathway_id
                    )
                    edge_list.append(
                        (None, pathway1_id, pathway2_id, label, {})
                    )

                    index += 1

            if self.early_stopping and index >= self.early_stopping:
                break

        # write pathway orthology edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Pathway_orthology.csv"
                )
            else:
                full_path = os.path.join(os.getcwd(), "Pathway_orthology.csv")

            df_list = [
                {
                    "pathway1_id": pathway1_id,
                    "pathway2_id": pathway2_id,
                    "label": label,
                }
                for _, pathway1_id, pathway2_id, label, _ in edge_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Pathway orthology edge data is written: {full_path}")

        return edge_list

    def prepare_mondo_mappings(self) -> None:

        logger.debug(
            "Started preparing MONDO mappings to other disease databases"
        )

        mondo = ontology.ontology(
            ontology="mondo", fields=["is_obsolete", "obo_xref"]
        )

        self.mondo_mappings = collections.defaultdict(dict)
        mapping_db_list = ["MESH", "OMIM", "ICD10CM"]

        for term in mondo:
            if (
                not term.is_obsolete
                and term.obo_id
                and "MONDO" in term.obo_id
                and term.obo_xref
            ):
                for xref in term.obo_xref:
                    if xref.get("database") in mapping_db_list:
                        db = mapping_db_list[
                            mapping_db_list.index(xref.get("database"))
                        ]
                        self.mondo_mappings[db][xref["id"]] = term.obo_id

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

    def ensure_iterable(self, element):
        return element if isinstance(element, (list, tuple, set)) else [element]

    def set_node_and_edge_types(self, node_types, edge_types):
        self.edge_types = edge_types if edge_types is not None else list(PathwayEdgeType)
        self.node_types = node_types if node_types is not None else list(PathwayNodeType)

    def set_node_fields(self, pathway_node_fields):
        if pathway_node_fields:
            self.pathway_node_fields = [
                field.value for field in pathway_node_fields
            ]
        else:
            self.pathway_node_fields = [
                field.value for field in PathwayNodeField
            ]

    def set_edge_fields(self, 
                        protein_pathway_edge_fields, 
                        disease_pathway_edge_fields,
                        drug_pathway_edge_fields):
        if protein_pathway_edge_fields:
            self.protein_pathway_edge_fields = [
                field.value for field in protein_pathway_edge_fields
            ]
        else:
            self.protein_pathway_edge_fields = [
                field.value for field in ProteinPathwayEdgeField
            ]
        if disease_pathway_edge_fields:
            self.disease_pathway_edge_fields = [
                field.value for field in disease_pathway_edge_fields
            ]
        else:
            self.disease_pathway_edge_fields = [
                field.value for field in DiseasePathwayEdgeField
            ]
        if drug_pathway_edge_fields:
            self.drug_pathway_edge_fields = [
                field.value for field in drug_pathway_edge_fields
            ]
        else:
            self.drug_pathway_edge_fields = [
                field.value for field in DrugPathwayEdgeField
            ]
