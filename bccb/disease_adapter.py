from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import (
    pathophenodb,
    ctdbase,
    clinvar,
    # disgenet,
    chembl,
    diseases,
    opentargets,
    drugbank,
    uniprot,
    unichem,
    humsavar,
)
from . import kegg_local
from . import disgenet_local as disgenet
import json
import os
import h5py

from pypath.inputs import ontology
from pypath.formats import obo
from pypath.utils import mapping
from pypath.share import cache

from typing import Union
from pydantic import BaseModel, DirectoryPath, FilePath, EmailStr, validate_call
from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
import collections
from biocypher._logger import logger

from enum import Enum, EnumMeta, auto

import numpy as np
import pandas as pd

logger.debug(f"Loading module {__name__}.")


class DiseaseEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class DiseaseNodeField(Enum, metaclass=DiseaseEnumMeta):
    NAME = "name"
    SYNONYMS = "synonyms"

    # xrefs
    UMLS = "UMLS"
    DOID = "DOID"
    MESH = "MESH"
    OMIM = "OMIM"
    EFO = "EFO"
    ORPHANET = "Orphanet"
    HP = "HP"  # Human Phenotype Ontology
    ICD10CM = "ICD10CM"
    NCIT = "NCIT"
    ICD9 = "ICD9"
    MEDDRA = "MedDRA"

    # embedding
    DOC2VEC_EMBEDDING = "doc2vec_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def get_xrefs(cls):
        return [
            cls.UMLS.value,
            cls.DOID.value,
            cls.MESH.value,
            cls.OMIM.value,
            cls.EFO.value,
            cls.ORPHANET.value,
            cls.HP.value,
            cls.ICD10CM.value,
            cls.NCIT.value,
            cls.ICD9.value,
            cls.MEDDRA.value,
        ]


class DiseaseEdgeType(Enum, metaclass=DiseaseEnumMeta):
    MONDO_HIERARCHICAL_RELATIONS = auto()
    ORGANISM_TO_DISEASE = auto()
    GENE_TO_DISEASE = auto()
    DISEASE_TO_DRUG = auto()
    DISEASE_TO_DISEASE = auto()
    DISEASE_COMOBORDITIY = auto()


class GENE_TO_DISEASE_INTERACTION_FIELD(Enum, metaclass=DiseaseEnumMeta):
    SOURCE = "source"
    VARIANT_SOURCE = "variant_source"
    OPENTARGETS_SCORE = "opentargets_score"
    DISEASES_CONFIDENCE_SCORE = "diseases_confidence_score"
    ALLELE_ID = "allele_id"
    CLINICAL_SIGNIFICANCE = "clinical_significance"
    REVIEW_STATUS = "review_status"
    DBSNP_ID = "dbsnp_id"
    VARIATION_ID = "variation_id"
    PUBMED_IDS = "pubmed_ids"
    DISGENET_GENE_DISEASE_SCORE = "disgenet_gene_disease_score"
    DISGENET_VARIANT_DISEASE_SCORE = "disgenet_variant_disease_score"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class DISEASE_TO_DRUG_INTERACTION_FIELD(Enum, metaclass=DiseaseEnumMeta):
    SOURCE = "source"
    MAX_PHASE = "max_phase"
    PUBMED_IDS = "pubmed_ids"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class DISEASE_TO_DISEASE_INTERACTION_FIELD(Enum, metaclass=DiseaseEnumMeta):
    SOURCE = "source"
    DISGENET_JACCARD_GENES_SCORE = "disgenet_jaccard_genes_score"
    DISGENET_JACCARD_VARIANTS_SCORE = "disgenet_jaccard_variants_score"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class DiseaseModel(BaseModel):
    drugbank_user: EmailStr
    drugbank_passwd: str
    disease_node_fields: Union[list[DiseaseNodeField], None] = None
    edge_types: Union[list[DiseaseEdgeType], None] = None
    gene_disease_edge_fields: Union[
        list[GENE_TO_DISEASE_INTERACTION_FIELD], None
    ] = None
    disease_drug_edge_fields: Union[
        list[DISEASE_TO_DRUG_INTERACTION_FIELD], None
    ] = None
    disease_disease_edge_fields: Union[
        list[DISEASE_TO_DISEASE_INTERACTION_FIELD], None
    ] = None
    add_prefix: bool = True
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None


class Disease:
    """
    Class that downloads disease data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(
        self,
        drugbank_user: str,
        drugbank_passwd: str,
        disease_node_fields: Union[list[DiseaseNodeField], None] = None,
        edge_types: Union[list[DiseaseEdgeType], None] = None,
        gene_disease_edge_fields: Union[
            list[GENE_TO_DISEASE_INTERACTION_FIELD], None
        ] = None,
        disease_drug_edge_fields: Union[
            list[DISEASE_TO_DRUG_INTERACTION_FIELD], None
        ] = None,
        disease_disease_edge_fields: Union[
            list[DISEASE_TO_DISEASE_INTERACTION_FIELD], None
        ] = None,
        add_prefix: bool = True,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
    ):
        """
        Args:
            drugbank_user: drugbank username
            drugbank_passwd: drugbank password
            disease_node_fields: disease node fields that will be included in graph, if defined it must be values of elements from DiseaseNodeField enum class (not the names)
            gene_disease_edge_fields: Gene-Disease edge fields that will be included in graph, if defined it must be values of elements from GENE_TO_DISEASE_INTERACTION_FIELD enum class (not the names)
            disease_drug_edge_fields: Disease-Drug edge fields that will be included in graph, if defined it must be values of elements from DISEASE_TO_DRUG_INTERACTION_FIELD enum class (not the names)
            disease_disease_edge_fields: Disease-Disease edge fields that will be included in graph, if defined it must be values of elements from DISEASE_TO_DISEASE_INTERACTION_FIELD enum class (not the names)
            edge_types: list of edge types that will be included in graph, if defined it must be elements (not values of elements) from DiseaseEdgeType enum class
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of output data
            export_csv: if True, export data as csv
            output_dir: Location of csv export if `export_csv` is True, if not defined it will be current directory
        """
        model = DiseaseModel(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            disease_node_fields=disease_node_fields,
            edge_types=edge_types,
            gene_disease_edge_fields=gene_disease_edge_fields,
            disease_drug_edge_fields=disease_drug_edge_fields,
            disease_disease_edge_fields=disease_disease_edge_fields,
            add_prefix=add_prefix,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
        ).model_dump()

        self.drugbank_user = model["drugbank_user"]
        self.drugbank_passwd = model["drugbank_passwd"]
        self.add_prefix = model["add_prefix"]
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set node and edge field
        self.set_node_and_edge_fields(
            disease_node_fields=model["disease_node_fields"],
            gene_disease_edge_fields=model["gene_disease_edge_fields"],
            disease_drug_edge_fields=model["disease_drug_edge_fields"],
            disease_disease_edge_fields=model["disease_disease_edge_fields"],
        )

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_disease_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
        doc2vec_embedding_path: FilePath | None = None,
        malacards_json_path: FilePath | None = None,
        malacards_related_diseases_json_path: FilePath | None = None,
    ):
        """
        Wrapper function to download disease data from various databases using pypath.
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

            t0 = time()

            self.download_mondo_data()
            self.prepare_mappings()
            self.download_pathophenodb_data()
            self.download_ctd_data()
            self.download_chembl_data()
            self.download_diseases_data()
            self.download_clinvar_data()
            self.download_disgenet_data()
            self.download_opentargets_data()
            self.download_malacards_data(malacards_json_path=malacards_json_path, malacards_related_diseases_json_path=malacards_related_diseases_json_path)
            self.download_kegg_data()
            self.download_humsavar_data()

            if DiseaseNodeField.DOC2VEC_EMBEDDING.value in self.disease_node_fields:
                self.retrieve_doc2vec_embedding(doc2vec_embedding_path=doc2vec_embedding_path)

            t1 = time()
            logger.info(
                f"All data is downloaded in {round((t1-t0) / 60, 2)} mins".upper()
            )

    def download_mondo_data(self) -> None:
        fields = ["is_obsolete"]
        if DiseaseNodeField.SYNONYMS.value in self.disease_node_fields:
            fields.append("obo_synonym")

        if set(DiseaseNodeField.get_xrefs()).intersection(
            set(self.disease_node_fields)
        ):
            fields.append("obo_xref")

        t0 = time()

        self.mondo = ontology.ontology(ontology="mondo", fields=fields)

        t1 = time()
        logger.info(
            f"Mondo data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

        if DiseaseEdgeType.MONDO_HIERARCHICAL_RELATIONS in self.edge_types:
            t0 = time()

            mondo_obo_reader = obo.Obo(
                "http://purl.obolibrary.org/obo/mondo.obo"
            )

            mondo_obo_reader.parent_terms()

            self.mondo_hierarchical_relations = {
                k: v for k, v in mondo_obo_reader.parents.items() if v
            }

            t1 = time()
            logger.info(
                f"Mondo hierarchical relations data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_pathophenodb_data(self) -> None:
        if DiseaseEdgeType.ORGANISM_TO_DISEASE in self.edge_types:
            t0 = time()

            self.pathopheno_organism_disease_int = (
                pathophenodb.disease_pathogen_interactions()
            )

            t1 = time()
            logger.info(
                f"PathophenoDB organism-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_ctd_data(self) -> None:
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()

            self.ctdbase_gda = ctdbase.ctdbase_relations(
                relation_type="gene_disease"
            )

            t1 = time()
            logger.info(
                f"CTD gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

        if DiseaseEdgeType.DISEASE_TO_DRUG in self.edge_types:
            t0 = time()

            self.ctdbase_cd = ctdbase.ctdbase_relations(
                relation_type="chemical_disease"
            )

            self.drugbank_data = drugbank.DrugbankFull(
                user=self.drugbank_user, passwd=self.drugbank_passwd
            )
            drugbank_drugs_detailed = self.drugbank_data.drugbank_drugs_full(
                fields=["cas_number"]
            )
            self.cas_to_drugbank = {
                drug.cas_number: drug.drugbank_id
                for drug in drugbank_drugs_detailed
                if drug.cas_number
            }

            t1 = time()
            logger.info(
                f"CTD chemical-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_chembl_data(self) -> None:
        if DiseaseEdgeType.DISEASE_TO_DRUG in self.edge_types:
            t0 = time()

            self.chembl_disease_drug = chembl.chembl_drug_indications()

            self.chembl_to_drugbank = {
                k: list(v)[0]
                for k, v in unichem.unichem_mapping(
                    "chembl", "drugbank"
                ).items()
            }

            t1 = time()
            logger.info(
                f"CHEMBL drug indication data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_kegg_data(self) -> None:
        if DiseaseEdgeType.DISEASE_TO_DRUG in self.edge_types:
            t0 = time()

            self.kegg_drug_disease = kegg_local.drug_to_disease()

            if not hasattr(self, "drugbank_data"):
                self.drugbank_data = drugbank.DrugbankFull(
                    user=self.drugbank_user, passwd=self.drugbank_passwd
                )

            drugbank_drugs_external_ids = (
                self.drugbank_data.drugbank_external_ids_full()
            )
            self.kegg_drug_to_drugbank = {
                v.get("KEGG Drug"): k
                for k, v in drugbank_drugs_external_ids.items()
                if v.get("KEGG Drug")
            }

            kegg_disease_ids = kegg_local._Disease()._data.keys()

            self.kegg_diseases_mappings = {}
            for dis in kegg_disease_ids:
                try:
                    result = kegg_local.get_diseases(dis)
                    self.kegg_diseases_mappings[dis] = result[0].db_links
                except (IndexError, UnicodeDecodeError) as e:
                    logger.debug(f"{dis} is not available due to {e}")

            t1 = time()
            logger.info(
                f"KEGG drug indication data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()

            self.kegg_gene_disease = kegg_local.gene_to_disease(org="hsa")

            self.kegg_gene_id_to_entrez = (
                kegg_local.kegg_gene_id_to_ncbi_gene_id("hsa")
            )

            if not hasattr(self, "kegg_diseases_mappings"):
                kegg_disease_ids = kegg_local._Disease()._data.keys()

                self.kegg_diseases_mappings = {}
                for dis in kegg_disease_ids:
                    result = kegg_local.get_diseases(dis)
                    self.kegg_diseases_mappings[dis] = result[0].db_links

            t1 = time()
            logger.info(
                f"KEGG gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_opentargets_data(self) -> None:
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()

            self.opentargets_direct = opentargets.opentargets_direct_score()

            uniprot_to_entrez = uniprot.uniprot_data(
                "xref_geneid", "9606", True
            )
            self.uniprot_to_entrez = {
                k: v.strip(";").split(";")[0]
                for k, v in uniprot_to_entrez.items()
            }

            if not hasattr(self, "ensembl_gene_to_uniprot"):
                uniprot_to_ensembl = uniprot.uniprot_data(
                    "xref_ensembl", "9606", True
                )

                self.ensembl_gene_to_uniprot = {
                    self.ensembl_transcript_to_ensembl_gene(ensts): uniprot_id
                    for uniprot_id, ensts in uniprot_to_ensembl.items()
                    if self.ensembl_transcript_to_ensembl_gene(ensts)
                }

            t1 = time()
            logger.info(
                f"Open Targets direct gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_diseases_data(self) -> None:
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()

            self.diseases_knowledge = diseases.knowledge_filtered()
            self.diseases_experimental = diseases.experiments_filtered()

            if not hasattr(self, "ensembl_gene_to_uniprot"):
                uniprot_to_ensembl = uniprot.uniprot_data(
                    "xref_ensembl", "9606", True
                )

                self.ensembl_gene_to_uniprot = {
                    self.ensembl_transcript_to_ensembl_gene(ensts): uniprot_id
                    for uniprot_id, ensts in uniprot_to_ensembl.items()
                    if self.ensembl_transcript_to_ensembl_gene(ensts)
                }

            self.ensembl_protein_to_uniprot = {}
            for k, v in self.ensembl_gene_to_uniprot.items():
                ensps = self.ensembl_gene_to_ensembl_protein(k)
                for p in ensps:
                    self.ensembl_protein_to_uniprot[p] = v

            t1 = time()
            logger.info(
                f"DISEASES gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_clinvar_data(self) -> None:
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()

            # with curl.cache_off():
            self.clinvar_variant_disease = clinvar.clinvar_raw()

            self.clinvar_citation = clinvar.clinvar_citations()

            t1 = time()
            logger.info(
                f"Clinvar variant-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_humsavar_data(self) -> None:
        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()

            self.humsavar_data = humsavar.uniprot_variants()

            if not hasattr(self, "uniprot_to_entrez"):
                uniprot_to_entrez = uniprot.uniprot_data(
                    "xref_geneid", "9606", True
                )
                self.uniprot_to_entrez = {
                    k: v.strip(";").split(";")[0]
                    for k, v in uniprot_to_entrez.items()
                }

            t1 = time()
            logger.info(
                f"Humsavar variant-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_disgenet_data(self) -> None:
        if not hasattr(self, "disgenet_id_mappings_dict"):
            self.prepare_disgenet_id_mappings()

        if not hasattr(self, "uniprot_to_entrez"):
            uniprot_to_entrez = uniprot.uniprot_data(
                "xref_geneid", "9606", True
            )
            self.uniprot_to_entrez = {
                k: v.strip(";").split(";")[0]
                for k, v in uniprot_to_entrez.items()
            }

        self.gene_symbol_to_uniprot = {}
        for k, v in uniprot.uniprot_data("gene_names", "9606", True).items():
            for symbol in v.split(" "):
                self.gene_symbol_to_uniprot[symbol] = k

        if DiseaseEdgeType.DISEASE_TO_DISEASE in self.edge_types:
            t0 = time()

            if not hasattr(self, "disgenet_api"):
                self.disgenet_api = disgenet.DisgenetApi()

            if not hasattr(self, "disgenet_disease_ids"):
                self.disgenet_disease_ids = (
                    disgenet.disease_id_mappings().keys()
                )

            if not hasattr(self, "disgenet_id_mappings_dict"):
                self.prepare_disgenet_id_mappings()

            self.disgenet_dda_gene = []
            self.disgenet_dda_variant = []
            for disease_id in tqdm(self.disgenet_disease_ids):
                try:
                    self.disgenet_dda_gene.extend(
                        self.disgenet_api.get_ddas_that_share_genes(disease_id)
                    )
                    self.disgenet_dda_variant.extend(
                        self.disgenet_api.get_ddas_that_share_variants(
                            disease_id
                        )
                    )
                except TypeError:
                    logger.debug(f"{disease_id} not available")

            t1 = time()
            logger.info(
                f"Disgenet disease-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            t0 = time()

            if not hasattr(self, "disgenet_api"):
                self.disgenet_api = disgenet.DisgenetApi()

            if not hasattr(self, "disgenet_disease_ids") or not hasattr(
                self, "disgenet_id_mappings_dict"
            ):
                self.disgenet_disease_ids = (
                    disgenet.disease_id_mappings().keys()
                )

            if not hasattr(self, "disgenet_id_mappings_dict"):
                self.prepare_disgenet_id_mappings()

            if not hasattr(self, "uniprot_to_entrez"):
                uniprot_to_entrez = uniprot.uniprot_data(
                    "xref_geneid", "9606", True
                )
                self.uniprot_to_entrez = {
                    k: v.strip(";").split(";")[0]
                    for k, v in uniprot_to_entrez.items()
                }

            self.gene_symbol_to_uniprot = {}
            for k, v in uniprot.uniprot_data(
                "gene_names", "9606", True
            ).items():
                for symbol in v.split(" "):
                    self.gene_symbol_to_uniprot[symbol] = k

            self.disgenet_gda = []
            self.disgenet_vda = []
            for disease_id in tqdm(self.disgenet_disease_ids):
                try:
                    self.disgenet_gda.extend(
                        self.disgenet_api.get_gdas_by_diseases(disease_id)
                    )
                    self.disgenet_vda.extend(
                        self.disgenet_api.get_vdas_by_diseases(disease_id)
                    )

                except (TypeError, ValueError) as e:
                    logger.debug(f"{disease_id} not available")

            t1 = time()
            logger.info(
                f"Disgenet gene-disease interaction data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_malacards_data(self, 
                                malacards_json_path: FilePath | None = None, 
                                malacards_related_diseases_json_path: FilePath | None = None) -> None:

        if DiseaseEdgeType.DISEASE_COMOBORDITIY in self.edge_types:
            t0 = time()
            if not malacards_json_path:
                malacards_json_path = os.path.join(
                    cache.get_cachedir(), "MalaCards.json"
                )
            with open(malacards_json_path, encoding="utf-8") as file:
                file_content = file.read()

            malacards_external_ids = json.loads(file_content)

            self.prepare_malacards_mondo_mappings(malacards_external_ids)

            self.malacards_disease_slug_to_malacards_id = {
                entry["DiseaseSlug"]: entry["McId"]
                for entry in malacards_external_ids
            }
            if not malacards_related_diseases_json_path:
                malacards_related_diseases_json_path = os.path.join(
                    cache.get_cachedir(), "MalaCardsRelatedDiseases.json"
                )

            with open(
                malacards_related_diseases_json_path, encoding="utf-8"
            ) as file:
                file_content = file.read()

            self.disease_comorbidity = json.loads(file_content)

            t1 = time()
            logger.info(
                f"Malacards disease comorbidity data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def retrieve_doc2vec_embedding(self, doc2vec_embedding_path: FilePath) -> None:
        logger.info("Retrieving Doc2Vec disease embeddings.")

        self.disease_id_to_doc2vec_embedding = {}
        with h5py.File(doc2vec_embedding_path, "r") as f:
            for mondo_id, embedding in f.items():
                self.hpo_id_to_cada_embedding[mondo_id] = np.array(embedding).astype(np.float16)

    def prepare_mappings(self) -> None:
        """
        Prepare disease database identifier mappings from MONDO to other disease databases
        """
        self.mondo_mappings = collections.defaultdict(dict)
        mapping_db_list = [
            "UMLS",
            "DOID",
            "MESH",
            "OMIM",
            "EFO",
            "Orphanet",
            "HP",
            "ICD10CM",
            "NCIT",
        ]

        if not hasattr(self, "mondo"):
            self.download_mondo_data()

        for term in self.mondo:

            if (
                not term.is_obsolete
                and term.obo_id
                and "MONDO" in term.obo_id
                and hasattr(term, "obo_xref")
                and term.obo_xref
            ):
                for xref in term.obo_xref:
                    if xref.get("database") in mapping_db_list:
                        db = mapping_db_list[
                            mapping_db_list.index(xref.get("database"))
                        ]
                        self.mondo_mappings[db][xref["id"]] = term.obo_id

    def prepare_disgenet_id_mappings(self):
        """
        Prepare disgenet id mappings
        """

        disgenet_id_mappings = disgenet.disease_id_mappings()

        selected_dbs = [
            "DO",
            "EFO",
            "HPO",
            "MONDO",
            "MSH",
            "NCI",
            "ICD10CM",
            "OMIM",
        ]
        self.disgenet_id_mappings_dict = collections.defaultdict(dict)

        for disg_id, mappings in disgenet_id_mappings.items():
            map_dict = {}
            for m in mappings.vocabularies:
                if (
                    m.vocabulary in selected_dbs
                    and m.vocabulary not in map_dict.keys()
                ):
                    map_dict[m.vocabulary] = (
                        m.code.split(":")[1] if ":" in m.code else m.code
                    )

            self.disgenet_id_mappings_dict[disg_id] = map_dict

    def prepare_malacards_mondo_mappings(self, malacards_external_ids):

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()

        malacards_dbs_to_mondo_dbs = {
            "OMIM®": "OMIM",
            "Disease Ontology": "DOID",
            "UMLS": "UMLS",
            "MeSH": "MESH",
            "NCIt": "NCIT",
            "EFO": "EFO",
            "Orphanet": "Orphanet",
            "ICD10": "ICD10CM",
        }

        self.malacards_id_to_mondo_id = {}
        for entry in malacards_external_ids:
            if entry["ExternalIds"]:
                for external_id in entry["ExternalIds"]:

                    if (
                        external_id["Source"]
                        in malacards_dbs_to_mondo_dbs.keys()
                    ):
                        malacards_external_id = external_id["SourceAccession"]

                        if external_id["Source"] == "Disease Ontology":
                            malacards_external_id = malacards_external_id.split(
                                ":"
                            )[1]
                        if external_id["Source"] == "EFO":
                            malacards_external_id = malacards_external_id.split(
                                "_"
                            )[1]
                        if external_id["Source"] == "Orphanet":
                            malacards_external_id = (
                                malacards_external_id.replace("ORPHA", "")
                            )

                        db = malacards_dbs_to_mondo_dbs[external_id["Source"]]

                        if self.mondo_mappings[db].get(malacards_external_id):
                            self.malacards_id_to_mondo_id[entry["McId"]] = (
                                self.mondo_mappings[db].get(
                                    malacards_external_id
                                )
                            )
                            break

    def process_ctd_chemical_disease(self) -> pd.DataFrame:

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()

        if not hasattr(self, "ctdbase_cd"):
            self.download_ctd_data()

        logger.debug("Started processing CTD chemical-disease data")
        t0 = time()

        df_list = []
        for interaction in tqdm(self.ctdbase_cd):
            if (
                interaction.CasRN
                and interaction.DirectEvidence
                and interaction.DirectEvidence == "therapeutic"
                and interaction.PubMedIDs
                and self.cas_to_drugbank.get(interaction.CasRN)
            ):
                db = interaction.DiseaseID.split(":")[0]
                disease_id = interaction.DiseaseID.split(":")[1]
                if self.mondo_mappings[db].get(disease_id, None):
                    disease_id = self.mondo_mappings[db].get(disease_id)
                    drug_id = self.cas_to_drugbank.get(interaction.CasRN)

                    if isinstance(interaction.PubMedIDs, list):
                        pubmed_ids = "|".join(interaction.PubMedIDs)
                    else:
                        pubmed_ids = interaction.PubMedIDs

                    df_list.append((disease_id, drug_id, pubmed_ids))

        df = pd.DataFrame(
            df_list, columns=["disease_id", "drug_id", "pubmed_ids"]
        )
        df["source"] = "CTD"

        df = df.groupby(
            ["disease_id", "drug_id"], sort=False, as_index=False
        ).aggregate(
            {
                "disease_id": "first",
                "drug_id": "first",
                "pubmed_ids": self.merge_source_column,
                "source": "first",
            }
        )
        t1 = time()
        logger.info(
            f"CTD chemical-disease interaction data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_chembl_drug_indication(self) -> pd.DataFrame:

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "chembl_disease_drug"):
            self.download_chembl_data()

        logger.debug("Started processing CHEMBL drug indication data")
        t0 = time()

        df_list = []
        for dd in tqdm(self.chembl_disease_drug):
            if (
                dd.efo_id
                and self.chembl_to_drugbank.get(dd.molecule_chembl)
                and dd.efo_id.split(":")[0]
                in list(self.mondo_mappings.keys()) + ["MONDO"]
                and dd.max_phase > 0.0
            ):
                db = dd.efo_id.split(":")[0]
                disease_id = dd.efo_id.split(":")[1]
                drug_id = self.chembl_to_drugbank.get(dd.molecule_chembl)

                if db == "MONDO":
                    df_list.append((dd.efo_id, drug_id, dd.max_phase))
                elif self.mondo_mappings[db].get(disease_id):
                    disease_id = self.mondo_mappings[db].get(disease_id)
                    df_list.append((disease_id, drug_id, dd.max_phase))

        df = pd.DataFrame(
            df_list, columns=["disease_id", "drug_id", "max_phase"]
        )
        df["source"] = "ChEMBL"

        df.sort_values(
            by="max_phase", ascending=False, ignore_index=True, inplace=True
        )

        df.drop_duplicates(
            subset=["disease_id", "drug_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"CHEMBL drug indication data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_kegg_drug_indication(self) -> pd.DataFrame:
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "kegg_drug_disease"):
            self.download_kegg_data()

        logger.debug("Started processing CHEMBL drug indication data")
        t0 = time()

        kegg_dbs_to_mondo_dbs = {
            "MeSH": "MESH",
            "OMIM": "OMIM",
            "ICD-10": "ICD10CM",
        }

        df_list = []
        for drug, kegg_diseases in self.kegg_drug_disease.items():
            if self.kegg_drug_to_drugbank.get(drug):
                for interaction in kegg_diseases.DiseaseEntries:
                    disease_id = None
                    if self.kegg_diseases_mappings.get(interaction.disease_id):
                        found = False

                        for db in kegg_dbs_to_mondo_dbs.keys():
                            if found:
                                break

                            if self.kegg_diseases_mappings[
                                interaction.disease_id
                            ].get(db):
                                for ref in self.ensure_iterable(
                                    self.kegg_diseases_mappings[
                                        interaction.disease_id
                                    ][db]
                                ):
                                    if found:
                                        break

                                    if self.mondo_mappings[
                                        kegg_dbs_to_mondo_dbs[db]
                                    ].get(ref):
                                        disease_id = self.mondo_mappings[
                                            kegg_dbs_to_mondo_dbs[db]
                                        ].get(ref)
                                        found = True

                    if disease_id:
                        drug_id = self.kegg_drug_to_drugbank.get(drug)
                        df_list.append((disease_id, drug_id))

        df = pd.DataFrame(
            df_list,
            columns=[
                "disease_id",
                "drug_id",
            ],
        )
        df["source"] = "KEGG"

        # Doesnt look necessary
        df.drop_duplicates(
            subset=["disease_id", "drug_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"KEGG drug indication data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_opentargets_gene_disease(self) -> pd.DataFrame:
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "opentargets_direct"):
            self.download_opentargets_data()

        logger.debug(
            "Started processing OpenTargets gene-disease interaction data"
        )
        t0 = time()

        df_list = []
        for entry in self.opentargets_direct:
            disease = entry["diseaseId"]
            db = entry["diseaseId"].split("_")[0]
            target = entry["targetId"]
            if self.uniprot_to_entrez.get(
                self.ensembl_gene_to_uniprot.get(target)
            ) and self.mondo_mappings[db].get(disease.split("_")[1]):
                score = round(entry["score"], 3)
                if score != 0.0:
                    gene_id = self.uniprot_to_entrez.get(
                        self.ensembl_gene_to_uniprot.get(target)
                    )
                    disease_id = self.mondo_mappings[db].get(
                        disease.split("_")[1]
                    )
                    df_list.append((gene_id, disease_id, score))

        df = pd.DataFrame(
            df_list, columns=["gene_id", "disease_id", "opentargets_score"]
        )
        df["source"] = "Open Targets"

        df.sort_values(
            by="opentargets_score",
            ascending=False,
            ignore_index=True,
            inplace=True,
        )

        df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"OpenTargets gene-disease interaction data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_diseases_gene_disease(self) -> pd.DataFrame:

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "diseases_knowledge") or not hasattr(
            self, "diseases_experimental"
        ):
            self.download_diseases_data()

        logger.debug("Started processing DISEASES gene-disease data")
        t0 = time()

        df_list = []
        for dg in tqdm(self.diseases_knowledge):
            if self.ensembl_protein_to_uniprot.get(
                dg.gene_id
            ) and self.mondo_mappings["DOID"].get(dg.disease_id.split(":")[1]):

                gene_id = self.ensembl_protein_to_uniprot.get(dg.gene_id)
                disease_id = self.mondo_mappings["DOID"].get(
                    dg.disease_id.split(":")[1]
                )
                df_list.append((gene_id, disease_id))

        diseases_knowledge_df = pd.DataFrame(
            df_list,
            columns=[
                "gene_id",
                "disease_id",
            ],
        )
        diseases_knowledge_df["source"] = "DISEASES Knowledge"

        # drop duplicates
        diseases_knowledge_df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        df_list = []
        for dg in tqdm(self.diseases_experimental):
            if self.ensembl_protein_to_uniprot.get(
                dg.gene_id
            ) and self.mondo_mappings["DOID"].get(dg.disease_id.split(":")[1]):
                gene_id = self.ensembl_protein_to_uniprot.get(dg.gene_id)
                disease_id = self.mondo_mappings["DOID"].get(
                    dg.disease_id.split(":")[1]
                )
                score = float(dg.confidence)
                df_list.append((gene_id, disease_id, score))

        diseases_experimental_df = pd.DataFrame(
            df_list,
            columns=["gene_id", "disease_id", "diseases_confidence_score"],
        )
        diseases_experimental_df["source"] = "DISEASES Experimental"

        diseases_experimental_df.sort_values(
            by="diseases_confidence_score",
            ascending=False,
            ignore_index=True,
            inplace=True,
        )

        # drop duplicates
        diseases_experimental_df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"DISEASES gene-disease data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return diseases_knowledge_df, diseases_experimental_df

    def process_clinvar_gene_disease(self) -> pd.DataFrame:
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "clinvar_variant_disease"):
            self.download_clinvar_data()

        logger.debug("Started processing Clinvar variant-disease data")
        t0 = time()

        selected_clinical_significances = [
            "Pathogenic",
            "Likely pathogenic",
            "Pathogenic/Likely pathogenic",
        ]
        selected_review_status = [
            "criteria provided, multiple submitters, no conflicts",
            "reviewed by expert panel",
            "practice guideline",
        ]
        review_status_dict = {
            "criteria provided, multiple submitters, no conflicts": 2,
            "reviewed by expert panel": 3,
            "practice guideline": 4,
        }
        clinvar_dbs_to_mondo_dbs = {
            "MONDO": "MONDO",
            "OMIM": "OMIM",
            "Orphanet": "Orphanet",
            "HP": "HP",
            "MeSH": "MESH",
        }

        df_list = []
        for var in tqdm(self.clinvar_variant_disease):
            if (
                var.entrez
                and var.clinical_significance in selected_clinical_significances
                and var.review_status in selected_review_status
            ):
                diseases_set = set()
                for phe in var.phenotype_ids:
                    dbs_and_ids = list(dict.fromkeys(phe.split(":")).keys())

                    if len(dbs_and_ids) > 2:
                        dbs_and_ids = phe.split(":")[1:]

                    if dbs_and_ids[0] in clinvar_dbs_to_mondo_dbs.keys():
                        if dbs_and_ids[0] == "MONDO":
                            diseases_set.add("MONDO:" + dbs_and_ids[1])
                        else:
                            if self.mondo_mappings[
                                clinvar_dbs_to_mondo_dbs[dbs_and_ids[0]]
                            ].get(dbs_and_ids[1]):
                                diseases_set.add(
                                    self.mondo_mappings[
                                        clinvar_dbs_to_mondo_dbs[dbs_and_ids[0]]
                                    ].get(dbs_and_ids[1])
                                )

                if diseases_set:
                    review_status = review_status_dict[var.review_status]
                    dbsnp_id = var.rs
                    if dbsnp_id:
                        dbsnp_id = "rs" + str(dbsnp_id)

                    for d in diseases_set:
                        df_list.append(
                            (
                                var.entrez,
                                d,
                                var.allele,
                                var.clinical_significance,
                                review_status,
                                dbsnp_id,
                                var.variation_id,
                            )
                        )

        df = pd.DataFrame(
            df_list,
            columns=[
                "gene_id",
                "disease_id",
                "allele_id",
                "clinical_significance",
                "review_status",
                "dbsnp_id",
                "variation_id",
            ],
        )
        df["source"] = "Clinvar"
        df["variant_source"] = "Clinvar"

        df_list = []
        for cit in self.clinvar_citation:
            if cit.citation_source in ["PubMed", "PubMedCentral"]:
                df_list.append((cit.allele, cit.variation_id, cit.citation_id))

        clinvar_citation_df = pd.DataFrame(
            df_list, columns=["allele_id", "variation_id", "pubmed_ids"]
        )
        clinvar_citation_df = clinvar_citation_df.groupby(
            ["allele_id", "variation_id"], sort=False, as_index=False
        ).aggregate(
            {
                "allele_id": "first",
                "variation_id": "first",
                "pubmed_ids": self.merge_source_column,
            }
        )

        df = df.merge(
            clinvar_citation_df, how="left", on=["allele_id", "variation_id"]
        )

        df.sort_values(
            by="review_status", ascending=False, ignore_index=True, inplace=True
        )

        df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"Clinvar variant-disease data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_humsavar_gene_disease(self) -> pd.DataFrame:
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "humsavar_data"):
            self.download_humsavar_data()

        logger.debug("Started processing Humsavar variant-disease data")
        t0 = time()

        df_list = []
        for protein, variant_set in tqdm(self.humsavar_data.items()):
            for variant in variant_set:
                if (
                    variant.variant_category == "LP/P"
                    and variant.disease_omim
                    and variant.dbsnp
                    and self.uniprot_to_entrez.get(protein)
                    and self.mondo_mappings["OMIM"].get(
                        variant.disease_omim.split(":")[1]
                    )
                ):
                    gene_id = self.uniprot_to_entrez.get(protein)
                    disease_id = self.mondo_mappings["OMIM"].get(
                        variant.disease_omim.split(":")[1]
                    )
                    df_list.append((gene_id, disease_id, variant.dbsnp))

        df = pd.DataFrame(
            df_list, columns=["gene_id", "disease_id", "dbsnp_id"]
        )
        df["source"] = "Humsavar"
        df["variant_source"] = "Humsavar"

        df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"Humsavar variant-disease data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_kegg_gene_disease(self) -> pd.DataFrame:
        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "kegg_gene_disease"):
            self.download_kegg_data()

        logger.debug("Started processing KEGG gene-disease data")
        t0 = time()

        kegg_dbs_to_mondo_dbs = {
            "MeSH": "MESH",
            "OMIM": "OMIM",
            "ICD-10": "ICD10CM",
        }

        df_list = []
        for gene, kegg_diseases in self.kegg_gene_disease.items():
            if self.kegg_gene_id_to_entrez.get(gene):
                for interaction in kegg_diseases.DiseaseEntries:
                    disease_id = None
                    if self.kegg_diseases_mappings.get(interaction.disease_id):
                        found = False

                        for db in kegg_dbs_to_mondo_dbs.keys():
                            if found:
                                break

                            if self.kegg_diseases_mappings[
                                interaction.disease_id
                            ].get(db):
                                for ref in self.ensure_iterable(
                                    self.kegg_diseases_mappings[
                                        interaction.disease_id
                                    ][db]
                                ):
                                    if found:
                                        break

                                    if self.mondo_mappings[
                                        kegg_dbs_to_mondo_dbs[db]
                                    ].get(ref):
                                        disease_id = self.mondo_mappings[
                                            kegg_dbs_to_mondo_dbs[db]
                                        ].get(ref)
                                        found = True

                    if disease_id:
                        gene_id = self.kegg_gene_id_to_entrez.get(gene)
                        df_list.append((gene_id, disease_id))

        df = pd.DataFrame(df_list, columns=["gene_id", "disease_id"])
        df["source"] = "KEGG"

        df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"KEGG gene-disease data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_disgenet_gene_disease(self) -> pd.DataFrame:

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "disgenet_id_mappings_dict"):
            self.prepare_disgenet_id_mappings()
        if not hasattr(self, "disgenet_gda") or not hasattr(
            self, "disgenet_vda"
        ):
            self.download_disgenet_data()

        logger.debug("Started processing Disgenet gene-disease data")
        t0 = time()

        df_list = []
        for gda in self.disgenet_gda:
            diseaseid = self.map_disgenet_disease_id_to_mondo_id(
                gda.diseaseid, return_pandas_none=False
            )

            if diseaseid and gda.geneid:
                df_list.append((str(gda.geneid), diseaseid, gda.score))

        disgenet_gda_df = pd.DataFrame(
            df_list,
            columns=["gene_id", "disease_id", "disgenet_gene_disease_score"],
        )
        disgenet_gda_df["source"] = "Disgenet Gene-Disease"

        # DOES NOT LOOK NECESSARY
        disgenet_gda_df.sort_values(
            by="disgenet_gene_disease_score",
            ascending=False,
            ignore_index=True,
            inplace=True,
        )
        disgenet_gda_df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        df_list = []
        for vda in self.disgenet_vda:
            if vda.gene_symbol and self.uniprot_to_entrez.get(
                self.gene_symbol_to_uniprot.get(vda.gene_symbol)
            ):
                diseaseid = self.map_disgenet_disease_id_to_mondo_id(
                    vda.diseaseid, return_pandas_none=False
                )
                gene_id = str(
                    self.uniprot_to_entrez.get(
                        self.gene_symbol_to_uniprot.get(vda.gene_symbol)
                    )
                )
                if diseaseid and gene_id:
                    df_list.append(
                        (gene_id, diseaseid, vda.score, vda.variantid)
                    )

        disgenet_vda_df = pd.DataFrame(
            df_list,
            columns=[
                "gene_id",
                "disease_id",
                "disgenet_variant_disease_score",
                "dbsnp_id",
            ],
        )
        disgenet_vda_df["source"] = "Disgenet Variant-Disease"
        disgenet_vda_df["variant_source"] = "Disgenet Variant-Disease"

        disgenet_vda_df.sort_values(
            by="disgenet_variant_disease_score",
            ascending=False,
            ignore_index=True,
            inplace=True,
        )
        disgenet_vda_df.drop_duplicates(
            subset=["gene_id", "disease_id"], ignore_index=True, inplace=True
        )

        t1 = time()
        logger.info(
            f"Disgenet gene-disease data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return disgenet_gda_df, disgenet_vda_df

    def process_disgenet_disease_disease(self) -> pd.DataFrame:

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()
        if not hasattr(self, "disgenet_id_mappings_dict"):
            self.prepare_disgenet_id_mappings()
        if not hasattr(self, "disgenet_dda_gene") or not hasattr(
            self, "disgenet_dda_variant"
        ):
            self.download_disgenet_data()

        logger.debug("Started processing Disgenet disease-disease data")
        t0 = time()

        df_list = []
        # DISEASE-DISEASE BY GENE
        for dda in self.disgenet_dda_gene:
            if round(dda.jaccard_genes, 3) != 0.0:
                diseaseid1 = self.map_disgenet_disease_id_to_mondo_id(
                    dda.diseaseid1, return_pandas_none=False
                )
                diseaseid2 = self.map_disgenet_disease_id_to_mondo_id(
                    dda.diseaseid2, return_pandas_none=False
                )

                if diseaseid1 and diseaseid2:
                    df_list.append(
                        (
                            diseaseid1,
                            diseaseid2,
                            round(dda.jaccard_genes, 3),
                        )
                    )

        disgenet_dda_gene_df = pd.DataFrame(
            df_list,
            columns=[
                "disease_id1",
                "disease_id2",
                "disgenet_jaccard_genes_score",
            ],
        )
        disgenet_dda_gene_df["source"] = "Disgenet Disease-Disease Gene"

        disgenet_dda_gene_df.sort_values(
            by="disgenet_jaccard_genes_score",
            ascending=False,
            ignore_index=True,
            inplace=True,
        )
        disgenet_dda_gene_df = disgenet_dda_gene_df[
            ~disgenet_dda_gene_df[["disease_id1", "disease_id2"]]
            .apply(frozenset, axis=1)
            .duplicated()
        ].reset_index(drop=True)

        df_list = []
        # DISEASE-DISEASE BY VARIANT
        for dda in self.disgenet_dda_variant:
            if round(dda.jaccard_variants, 3) != 0.0:
                diseaseid1 = self.map_disgenet_disease_id_to_mondo_id(
                    dda.diseaseid1, return_pandas_none=False
                )
                diseaseid2 = self.map_disgenet_disease_id_to_mondo_id(
                    dda.diseaseid2, return_pandas_none=False
                )

                if diseaseid1 and diseaseid2:
                    df_list.append(
                        (
                            diseaseid1,
                            diseaseid2,
                            round(dda.jaccard_variants, 3),
                        )
                    )

        disgenet_dda_variant_df = pd.DataFrame(
            df_list,
            columns=[
                "disease_id1",
                "disease_id2",
                "disgenet_jaccard_variants_score",
            ],
        )
        disgenet_dda_variant_df["source"] = "Disgenet Disease-Disease Variant"

        disgenet_dda_variant_df.sort_values(
            by="disgenet_jaccard_variants_score",
            ascending=False,
            ignore_index=True,
            inplace=True,
        )
        disgenet_dda_variant_df = disgenet_dda_variant_df[
            ~disgenet_dda_variant_df[["disease_id1", "disease_id2"]]
            .apply(frozenset, axis=1)
            .duplicated()
        ].reset_index(drop=True)

        t1 = time()
        logger.info(
            f"Disgenet disease-disease data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return disgenet_dda_gene_df, disgenet_dda_variant_df

    def process_malacards_disease_comorbidity(self) -> pd.DataFrame:
        if not hasattr(self, "disease_comorbidity"):
            self.download_malacards_data()

        logger.info("Started processing Malacards disease comorbidity data")
        t0 = time()

        df_list = []
        for disease in self.disease_comorbidity:
            if disease["Comorbidities"] and self.malacards_id_to_mondo_id.get(
                disease["McId"]
            ):
                for comorbidity in disease["Comorbidities"]:
                    if self.malacards_id_to_mondo_id.get(
                        self.malacards_disease_slug_to_malacards_id.get(
                            comorbidity["DiseaseSlug"]
                        )
                    ):
                        disease1 = self.malacards_id_to_mondo_id.get(
                            disease["McId"]
                        )
                        disease2 = self.malacards_id_to_mondo_id.get(
                            self.malacards_disease_slug_to_malacards_id.get(
                                comorbidity["DiseaseSlug"]
                            )
                        )
                        df_list.append((disease1, disease2))

        comorbidity_df = pd.DataFrame(df_list, columns=["disease1", "disease2"])

        comorbidity_df = comorbidity_df[
            comorbidity_df["disease1"].ne(comorbidity_df["disease2"])
        ]

        comorbidity_df = comorbidity_df[
            ~comorbidity_df[["disease1", "disease2"]]
            .apply(frozenset, axis=1)
            .duplicated()
        ].reset_index(drop=True)

        t1 = time()
        logger.info(
            f"Malacards disease comorbidity data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        # write disease-disease comorbidity edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Disease_to_disease_comorbidity_edge.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Disease_to_disease_comorbidity_edge.csv"
                )

            comorbidity_df.to_csv(full_path, index=False)
            logger.info(
                f"Disease-Disease comorbidity edge data is written: {full_path}"
            )

        return comorbidity_df

    def merge_disease_drug_edge_data(self) -> pd.DataFrame:
        # Prepare dataframes for merging
        ctd_df = self.process_ctd_chemical_disease()

        chembl_df = self.process_chembl_drug_indication()

        kegg_df = self.process_kegg_drug_indication()

        logger.debug("Started merging disease-drug data")
        t0 = time()

        # MERGE CHEMBL AND CTD
        merged_df = chembl_df.merge(
            ctd_df, how="outer", on=["disease_id", "drug_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        logger.debug("CHEMBL and CTD disease-drug data is merged")

        # MERGE CHEMBL+CTD AND KEGG
        merged_df = merged_df.merge(
            kegg_df, how="outer", on=["disease_id", "drug_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        t1 = time()
        logger.info(
            f"Disease-drug edge data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        # write disease-drug edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Disease_to_drug_edge.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Disease_to_drug_edge.csv"
                )

            merged_df.to_csv(full_path, index=False)
            logger.info(f"Disease-Drug edge data is written: {full_path}")

        return merged_df

    def merge_gene_disease_edge_data(self) -> pd.DataFrame:

        opentargets_df = self.process_opentargets_gene_disease()

        diseases_knowledge_df, diseases_experimental_df = (
            self.process_diseases_gene_disease()
        )

        humsavar_df = self.process_humsavar_gene_disease()

        kegg_df = self.process_kegg_gene_disease()

        disgenet_gda_df, disgenet_vda_df = self.process_disgenet_gene_disease()

        clinvar_df = self.process_clinvar_gene_disease()

        logger.debug("Started merging gene-disease data")
        t0 = time()

        # MERGE DISEASES KNOWLEDGE AND EXPERIMENTAL DATA
        diseases_df = pd.merge(
            diseases_knowledge_df,
            diseases_experimental_df,
            how="outer",
            on=["gene_id", "disease_id"],
        )

        diseases_df["source"] = diseases_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        diseases_df.drop(columns=["source_x", "source_y"], inplace=True)

        logger.debug(
            "DISEASES Knowledge and Experimental gene-disease data is merged"
        )

        # MERGE DISGENET GENE-DISEASE AND VARIANT-DISEASE ASSOCIATION DATA
        disgenet_df = pd.merge(
            disgenet_gda_df,
            disgenet_vda_df,
            how="outer",
            on=["gene_id", "disease_id"],
        )

        disgenet_df["source"] = disgenet_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        disgenet_df.drop(columns=["source_x", "source_y"], inplace=True)

        logger.debug(
            "Disgenet gene-disease and Disgenet variant-disease data is merged"
        )

        # MERGE OPENTARGETS AND DISEASES DATA
        merged_df = opentargets_df.merge(
            diseases_df, how="outer", on=["gene_id", "disease_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        logger.debug("Opentargets and DISEASES gene-disease data is merged")

        # MERGE OPENTARGETS+DISEASES AND KEGG DATA
        merged_df = merged_df.merge(
            kegg_df, how="outer", on=["gene_id", "disease_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        logger.debug(
            "Opentargets+DISEASES and KEGG gene-disease data is merged"
        )

        # MERGE OPENTARGETS+DISEASES+KEGG AND CLINVAR DATA
        merged_df = merged_df.merge(
            clinvar_df, how="outer", on=["gene_id", "disease_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        logger.debug(
            "Opentargets+DISEASES+KEGG and Clinvar gene-disease data is merged"
        )

        # MERGE OPENTARGETS+DISEASES+KEGG+CLINVAR AND HUMSAVAR DATA
        merged_df = merged_df.merge(
            humsavar_df, how="outer", on=["gene_id", "disease_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        # merge variant_source column
        merged_df["variant_source"] = merged_df[
            ["variant_source_x", "variant_source_y"]
        ].apply(self.merge_source_column, axis=1)

        merged_df.drop(
            columns=["variant_source_x", "variant_source_y"], inplace=True
        )

        # merge dbsnp_id column
        merged_df["dbsnp_id"] = merged_df[["dbsnp_id_x", "dbsnp_id_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["dbsnp_id_x", "dbsnp_id_y"], inplace=True)

        merged_df.replace("", np.nan, inplace=True)

        logger.debug(
            "Opentargets+DISEASES+KEGG+Clinvar and Humsavar gene-disease data is merged"
        )

        # MERGE OPENTARGETS+DISEASES+KEGG+CLINVAR+HUMSAVAR AND DISGENET DATA
        merged_df = merged_df.merge(
            disgenet_df, how="outer", on=["gene_id", "disease_id"]
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        # merge variant_source column
        merged_df["variant_source"] = merged_df[
            ["variant_source_x", "variant_source_y"]
        ].apply(self.merge_source_column, axis=1)

        merged_df.drop(
            columns=["variant_source_x", "variant_source_y"], inplace=True
        )

        # merge dbsnp_id column
        merged_df["dbsnp_id"] = merged_df[["dbsnp_id_x", "dbsnp_id_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["dbsnp_id_x", "dbsnp_id_y"], inplace=True)

        merged_df.replace("", np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"Gene-disease edge data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        # write gene-disease edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Gene_to_disease_edge.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Gene_to_disease_edge.csv"
                )

            merged_df.to_csv(full_path, index=False)
            logger.info(f"Gene-Disease edge data is written: {full_path}")

        return merged_df

    def merge_disease_disease_edge_data(self) -> pd.DataFrame:
        disgenet_dda_gene_df, disgenet_dda_variant_df = (
            self.process_disgenet_disease_disease()
        )

        logger.debug("Started merging disease-disease data")
        t0 = time()

        merged_df = disgenet_dda_gene_df.merge(
            disgenet_dda_variant_df,
            how="outer",
            on=["disease_id1", "disease_id2"],
        )

        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )

        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        merged_df = merged_df[
            merged_df["disease_id1"].ne(merged_df["disease_id2"])
        ]

        t1 = time()
        logger.info(
            f"Disease-disease edge data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        # write disease-disease association edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Disease_to_disease_association_edge.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Disease_to_disease_association_edge.csv"
                )

            merged_df.to_csv(full_path, index=False)
            logger.info(
                f"Disease-Disease association edge data is written: {full_path}"
            )

        return merged_df

    @validate_call
    def get_nodes(self, label: str = "disease") -> list[tuple]:
        if not hasattr(self, "mondo"):
            self.download_mondo_data()

        logger.debug("Started writing Disease nodes")

        node_list = []

        xref_dbs = list(
            set(DiseaseNodeField.get_xrefs()).intersection(
                self.disease_node_fields
            )
        )

        for index, term in enumerate(tqdm(self.mondo)):
            if not term.is_obsolete and term.obo_id and "MONDO" in term.obo_id:
                disease_id = self.add_prefix_to_id(
                    prefix="MONDO", identifier=term.obo_id
                )
                props = {}

                if (
                    DiseaseNodeField.NAME.value in self.disease_node_fields
                    and term.label
                ):
                    props[DiseaseNodeField.NAME.value] = term.label.replace(
                        "'", "^"
                    ).replace("|", ",")

                if (
                    DiseaseNodeField.SYNONYMS.value in self.disease_node_fields
                    and term.obo_synonym
                ):
                    synonym_set = {
                        syn["name"].replace("'", "^").replace("|", ",")
                        for syn in term.obo_synonym
                    }
                    props[DiseaseNodeField.SYNONYMS.value] = list(synonym_set)

                if xref_dbs and term.obo_xref:
                    for xref in term.obo_xref:
                        if xref["database"] in xref_dbs:
                            props[xref["database"].lower()] = xref["id"]

                if DiseaseNodeField.DOC2VEC_EMBEDDING.value in self.disease_node_fields and self.disease_id_to_doc2vec_embedding.get(term.obo_id) is not None:
                    props[DiseaseNodeField.DOC2VEC_EMBEDDING.value] = [str(emb) for emb in self.disease_id_to_doc2vec_embedding[term.obo_id]]


                node_list.append((disease_id, label, props))

                if self.early_stopping and index == self.early_stopping:
                    break

        # write node data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Disease.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Disease.csv")

            df_list = [
                {"disease_id": _id} | props for _id, _, props in node_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Disease node data is written: {full_path}")

        return node_list

    @validate_call
    def get_edges(self, 
                  disease_comorbidity_label: str = "disease_is_comorbid_with_disease",
                  disease_to_disase_label: str = "disease_is_associated_with_disease",
                  disease_to_drug_label: str = "disease_is_treated_by_drug",
                  gene_disease_label: str = "gene_is_related_to_disease",
                  mondo_hierarchy_label: str = "disease_is_a_disease",
                  organism_to_disease_label: str = "organism_causes_disease"
                  ) -> list[tuple]:
        edge_list = []

        if DiseaseEdgeType.DISEASE_COMOBORDITIY in self.edge_types:
            edge_list.extend(self.get_disease_comorbidity_edges(disease_comorbidity_label))

        if DiseaseEdgeType.DISEASE_TO_DISEASE in self.edge_types:
            edge_list.extend(self.get_disease_disease_edges(disease_to_disase_label))

        if DiseaseEdgeType.DISEASE_TO_DRUG in self.edge_types:
            edge_list.extend(self.get_disease_drug_edges(disease_to_drug_label))

        if DiseaseEdgeType.GENE_TO_DISEASE in self.edge_types:
            edge_list.extend(self.get_gene_disease_edges(gene_disease_label))

        if DiseaseEdgeType.MONDO_HIERARCHICAL_RELATIONS in self.edge_types:
            edge_list.extend(self.get_mondo_hiererchical_edges(mondo_hierarchy_label))

        if DiseaseEdgeType.ORGANISM_TO_DISEASE in self.edge_types:
            edge_list.extend(self.get_organism_disease_edges(organism_to_disease_label))

        return edge_list

    @validate_call
    def get_mondo_hiererchical_edges(
        self, label: str = "disease_is_a_disease"
    ) -> list[tuple]:
        if not hasattr(self, "mondo_hierarchical_relations"):
            self.download_mondo_data()

        logger.debug("Started writing Mondo hiererchical edges")

        edge_list = []

        counter = 0
        for source, target_list in tqdm(
            self.mondo_hierarchical_relations.items()
        ):
            source_id = self.add_prefix_to_id(prefix="MONDO", identifier=source)
            for target in target_list:
                target_id = self.add_prefix_to_id(
                    prefix="MONDO", identifier=target
                )
                edge_list.append((None, source_id, target_id, label, {}))

                counter += 1

            if self.early_stopping and counter >= self.early_stopping:
                break

        # write hiererchical edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Disease_hiererchical_edges.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Disease_hiererchical_edges.csv"
                )

            df_list = [
                {"child_id": child, "parent_id": parent, "label": label}
                for _, child, parent, label, _ in edge_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Mondo hiererchical edge data is written: {full_path}")

        return edge_list

    @validate_call
    def get_organism_disease_edges(
        self, label: str = "organism_causes_disease"
    ) -> list[tuple]:

        if not hasattr(self, "pathopheno_organism_disease_int"):
            self.download_pathophenodb_data()

        if not hasattr(self, "mondo_mappings"):
            self.prepare_mappings()

        logger.debug("Started writing organism-disease edges")

        edge_list = []

        for index, interaction in enumerate(
            tqdm(self.pathopheno_organism_disease_int)
        ):
            if (
                interaction.evidence == "manual assertion"
                and self.mondo_mappings["DOID"].get(
                    interaction.disease_id.split(":")[1], None
                )
            ):

                disease_id = self.add_prefix_to_id(
                    prefix="MONDO",
                    identifier=self.mondo_mappings["DOID"].get(
                        interaction.disease_id.split(":")[1]
                    ),
                )
                organism_id = self.add_prefix_to_id(
                    prefix="ncbitaxon", identifier=interaction.pathogen_taxid
                )

                edge_list.append((None, organism_id, disease_id, label, {}))

            if self.early_stopping and index == self.early_stopping:
                break

        # write hiererchical edge data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Organism_to_disease_edge.csv"
                )
            else:
                full_path = os.path.join(
                    os.getcwd(), "Organism_to_disease_edge.csv"
                )

            df_list = [
                {"organism_id": organism, "disease_id": disease, "label": label}
                for _, organism, disease, label, _ in edge_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Organism-Disease edge data is written: {full_path}")

        return edge_list

    @validate_call
    def get_disease_drug_edges(
        self, label: str = "disease_is_treated_by_drug"
    ) -> list[tuple]:

        disease_drug_edges_df = self.merge_disease_drug_edge_data()

        logger.debug("Started writing disease-drug edges")

        edge_list = []
        for index, row in tqdm(
            disease_drug_edges_df.iterrows(),
            total=disease_drug_edges_df.shape[0],
        ):
            _dict = row.to_dict()

            disease_id = self.add_prefix_to_id(
                prefix="MONDO", identifier=_dict["disease_id"]
            )
            drug_id = self.add_prefix_to_id(
                prefix="drugbank", identifier=_dict["drug_id"]
            )

            del _dict["disease_id"], _dict["drug_id"]

            props = {}
            for k, v in _dict.items():
                if k in self.disease_drug_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v

            edge_list.append((None, disease_id, drug_id, label, props))

            if self.early_stopping and index == self.early_stopping:
                break

        return edge_list

    @validate_call
    def get_gene_disease_edges(
        self, label: str = "gene_is_related_to_disease"
    ) -> list[tuple]:

        gene_disease_edges_df = self.merge_gene_disease_edge_data()

        logger.debug("Started writing gene-disease edges")

        edge_list = []

        for index, row in tqdm(
            gene_disease_edges_df.iterrows(),
            total=gene_disease_edges_df.shape[0],
        ):
            _dict = row.to_dict()

            gene_id = self.add_prefix_to_id(
                prefix="ncbigene", identifier=_dict["gene_id"]
            )
            disease_id = self.add_prefix_to_id(
                prefix="MONDO", identifier=_dict["disease_id"]
            )

            del _dict["disease_id"], _dict["gene_id"]

            props = {}
            for k, v in _dict.items():
                if k in self.gene_disease_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    elif k == "review_status":
                        props[k] = int(v)
                    else:
                        props[k] = v

            edge_list.append((None, gene_id, disease_id, label, props))

            if self.early_stopping and index == self.early_stopping:
                break

        return edge_list

    @validate_call
    def get_disease_disease_edges(
        self, label: str = "disease_is_associated_with_disease"
    ) -> list[tuple]:
        disease_disease_edges_df = self.merge_disease_disease_edge_data()

        logger.debug("Started writing disease-disease edges")

        edge_list = []
        for index, row in tqdm(
            disease_disease_edges_df.iterrows(),
            total=disease_disease_edges_df.shape[0],
        ):
            _dict = row.to_dict()

            disease_id1 = self.add_prefix_to_id(
                prefix="MONDO", identifier=_dict["disease_id1"]
            )
            disease_id2 = self.add_prefix_to_id(
                prefix="MONDO", identifier=_dict["disease_id2"]
            )

            del _dict["disease_id1"], _dict["disease_id2"]

            props = {}
            for k, v in _dict.items():
                if k in self.disease_disease_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v

            edge_list.append((None, disease_id1, disease_id2, label, props))

            if self.early_stopping and index == self.early_stopping:
                break

        return edge_list

    @validate_call
    def get_disease_comorbidity_edges(
        self, label: str = "disease_is_comorbid_with_disease"
    ) -> list[tuple]:
        comorbidity_edges_df = self.process_malacards_disease_comorbidity()

        logger.debug("Started writing disease comorbidity edges")

        edge_list = []

        props = {}
        for index, row in tqdm(
            comorbidity_edges_df.iterrows(), total=comorbidity_edges_df.shape[0]
        ):
            disease_id1 = self.add_prefix_to_id(
                prefix="MONDO", identifier=row["disease1"]
            )
            disease_id2 = self.add_prefix_to_id(
                prefix="MONDO", identifier=row["disease2"]
            )

            edge_list.append((None, disease_id1, disease_id2, label, props))

            if self.early_stopping and index == self.early_stopping:
                break

        return edge_list

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

    def map_disgenet_disease_id_to_mondo_id(
        self, disgenet_id, return_pandas_none=True
    ):

        disgenet_dbs_to_mondo_dbs = {
            "DO": "DOID",
            "EFO": "EFO",
            "HPO": "HP",
            "MSH": "MESH",
            "NCI": "NCIT",
            "ICD10CM": "ICD10CM",
            "OMIM": "OMIM",
        }

        diseaseid = self.mondo_mappings["UMLS"].get(disgenet_id)

        if not diseaseid and self.disgenet_id_mappings_dict.get(disgenet_id):
            if self.disgenet_id_mappings_dict.get(disgenet_id).get("MONDO"):
                diseaseid = "MONDO:" + self.disgenet_id_mappings_dict.get(
                    disgenet_id
                ).get("MONDO")
            else:
                map_dict = self.disgenet_id_mappings_dict.get(disgenet_id)
                for db, map_v in map_dict.items():
                    if self.mondo_mappings[disgenet_dbs_to_mondo_dbs[db]].get(
                        map_v
                    ):
                        diseaseid = self.mondo_mappings[
                            disgenet_dbs_to_mondo_dbs[db]
                        ].get(map_v)
                        break

        if diseaseid:
            return diseaseid
        elif return_pandas_none:
            return np.nan
        else:
            return None

    def ensembl_transcript_to_ensembl_gene(self, enst_ids) -> str | None:
        if not (
            enst_id_list := [
                e.split(" ")[0].split(".")[0] for e in enst_ids.split(";") if e
            ]
        ):
            return None
        if ensg_id_list := {
            list(mapping.map_name(_id, "enst_biomart", "ensg_biomart"))[0]
            for _id in enst_id_list
            if mapping.map_name(_id, "enst_biomart", "ensg_biomart")
        }:
            return list(ensg_id_list)[0]
        else:
            return None

    def ensembl_gene_to_ensembl_protein(self, ensg_id):
        return mapping.map_name(ensg_id, "ensg_biomart", "ensp_biomart")

    def set_edge_types(self, edge_types):
        self.edge_types = edge_types or list(DiseaseEdgeType)

    def set_node_and_edge_fields(
        self,
        disease_node_fields,
        gene_disease_edge_fields,
        disease_drug_edge_fields,
        disease_disease_edge_fields,
    ):
        if disease_node_fields:
            self.disease_node_fields = [
                field.value for field in disease_node_fields
            ]
        else:
            self.disease_node_fields = [
                field.value for field in DiseaseNodeField
            ]

        if gene_disease_edge_fields:
            self.gene_disease_edge_fields = [
                field.value for field in gene_disease_edge_fields
            ]
        else:
            self.gene_disease_edge_fields = [
                field.value for field in GENE_TO_DISEASE_INTERACTION_FIELD
            ]

        if disease_drug_edge_fields:
            self.disease_drug_edge_fields = [
                field.value for field in disease_drug_edge_fields
            ]
        else:
            self.disease_drug_edge_fields = [
                field.value for field in DISEASE_TO_DRUG_INTERACTION_FIELD
            ]

        if disease_disease_edge_fields:
            self.disease_disease_edge_fields = [
                field.value for field in disease_disease_edge_fields
            ]
        else:
            self.disease_disease_edge_fields = [
                field.value for field in DISEASE_TO_DISEASE_INTERACTION_FIELD
            ]
