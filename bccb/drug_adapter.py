from __future__ import annotations

from pypath.share import curl, settings, common
from pypath.inputs import (
    drugbank,
    drugcentral,
    stitch,
    string,
    uniprot,
    dgidb,
    pharos,
    ctdbase,
    unichem,
    chembl,
    ddinter,
)
from . import kegg_local
from contextlib import ExitStack
from typing import Literal, Union, Optional
from bioregistry import normalize_curie
from tqdm import tqdm
from time import time

import os
import collections
import h5py
import gzip

import pandas as pd
import numpy as np

from biocypher._logger import logger

from enum import Enum, EnumMeta, auto
from pydantic import BaseModel, DirectoryPath, EmailStr, FilePath, validate_call

logger.debug(f"Loading module {__name__}.")


class DrugEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class DrugNodeField(Enum, metaclass=DrugEnumMeta):
    SMILES = "SMILES"
    INCHI = "InChI"
    INCHIKEY = "InChIKey"
    CAS = "cas_number"
    NAME = "name"
    GROUPS = "groups"
    GENERAL_REFERENCES = "general_references"
    ATC_CODES = "atc_codes"
    ZINC = "zinc"
    CHEMBL = "chembl"
    BINDINGDB = "bindingdb"
    CLINICALTRIALS = "clinicaltrials"
    CHEBI = "chebi"
    PUBCHEM = "pubchem"
    KEGG_DRUG = "KEGG Drug"
    RXCUI = "RxCUI"
    PHARMGKB = "PharmGKB"
    PDB = "PDB"
    DRUGCENTRAL = "Drugcentral"

    # read from h5 file
    SELFORMER_EMBEDDING = "selformer_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def get_drugbank_properties(cls):
        return [cls.SMILES.value, cls.INCHI.value, cls.INCHIKEY.value]

    @classmethod
    def get_drugbank_external_fields(cls):
        return [
            cls.KEGG_DRUG.value,
            cls.RXCUI.value,
            cls.PHARMGKB.value,
            cls.PDB.value,
            cls.DRUGCENTRAL.value,
        ]

    @classmethod
    def get_unichem_mapping_fields(cls):
        return [
            cls.ZINC.value,
            cls.CHEMBL.value,
            cls.BINDINGDB.value,
            cls.CLINICALTRIALS.value,
            cls.CHEBI.value,
            cls.PUBCHEM.value,
        ]

    @classmethod
    def get_drugbank_core_fields(cls):
        return [
            cls.CAS.value,
            cls.NAME.value,
            cls.GROUPS.value,
            cls.GENERAL_REFERENCES.value,
            cls.ATC_CODES.value,
        ]


class DrugDTIEdgeField(Enum, metaclass=DrugEnumMeta):
    SOURCE = "source"
    MECHANISM_OF_ACTION_TYPE = "mechanism_of_action_type"
    MECHANISM_OF_ACTION = "mechanism_of_action"
    REFERENCES = "references"
    KNOWN_ACTION = "known_action"
    DGIDB_SCORE = "dgidb_score"
    ACTIVITY_VALUE = "activity_value"
    ACTIVITY_TYPE = "activity_type"
    PCHEMBL = "pchembl"
    CONFIDENCE_SCORE = "confidence_score"
    DISEASE_EFFICACY = "disease_efficacy"
    DIRECT_INTERACTION = "direct_interaction"
    STITCH_COMBINED_SCORE = "stitch_combined_score"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class DrugDDIEdgeField(Enum, metaclass=DrugEnumMeta):
    SOURCE = "source"
    RECOMMENDATION = "recommendation"
    INTERACTION_LEVEL = "interaction_level"
    INTERACTION_TYPE = "interaction_type"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class DrugDGIEdgeField(Enum, metaclass=DrugEnumMeta):
    SOURCE = "source"
    ACTIVITY_TYPE = "activity_type"
    REFERENCES = "references"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class DrugEdgeType(Enum, metaclass=DrugEnumMeta):
    DRUG_DRUG_INTERACTION = auto()
    DRUG_TARGET_INTERACTION = auto()
    DRUG_GENE_INTERACTION = auto()


class DrugModel(BaseModel):
    drugbank_user: EmailStr
    drugbank_passwd: str
    node_fields: Union[list[DrugNodeField], None] = None
    dti_edge_fields: Union[list[DrugDTIEdgeField], None] = None
    ddi_edge_fields: Union[list[DrugDDIEdgeField], None] = None
    dgi_edge_fields: Union[list[DrugDGIEdgeField], None] = None
    edge_types: Union[list[DrugEdgeType], None] = None
    add_prefix: bool = True
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None


class Drug:
    """
    WARNING: Drugbank, Drug Central and STITCH database download urls contain
        version number/date. Please update this urls (or request an update) in
        resources.urls module of the pypath library before using this script
        in order to access the latest versions of data.

    Class that downloads drug data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(
        self,
        drugbank_user: EmailStr,
        drugbank_passwd: str,
        node_fields: Optional[Union[list[DrugNodeField], None]] = None,
        dti_edge_fields: Optional[Union[list[DrugDTIEdgeField], None]] = None,
        ddi_edge_fields: Optional[Union[list[DrugDDIEdgeField], None]] = None,
        dgi_edge_fields: Optional[Union[list[DrugDGIEdgeField], None]] = None,
        edge_types: Optional[Union[list[DrugEdgeType], None]] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
        export_csv: Optional[bool] = False,
        output_dir: Optional[DirectoryPath | None] = None,
    ):
        """
        Args:
            drugbank_user: drugbank username
            drugbank_passwd: drugbank password
            node_fields: drug node fields that will be included in graph, if defined it must be values of elements from DrugNodeField enum class (not the names)
            dti_edge_fields: Drug-Target edge fields that will be included in graph, if defined it must be values of elements from DrugDTIEdgeField enum class (not the names)
            ddi_edge_fields: Drug-Drug edge fields that will be included in graph, if defined it must be values of elements from DrugDDIEdgeField enum class (not the names)
            dgi_edge_fields: Drug-Gene edge fields that will be included in graph, if defined it must be values of elements from DrugDGIEdgeField enum class (not the names)
            edge_types: list of edge types that will be included in graph, if defined it must be elements (not values of elements) from DrugEdgeType enum class
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of output data
            export_csv: if True, export data as csv
            output_dir: Location of csv export if `export_csv` is True, if not defined it will be current directory
        """

        model = DrugModel(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            node_fields=node_fields,
            dti_edge_fields=dti_edge_fields,
            ddi_edge_fields=ddi_edge_fields,
            dgi_edge_fields=dgi_edge_fields,
            edge_types=edge_types,
            add_prefix=add_prefix,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
        ).model_dump()

        self.user = model["drugbank_user"]
        self.passwd = model["drugbank_passwd"]
        self.add_prefix = model["add_prefix"]
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]

        self.swissprots = set(
            uniprot._all_uniprots(organism="*", swissprot=True)
        )

        # set node fields
        self.set_node_fields(node_fields=model["node_fields"])

        # set edge fields
        self.set_edge_fields(
            model["dti_edge_fields"],
            model["ddi_edge_fields"],
            model["dgi_edge_fields"],
        )

        # set edge types
        self.set_edge_types(model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_drug_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 6,
        selformer_embedding_path: FilePath | None = None
    ):
        """
        Wrapper function to download drug data from various databases using pypath.

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

            self.download_drugbank_node_data()

            if DrugNodeField.SELFORMER_EMBEDDING.value in self.node_fields:
                self.retrieve_selformer_embeddings(selformer_embedding_path)

            if DrugEdgeType.DRUG_TARGET_INTERACTION in self.edge_types:
                self.download_chembl_dti_data()
                self.download_drugbank_dti_data()
                self.download_pharos_dti_data()
                self.download_dgidb_dti_data()
                self.download_kegg_dti_data()
                self.download_stitch_dti_data()

            if DrugEdgeType.DRUG_DRUG_INTERACTION in self.edge_types:
                self.download_kegg_ddi_data()
                self.download_ddinter_ddi_data()

            if DrugEdgeType.DRUG_GENE_INTERACTION in self.edge_types:
                self.download_ctd_data()

            t1 = time()
            logger.info(
                f"All data is downloaded in {round((t1-t0) / 60, 2)} mins".upper()
            )

    def process_drug_data(self):
        """
        An encompassing function that include all data processing functions
        """

        t0 = time()

        self.process_drugbank_node_data()

        if DrugEdgeType.DRUG_TARGET_INTERACTION in self.edge_types:
            self.process_drugbank_dti_data()
            self.process_chembl_dti_data()
            self.process_pharos_dti_data()
            self.process_dgidb_dti_data()
            self.process_stitch_dti_data()
            self.process_kegg_dti_data()

        if DrugEdgeType.DRUG_DRUG_INTERACTION in self.edge_types:
            self.process_kegg_ddi_data()
            self.process_ddinter_ddi_data()

        if DrugEdgeType.DRUG_GENE_INTERACTION in self.edge_types:
            self.process_ctd_data()

        t1 = time()
        logger.info(
            f"All data is processed in {round((t1-t0) / 60, 2)} mins".upper()
        )

    def download_drugbank_node_data(self):

        # define  fields belong to drugbank or unichem
        fields_list = DrugNodeField.get_drugbank_core_fields()
        unichem_external_fields_list = (
            DrugNodeField.get_unichem_mapping_fields()
        )
        drugbank_external_fields_list = (
            DrugNodeField.get_drugbank_external_fields()
        )
        drugbank_property_fields_list = DrugNodeField.get_drugbank_properties()

        self.unichem_external_fields = []
        self.drugbank_external_fields = []
        fields = []
        self.drugbank_property_fields = []
        self.add_inchi = False
        self.add_inchikey = False

        # TO DO: Unflag drugbank properties. Make them more proper
        for f in self.node_fields:
            if f in drugbank_property_fields_list:
                self.drugbank_property_fields.append(f)
            elif f in fields_list:
                fields.append(f)
            elif f in unichem_external_fields_list:
                self.unichem_external_fields.append(f)
            elif f in drugbank_external_fields_list:
                self.drugbank_external_fields.append(f)

        logger.debug("Downloading Drugbank drug node data")
        t0 = time()

        # Drugbank Object
        self.drugbank_data = drugbank.DrugbankFull(
            user=self.user, passwd=self.passwd
        )

        # external ids
        self.drugbank_drugs_external_ids = (
            self.drugbank_data.drugbank_external_ids_full()
        )

        # inchi, inchikey and smiles
        if self.drugbank_property_fields:
            self.drugbank_properties = (
                self.drugbank_data.drugbank_properties_full()
            )

        # core properties
        self.drugbank_drugs_detailed = self.drugbank_data.drugbank_drugs_full(
            fields=fields
        )

        # combine all external id mappings to self.drug_mappings_dict
        self.get_external_database_mappings()

        t1 = time()
        logger.info(
            f"Drugbank drug node data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_drugbank_node_data(self) -> dict:

        if not hasattr(self, "drugbank_drugs_detailed"):
            self.download_drugbank_node_data()

        logger.debug("Processing Drugbank drug node data")
        t0 = time()

        drugbank_drugs = {}
        all_fields = (
            list(self.drugbank_drugs_detailed[0]._fields)
            + self.drugbank_external_fields
            + self.unichem_external_fields
        )

        for drug in self.drugbank_drugs_detailed:

            drugbank_id = drug.drugbank_id
            temp_dict = (
                drug._asdict() | self.drug_mappings_dict[drugbank_id]
                if drugbank_id in self.drug_mappings_dict
                else drug._asdict()
            )

            drugbank_drugs[drugbank_id] = {
                f: temp_dict.get(f, None) for f in all_fields
            }

            if self.drugbank_property_fields:
                for p in self.drugbank_property_fields:
                    drugbank_drugs[drugbank_id][p] = (
                        self.drugbank_properties.get(drugbank_id, {}).get(
                            p, None
                        )
                    )

            del drugbank_drugs[drugbank_id]["drugbank_id"]

        t1 = time()
        logger.info(
            f"Drugbank drug node data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        # write node data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Drug.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Drug.csv")

            drugs_df_list = []
            for drug, attributes in drugbank_drugs.items():
                props = {"drugbank_id": drug} | attributes
                drugs_df_list.append(props)

            df = pd.DataFrame.from_records(drugs_df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Drug node data is written: {full_path}")

        return drugbank_drugs
    
    @validate_call
    def retrieve_selformer_embeddings(self,
                                selformer_embedding_path: FilePath | None = None) -> None:
        
        logger.info("Retrieving SELFormer drug embeddings")

        self.drugbank_id_to_selformer_embedding = {}
        with h5py.File(selformer_embedding_path, "r") as f:
            for drug_id, embedding in tqdm(f.items(), total=len(f.keys())):
                self.drugbank_id_to_selformer_embedding[drug_id] = np.array(embedding).astype(np.float16)

    def download_drugbank_dti_data(self):

        # DTI edge data
        logger.debug("Downloading Drugbank DTI data")
        t0 = time()

        if not hasattr(self, "drugbank_data"):
            self.drugbank_data = drugbank.DrugbankFull(
                user=self.user, passwd=self.passwd
            )

        self.drugbank_dti = self.drugbank_data.drugbank_targets_full(
            fields=[
                "drugbank_id",
                "actions",
                "references",
                "known_action",
                "polypeptide",
            ]
        )

        t1 = time()
        logger.info(
            f"Drugbank DTI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def get_external_database_mappings(self):
        logger.debug("Createing external database mappings")

        if not hasattr(self, "unichem_external_fields"):
            self.unichem_external_fields = (
                DrugNodeField.get_unichem_mapping_fields()
            )

        if not hasattr(self, "drugbank_external_fields"):
            self.drugbank_external_fields = (
                DrugNodeField.get_drugbank_external_fields()
            )

        if not hasattr(self, "drugbank_drugs_external_ids"):
            if not hasattr(self, "drugbank_data"):
                self.drugbank_data = drugbank.DrugbankFull(
                    user=self.user, passwd=self.passwd
                )

            self.drugbank_drugs_external_ids = (
                self.drugbank_data.drugbank_external_ids_full()
            )

        if not hasattr(self, "drugbank_drugs_detailed"):
            if not hasattr(self, "drugbank_data"):
                self.drugbank_data = drugbank.DrugbankFull(
                    user=self.user, passwd=self.passwd
                )

            self.drugbank_drugs_detailed = (
                self.drugbank_data.drugbank_drugs_full(fields=["cas_number"])
            )

        # create dictionaries for every unichem external fields
        unichem_drugbank_to_zinc_mapping = unichem.unichem_mapping(
            "drugbank", "zinc"
        )
        unichem_drugbank_to_chembl_mapping = unichem.unichem_mapping(
            "chembl", "drugbank"
        )
        unichem_drugbank_to_chembl_mapping = {
            list(v)[0]: k for k, v in unichem_drugbank_to_chembl_mapping.items()
        }
        unichem_drugbank_to_bindingdb_mapping = unichem.unichem_mapping(
            "drugbank", "bindingdb"
        )
        unichem_drugbank_to_clinicaltrials_mapping = unichem.unichem_mapping(
            "drugbank", "clinicaltrials"
        )
        unichem_drugbank_to_chebi_mapping = unichem.unichem_mapping(
            "drugbank", "chebi"
        )
        unichem_drugbank_to_pubchem_mapping = unichem.unichem_mapping(
            "drugbank", "pubchem"
        )

        # store above dicts in a unified dict
        self.unichem_external_fields_dict = {
            "zinc": unichem_drugbank_to_zinc_mapping,
            "chembl": unichem_drugbank_to_chembl_mapping,
            "bindingdb": unichem_drugbank_to_bindingdb_mapping,
            "clinicaltrials": unichem_drugbank_to_clinicaltrials_mapping,
            "chebi": unichem_drugbank_to_chebi_mapping,
            "pubchem": unichem_drugbank_to_pubchem_mapping,
        }

        # arrange unichem dict for selected unichem fields
        for field in list(self.unichem_external_fields_dict.keys()):
            if field not in self.unichem_external_fields:
                del self.unichem_external_fields_dict[field]

        unichem_drugs_id_mappings = collections.defaultdict(dict)
        for k in self.drugbank_drugs_external_ids.keys():
            for (
                field_name,
                field_dict,
            ) in self.unichem_external_fields_dict.items():
                mapping = field_dict.get(k, None)
                if mapping and field_name != "chembl":
                    mapping = list(mapping)[0]

                unichem_drugs_id_mappings[k][field_name] = mapping

        # get drugcentral mappings
        self.cas_to_drugbank = {
            drug.cas_number: drug.drugbank_id
            for drug in self.drugbank_drugs_detailed
            if drug.cas_number
        }
        drugcentral_to_cas = drugcentral.drugcentral_mapping(
            id_type="drugcentral", target_id_type="cas"
        )

        # create drugbank-drugcentral, drugcentral-drugbank, chembl-drugbank mappings that will be used for the future processes
        chembl_to_drugbank = unichem.unichem_mapping("chembl", "drugbank")
        self.chembl_to_drugbank = {
            k: list(v)[0] for k, v in chembl_to_drugbank.items()
        }

        # create kegg-drugbank mapping
        self.kegg_to_drugbank = {
            v["KEGG Drug"]: k
            for k, v in self.drugbank_drugs_external_ids.items()
            if v.get("KEGG Drug", None)
        }

        self.drugcentral_to_drugbank = collections.defaultdict(list)
        self.drugbank_to_drugcentral = collections.defaultdict(None)
        for k, v in drugcentral_to_cas.items():
            if list(v)[0] in self.cas_to_drugbank and k:
                self.drugbank_to_drugcentral[
                    self.cas_to_drugbank[list(v)[0]]
                ] = k
                self.drugcentral_to_drugbank[k] = self.cas_to_drugbank[
                    list(v)[0]
                ]

                # add drugcentral id to drugbank_drugs_external_ids
                if (
                    self.cas_to_drugbank[list(v)[0]]
                    in self.drugbank_drugs_external_ids
                ):
                    self.drugbank_drugs_external_ids[
                        self.cas_to_drugbank[list(v)[0]]
                    ]["Drugcentral"] = k

        # create final external id mapping dict
        self.drug_mappings_dict = collections.defaultdict(dict)
        for k in self.drugbank_drugs_external_ids.keys():
            drugbank_mappings = {
                field: self.drugbank_drugs_external_ids[k].get(field, None)
                for field in self.drugbank_external_fields
            }
            unichem_mappings = unichem_drugs_id_mappings[k]
            self.drug_mappings_dict[k] = drugbank_mappings | unichem_mappings

    def process_drugbank_dti_data(self) -> pd.DataFrame:

        if not hasattr(self, "drugbank_dti"):
            self.download_drugbank_dti_data()

        def get_uniprot_ids(element):
            if element:
                if isinstance(element, tuple) and element[1] == "Swiss-Prot":
                    return [element[0]]
                elif isinstance(element, list):
                    _list = []
                    for x in element:
                        if x[1] == "Swiss-Prot":
                            _list.append(x[0])

                    return _list
            else:
                return None

        def aggregate_one_field(element, joiner="|"):
            return joiner.join(set(element)) if element else None

        logger.debug("Processing Drugbank DTI data")
        t0 = time()

        # create list instance for pandas dataframes
        df_list = []

        selected_fields = ["actions", "references", "known_action"]

        # process drugbank dti
        for dti in self.drugbank_dti:
            uniprot_ids = get_uniprot_ids(dti.polypeptide)

            if not uniprot_ids or not dti.drugbank_id:
                continue

            for _id in uniprot_ids:
                attributes = dti._asdict()

                _list = []

                _list.append(dti.drugbank_id)
                _list.append(_id)

                for field in selected_fields:

                    if field == "references":
                        if isinstance(attributes[field], list):
                            aggregated = aggregate_one_field(
                                [i for i in attributes[field] if i is not None]
                            )
                            _list.append(aggregated)
                        else:
                            _list.append(attributes[field])

                    elif field == "actions":
                        if isinstance(attributes[field], list):
                            _list.append(attributes[field][0])
                        else:
                            _list.append(attributes[field])

                    else:
                        _list.append(attributes[field])

                df_list.append(_list)

        # create pandas dataframe
        drugbank_dti_df = pd.DataFrame(
            df_list,
            columns=[
                "drugbank_id",
                "uniprot_id",
                "mechanism_of_action_type",
                "references",
                "known_action",
            ],
        )
        drugbank_dti_df.fillna(value=np.nan, inplace=True)
        drugbank_dti_df.replace("", np.nan, inplace=True)

        # add source
        drugbank_dti_df["source"] = "DrugBank"

        # Aggregate references field and remove duplicates
        drugbank_dti_df = (
            drugbank_dti_df.groupby(
                ["drugbank_id", "uniprot_id"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "drugbank_id": "first",
                    "uniprot_id": "first",
                    "mechanism_of_action_type": "first",
                    "references": self.aggregate_column_level,
                    "known_action": "first",
                    "source": "first",
                }
            )
            .replace("", np.nan)
        )

        drugbank_dti_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"Drugbank DTI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return drugbank_dti_df

    def download_dgidb_dti_data(self):
        """
        Wrapper function to download DGIdb DTI data using pypath

        It returns ChEMBL IDs as drug IDs and Entrez Gene IDs as protein IDs.
        """

        # edge data
        logger.debug("Downloading DGIdb DTI data")
        t0 = time()

        self.dgidb_dti = dgidb.dgidb_interactions()

        # map entrez gene ids to swissprot ids
        uniprot_to_entrez = uniprot.uniprot_data("xref_geneid", "*", True)
        self.entrez_to_uniprot = collections.defaultdict(list)
        for k, v in uniprot_to_entrez.items():
            for entrez_id in list(filter(None, v.split(";"))):
                self.entrez_to_uniprot[entrez_id].append(k)

        if not hasattr(self, "chembl_to_drugbank"):
            self.get_external_database_mappings()

        t1 = time()
        logger.info(
            f"Dgidb DTI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_dgidb_dti_data(self) -> pd.DataFrame:

        if not hasattr(self, "dgidb_dti"):
            self.download_dgidb_dti_data()

        logger.debug("Processing DGIdb DTI data")
        t0 = time()

        df_list = []

        for dti in self.dgidb_dti:
            if (
                dti.entrez
                and dti.drug_chembl
                and self.entrez_to_uniprot.get(dti.entrez)
                and self.chembl_to_drugbank.get(dti.drug_chembl.split(":")[1])
            ):
                pmid = "|".join(dti.pmid.split(",")) if dti.pmid else None
                df_list.append(
                    (
                        self.entrez_to_uniprot[dti.entrez][0],
                        self.chembl_to_drugbank[dti.drug_chembl.split(":")[1]],
                        dti.type,
                        dti.score,
                        pmid,
                    )
                )

        dgidb_dti_df = pd.DataFrame(
            df_list,
            columns=[
                "uniprot_id",
                "drugbank_id",
                "mechanism_of_action_type",
                "dgidb_score",
                "references",
            ],
        )

        dgidb_dti_df.fillna(value=np.nan, inplace=True)

        # add source
        dgidb_dti_df["source"] = "DGIdb"

        # sort by dgidb_score
        dgidb_dti_df = dgidb_dti_df.sort_values(
            by="dgidb_score", ignore_index=True, ascending=False
        )

        # remove pairs without drugbank ids
        dgidb_dti_df = dgidb_dti_df.dropna(
            subset="drugbank_id", axis=0
        ).reset_index(drop=True)

        # remove duplicates
        dgidb_dti_df = (
            dgidb_dti_df.groupby(
                ["drugbank_id", "uniprot_id"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "uniprot_id": "first",
                    "drugbank_id": "first",
                    "mechanism_of_action_type": self.get_middle_row,
                    "dgidb_score": "first",
                    "references": self.aggregate_column_level,
                    "source": "first",
                }
            )
            .replace("", np.nan)
        )

        dgidb_dti_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"Dgidb DTI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return dgidb_dti_df

    def download_kegg_dti_data(self, organism: str | list = "hsa"):
        """
        Wrapper function to download KEGG DTI and DDI data using pypath

        Args:
            organism: KEGG organism code. Default is "hsa" for human.
                    If None, it downloads DTI data for all organisms in KEGG.

        """

        # DTI
        if organism is None:
            organism = kegg_local._kegg_list("organism")

        organism = common.to_list(organism)

        logger.debug(
            f"Downloading KEGG DTI data for {len(organism)} organism(s)"
        )
        t0 = time()

        self.kegg_dti = set()

        for org_info in tqdm(organism):

            org = org_info if isinstance(org_info, str) else org_info[1]

            if org.endswith("gz"):
                continue

            organism_dti = kegg_local.drug_to_gene(org=org)

            if not organism_dti:
                logger.debug(
                    f"Skipped taxonomy abbreviation {org}. This is most likely due to the empty return in database."
                )
                continue

            logger.debug(
                f"Downloaded KEGG DTI data with taxonomy abbreviation {str(org)}"
            )

            organism_dti_count = 0
            for k, v in organism_dti.items():
                if self.kegg_to_drugbank.get(k, None) and not isinstance(
                    v, str
                ):
                    drugbank_id = self.kegg_to_drugbank[k]
                    for gene_entry in v.GeneEntries:
                        if isinstance(
                            gene_entry.uniprot_ids, (tuple, list, set)
                        ):
                            for uniprot_id in list(gene_entry.uniprot_ids):
                                if uniprot_id in self.swissprots:
                                    self.kegg_dti.add((drugbank_id, uniprot_id))
                                    organism_dti_count += 1
                        elif str(gene_entry.uniprot_ids) in self.swissprots:
                            self.kegg_dti.add(
                                (drugbank_id, gene_entry.uniprot_ids)
                            )
                            organism_dti_count += 1

            logger.debug(
                f"Total interaction count for {str(org)} is {organism_dti_count}"
            )

        self.kegg_dti = list(self.kegg_dti)

        t1 = time()
        logger.info(
            f"KEGG DTI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_kegg_dti_data(self) -> pd.DataFrame:
        if not hasattr(self, "kegg_dti"):
            self.download_kegg_dti_data()

        logger.debug("Processing KEGG DTI data")
        t0 = time()

        kegg_dti_df = pd.DataFrame(
            self.kegg_dti, columns=["drugbank_id", "uniprot_id"]
        )

        kegg_dti_df["source"] = "Kegg"

        t1 = time()
        logger.info(
            f"KEGG DTI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return kegg_dti_df

    def download_kegg_ddi_data(self, from_csv=False):

        # DDI
        logger.debug("Downloading KEGG DDI data, this may take around 12 hours")
        t0 = time()

        if from_csv:
            logger.info("Skipping to processing part")
        else:
            self.kegg_ddi_data = kegg_local.drug_to_drug()

        t1 = time()
        logger.info(
            f"KEGG DDI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_kegg_ddi_data(self, from_csv=False) -> pd.DataFrame:

        logger.debug("Processing KEGG DDI data")
        t0 = time()

        if from_csv:
            # even with cache it takes a lot of time to download kegg ddi data. To interfere that, processed
            # version of it saved in a csv file
            kegg_ddi_df = pd.read_csv("kegg_ddi_duplicate_removed_df.csv")
            kegg_ddi_df.rename(
                columns={"interaction_type": "recommendation"}, inplace=True
            )

        else:
            if not hasattr(self, "kegg_ddi_data"):
                self.download_kegg_ddi_data()

            if not hasattr(self, "kegg_to_drugbank"):
                self.get_external_database_mappings()

            kegg_ddi = set()

            for k, v in self.kegg_ddi_data.items():
                if self.kegg_to_drugbank.get(k, None) and not isinstance(
                    v, str
                ):
                    drug1_drugbank_id = self.kegg_to_drugbank[k]
                    if v.interactions:
                        for interaction in list(v.interactions):
                            if (
                                interaction.type == "drug"
                                and interaction.id
                                and self.kegg_to_drugbank.get(
                                    interaction.id, None
                                )
                            ):
                                drug2_drugbank_id = self.kegg_to_drugbank[
                                    interaction.id
                                ]

                                if interaction.contraindication:
                                    contraindication = "contraindication"
                                else:
                                    contraindication = ""

                                if interaction.precaution:
                                    precaution = "precaution"
                                else:
                                    precaution = ""

                                if contraindication and precaution:
                                    interaction_type = "|".join(
                                        [contraindication, precaution]
                                    )
                                else:
                                    interaction_type = (
                                        contraindication or precaution
                                    )

                                kegg_ddi.add(
                                    (
                                        drug1_drugbank_id,
                                        drug2_drugbank_id,
                                        interaction_type,
                                    )
                                )

            kegg_ddi_df = pd.DataFrame(
                list(kegg_ddi), columns=["drug1", "drug2", "recommendation"]
            )

            kegg_ddi_df.replace("", np.nan, inplace=True)

            kegg_ddi_df["source"] = "Kegg"

            kegg_ddi_df = kegg_ddi_df[
                ~kegg_ddi_df[["drug1", "drug2"]]
                .apply(frozenset, axis=1)
                .duplicated()
            ].reset_index(drop=True)

        t1 = time()
        logger.info(
            f"KEGG DDI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return kegg_ddi_df

    def download_ddinter_ddi_data(self):

        logger.debug("Downloading DDInter DDI data")
        t0 = time()

        ddinter_mappings = ddinter.ddinter_mappings()

        self.ddinter_to_drugbank = {
            mapping.ddinter: mapping.drugbank
            for mapping in ddinter_mappings
            if mapping.drugbank
        }

        self.ddinter_interactions = ddinter.ddinter_interactions()

        t1 = time()
        logger.info(
            f"DDInter DDI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_ddinter_ddi_data(self) -> pd.DataFrame:
        if not hasattr(self, "ddinter_interactions"):
            self.download_ddinter_ddi_data()

        logger.debug("Processing DDInter DDI data")
        t0 = time()

        df_list = []
        for interaction in self.ddinter_interactions:
            if self.ddinter_to_drugbank.get(
                interaction.drug1_id, None
            ) and self.ddinter_to_drugbank.get(interaction.drug2_id, None):
                if isinstance(interaction.level, tuple):
                    if len(interaction.level) > 1:
                        level = "|".join(list(interaction.level))
                    else:
                        level = list(interaction.level)[0]
                else:
                    level = interaction.level

                if isinstance(interaction.actions, tuple):
                    if len(interaction.actions) > 1:
                        actions = "|".join(list(interaction.actions))
                    else:
                        actions = list(interaction.actions)[0]
                else:
                    actions = interaction.actions

                df_list.append(
                    (
                        self.ddinter_to_drugbank[interaction.drug1_id],
                        self.ddinter_to_drugbank[interaction.drug2_id],
                        level,
                        actions,
                    )
                )

        ddinter_df = pd.DataFrame(
            df_list,
            columns=["drug1", "drug2", "interaction_level", "interaction_type"],
        )

        ddinter_df["source"] = "DDInter"

        ddinter_df = ddinter_df[
            ~ddinter_df[["drug1", "drug2"]]
            .apply(frozenset, axis=1)
            .duplicated()
        ].reset_index(drop=True)

        t1 = time()
        logger.info(
            f"DDInter DDI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return ddinter_df

    def download_pharos_dti_data(self):

        logger.debug("Downloading Pharos DTI data")
        t0 = time()

        self.pharos_dti = pharos.pharos_targets(ligands=True)

        t1 = time()
        logger.info(
            f"Pharos DTI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_pharos_dti_data(self) -> pd.DataFrame:
        if not hasattr(self, "pharos_dti"):
            self.download_pharos_dti_data()

        logger.debug("Processing Pharos DTI data")
        t0 = time()

        df_list = []
        for dti in self.pharos_dti:

            if dti["ligands"]:
                for lig in dti["ligands"]:
                    synonyms = dict(
                        [
                            (syn["name"], syn["value"].split(",")[0])
                            for syn in lig["synonyms"]
                            if syn["name"] in ["ChEMBL", "DrugCentral"]
                        ]
                    )

                    for act in lig["activities"]:
                        # add dtis that have activity value and type
                        if act["value"] and act["type"] and act["type"] != "-":

                            if act["pubs"]:
                                _list = [
                                    pub["pmid"]
                                    for pub in act["pubs"]
                                    if pub["__typename"] == "PubMed"
                                    and pub["pmid"]
                                ]
                                pubmeds = "|".join(set(_list))
                            else:
                                pubmeds = None

                            df_list.append(
                                (
                                    dti["uniprot"],
                                    act["type"],
                                    act["moa"],
                                    act["value"],
                                    pubmeds,
                                    self.drugcentral_to_drugbank.get(
                                        synonyms.get("DrugCentral"), None
                                    ),
                                )
                            )

        pharos_dti_df = pd.DataFrame(
            df_list,
            columns=[
                "uniprot_id",
                "activity_type",
                "mechanism_of_action_type",
                "pchembl",
                "references",
                "drugbank_id",
            ],
        )

        pharos_dti_df.fillna(value=np.nan, inplace=True)
        pharos_dti_df.replace("", np.nan, inplace=True)

        # add source
        pharos_dti_df["source"] = "Pharos"

        # Remove rows without drugbank id
        pharos_dti_df.dropna(axis=0, subset="drugbank_id", inplace=True)

        # Sort by activity_value
        pharos_dti_df.sort_values(by="pchembl", ignore_index=True, inplace=True)

        # For every drug-target pair instance the preprocess as follows:
        # - get middle row for activity_type and mechanism_of_action_type
        # - get median of activity_value
        # - aggregate all the references
        pharos_dti_df = (
            pharos_dti_df.groupby(
                ["uniprot_id", "drugbank_id"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "uniprot_id": "first",
                    "activity_type": self.get_middle_row,
                    "mechanism_of_action_type": self.get_middle_row,
                    "pchembl": self.get_median,
                    "references": self.aggregate_column_level,
                    "drugbank_id": "first",
                    "source": "first",
                }
            )
            .replace("", np.nan)
        )

        pharos_dti_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"Pharos DTI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return pharos_dti_df

    def download_chembl_dti_data(self) -> pd.DataFrame:

        logger.debug("Downloading Chembl DTI data")
        t0 = time()

        self.chembl_acts = chembl.chembl_activities(standard_relation="=")
        self.chembl_document_to_pubmed = chembl.chembl_documents()
        self.chembl_targets = chembl.chembl_targets()
        self.chembl_assays = chembl.chembl_assays()
        self.chembl_mechanisms = chembl.chembl_mechanisms()

        if not hasattr(self, "chembl_to_drugbank"):
            self.get_external_database_mappings()

        t1 = time()
        logger.info(
            f"Chembl DTI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_chembl_dti_data(self) -> pd.DataFrame:

        if not hasattr(self, "chembl_mechanisms"):
            self.download_chembl_dti_data()

        logger.debug("Processing Chembl DTI data")
        t0 = time()

        mechanism_dict = {i.chembl: i._asdict() for i in self.chembl_mechanisms}
        targets = [
            i for i in self.chembl_targets if i.accession in self.swissprots
        ]
        target_dict = {i.target_chembl_id: i.accession for i in targets}
        assay_dict = {
            i.assay_chembl_id: i
            for i in self.chembl_assays
            if i.assay_type == "B"
        }

        df_list = []

        for act in self.chembl_acts:
            # if activity is belong to a binding assay and if it have activity_type, activity_value and target uniprot id
            if act.assay_chembl in assay_dict and all(
                [
                    True if item else False
                    for item in [
                        act.standard_value,
                        act.standard_type,
                        target_dict.get(act.target_chembl),
                    ]
                ]
            ):

                df_list.append(
                    (
                        act.pchembl,
                        act.standard_value,
                        act.standard_type,
                        target_dict.get(act.target_chembl),
                        str(
                            self.chembl_document_to_pubmed.get(
                                act.document, None
                            )
                        ),
                        assay_dict[act.assay_chembl].confidence_score,
                        self.chembl_to_drugbank.get(act.chembl, None),
                        mechanism_dict.get(act.chembl, {}).get(
                            "action_type", None
                        ),
                        mechanism_dict.get(act.chembl, {}).get(
                            "direct_interaction", None
                        ),
                        mechanism_dict.get(act.chembl, {}).get(
                            "disease_efficacy", None
                        ),
                        mechanism_dict.get(act.chembl, {}).get(
                            "mechanism_of_action", None
                        ),
                    )
                )

        # create pandas dataframe
        chembl_cti_df = pd.DataFrame(
            df_list,
            columns=[
                "pchembl",
                "activity_value",
                "activity_type",
                "uniprot_id",
                "references",
                "confidence_score",
                "drugbank_id",
                "mechanism_of_action_type",
                "direct_interaction",
                "disease_efficacy",
                "mechanism_of_action",
            ],
        )

        chembl_cti_df.fillna(value=np.nan, inplace=True)
        chembl_cti_df.replace("None", np.nan, inplace=True)

        # add source
        chembl_cti_df["source"] = "ChEMBL"

        # SORT BY activity_value
        chembl_cti_df.sort_values(
            by="activity_value", ignore_index=True, inplace=True
        )

        chembl_dti_df = chembl_cti_df.dropna(
            subset=["drugbank_id"], axis=0
        ).reset_index(drop=True)
        
        chembl_dti_df = (
            chembl_dti_df.groupby(
                ["uniprot_id", "drugbank_id"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "pchembl": self.get_median,
                    "activity_value": self.get_median,
                    "activity_type": self.get_middle_row,
                    "uniprot_id": "first",
                    "references": self.aggregate_column_level,
                    "confidence_score": self.get_middle_row,
                    "drugbank_id": "first",
                    "mechanism_of_action_type": self.get_middle_row,
                    "direct_interaction": self.get_middle_row,
                    "disease_efficacy": self.get_middle_row,
                    "mechanism_of_action": self.get_middle_row,
                    "source": "first",
                }
            )
            .replace("", np.nan)
        )

        chembl_dti_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"Chembl DTI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return chembl_dti_df

    def download_ctd_data(self):

        logger.debug("Downloading CTD DGI data")
        t0 = time()

        self.ctd_dgi = ctdbase.ctdbase_relations(relation_type="chemical_gene")

        t1 = time()
        logger.info(
            f"CTD DGI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_ctd_data(self) -> pd.DataFrame:

        if not hasattr(self, "ctd_dgi"):
            self.download_ctd_data()

        if not hasattr(self, "cas_to_drugbank"):
            self.get_external_database_mappings()

        logger.debug("Processing CTD DGI data")
        t0 = time()

        df_list = []

        for cgi in self.ctd_dgi:
            if (
                cgi.GeneID
                and cgi.CasRN
                and any(
                    [
                        (
                            True
                            if item
                            in [
                                ["increases", "expression"],
                                ["decreases", "expression"],
                            ]
                            else False
                        )
                        for item in cgi.InteractionActions
                    ]
                )
                and self.cas_to_drugbank.get(cgi.CasRN, None)
            ):

                # if both (increases and decreases expression) of them occur in same InteractionActions field, don't add to list
                if (
                    len(
                        [
                            item
                            for item in cgi.InteractionActions
                            if item
                            in [
                                ["increases", "expression"],
                                ["decreases", "expression"],
                            ]
                        ]
                    )
                    > 1
                ):
                    continue

                if isinstance(cgi.PubMedIDs, list):
                    pmid = "|".join(cgi.PubMedIDs)
                else:
                    pmid = cgi.PubMedIDs

                interaction_actions = "_".join(
                    [
                        item
                        for item in cgi.InteractionActions
                        if item
                        in [
                            ["increases", "expression"],
                            ["decreases", "expression"],
                        ]
                    ][0]
                )

                df_list.append(
                    (
                        cgi.GeneID,
                        self.cas_to_drugbank.get(cgi.CasRN),
                        interaction_actions,
                        pmid,
                    )
                )

        ctd_cgi_df = pd.DataFrame(
            df_list,
            columns=["entrez_id", "drugbank_id", "action_type", "references"],
        )

        def detect_conflicting_action_type(element):
            # if decreases_expression and increases_expression occur in same drug-gene pair, the pair is probably a erroneous entry
            if len(set(element.dropna().values)) > 1:
                return np.nan
            else:
                return list(set(element.dropna().values))[0]

        ctd_cgi_df = (
            ctd_cgi_df.groupby(
                ["drugbank_id", "entrez_id"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "entrez_id": "first",
                    "drugbank_id": "first",
                    "action_type": detect_conflicting_action_type,
                    "references": "first",
                }
            )
            .replace("", np.nan)
        )

        ctd_cgi_df.dropna(subset="action_type", inplace=True)

        ctd_cgi_df["source"] = "CTD"

        t1 = time()
        logger.info(
            f"CTD DGI data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "DGI.csv")
            else:
                full_path = os.path.join(os.getcwd(), "DGI.csv")

            ctd_cgi_df.to_csv(full_path, index=False)
            logger.info(f"DGI data is written: {full_path}")

        return ctd_cgi_df

    @validate_call
    def download_stitch_dti_data(
        self,
        organism: str | list | None = "9606",
        score_threshold: (
            int
            | Literal[
                "highest_confidence",
                "high_confidence",
                "medium_confidence",
                "low_confidence",
            ]
        ) = "high_confidence",
        physical_interaction_score: bool = False,  # currently this arg doesnt work for organisms other than human. can be fixed if necessary.
    ):
        """
        Wrapper function to download STITCH DTI data using pypath

        Args:
            organism: Name or NCBI Taxonomy IDs of organisms
                If None, DTI for all organisms will be downloaded.
            score_threshold: Minimum required interaction score. user can use
                pre-defined confidence limits or can define a custom value.
            physical_interaction_score: If True, returns physical interaction scores of interactions.

        """

        if organism is None:
            organism = string.string_species()

        organism = common.to_list(organism)

        logger.debug("Started downloading STITCH data")
        t0 = time()

        # map string ids to swissprot ids
        uniprot_to_string = uniprot.uniprot_data("xref_string", "*", True)
        self.string_to_uniprot = collections.defaultdict(list)
        for k, v in uniprot_to_string.items():
            for string_id in list(filter(None, v.split(";"))):
                self.string_to_uniprot[string_id.split(".")[1]].append(k)

        # mapping to convert pubchem ids to drugbank ids
        if self.unichem_external_fields_dict.get("pubchem", None):

            self.pubchem_to_drugbank = {}
            for k, v in self.unichem_external_fields_dict["pubchem"].items():
                if len(v) > 1:
                    for value in list(v):
                        self.pubchem_to_drugbank[value] = k
                else:
                    self.pubchem_to_drugbank[list(v)[0]] = k
        else:

            self.pubchem_to_drugbank = {}
            for k, v in unichem.unichem_mapping("drugbank", "pubchem").items():
                if len(v) > 1:
                    for value in list(v):
                        self.pubchem_to_drugbank[value] = k
                else:
                    self.pubchem_to_drugbank[list(v)[0]] = k

        self.stitch_ints = []

        for tax in tqdm(organism):
            if str(tax) in {"36329"}:
                continue

            try:
                organism_stitch_ints = [
                    i
                    for i in stitch.stitch_links_interactions(
                        ncbi_tax_id=int(tax),
                        score_threshold=score_threshold,
                        physical_interaction_score=physical_interaction_score,
                    )
                    if i.partner_b in self.string_to_uniprot
                    and i.partner_a in self.pubchem_to_drugbank
                ]  # filter with swissprot ids

                logger.debug(
                    f"Downloaded STITCH data with taxonomy id {str(tax)}, interaction count is {len(organism_stitch_ints)}"
                )

                if organism_stitch_ints:
                    self.stitch_ints.extend(organism_stitch_ints)

            except (TypeError, gzip.BadGzipFile) as e:  #'NoneType' object is not an iterator
                logger.debug(
                    f"Error: {e}. Skipped tax id {tax}. This is most likely due to the empty file in database."
                )

        t1 = time()
        logger.info(
            f"STITCH data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_stitch_dti_data(self) -> pd.DataFrame:

        if not hasattr(self, "stitch_ints"):
            self.download_stitch_dti_data()

        logger.debug("Started processing STITCH data")
        t0 = time()

        df_list = [
            (
                self.pubchem_to_drugbank[dti.partner_a],
                self.string_to_uniprot[dti.partner_b][0],
                dti.combined_score,
            )
            for dti in self.stitch_ints
        ]
        stitch_dti_df = pd.DataFrame(
            df_list,
            columns=["drugbank_id", "uniprot_id", "stitch_combined_score"],
        )

        stitch_dti_df.fillna(value=np.nan, inplace=True)

        # add source
        stitch_dti_df["source"] = "STITCH"

        # sort by stitch_combined_score
        stitch_dti_df.sort_values(
            by="stitch_combined_score",
            ignore_index=True,
            inplace=True,
            ascending=False,
        )

        # remove duplicates
        stitch_dti_df = (
            stitch_dti_df.groupby(
                ["drugbank_id", "uniprot_id"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "drugbank_id": "first",
                    "uniprot_id": "first",
                    "stitch_combined_score": self.get_median,
                    "source": "first",
                }
            )
            .replace("", np.nan)
        )

        stitch_dti_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"STITCH data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return stitch_dti_df

    def merge_all_dtis(self) -> pd.DataFrame:

        drugbank_dti_df = self.process_drugbank_dti_data()
        chembl_dti_df = self.process_chembl_dti_data()
        pharos_dti_df = self.process_pharos_dti_data()
        dgidb_dti_df = self.process_dgidb_dti_data()
        stitch_dti_df = self.process_stitch_dti_data()
        kegg_dti_df = self.process_kegg_dti_data()

        logger.debug("Started merging Drugbank and Chembl DTI data")
        t0 = time()

        # merge drugbank and chembl dti
        drugbank_plus_chembl_dti_df = drugbank_dti_df.merge(
            chembl_dti_df, how="outer", on=["uniprot_id", "drugbank_id"]
        )

        # merge references
        drugbank_plus_chembl_dti_df["references"] = drugbank_plus_chembl_dti_df[
            ["references_x", "references_y"]
        ].apply(self.aggregate_column_level, axis=1)

        # merge sources
        drugbank_plus_chembl_dti_df["source"] = drugbank_plus_chembl_dti_df[
            ["source_x", "source_y"]
        ].apply(self.merge_source_column, axis=1)

        # merge mechanism of action types
        drugbank_plus_chembl_dti_df[
            "mechanism_of_action_type"
        ] = drugbank_plus_chembl_dti_df[
            ["mechanism_of_action_type_x", "mechanism_of_action_type_y"]
        ].apply(
            lambda x: (
                str(list(x.dropna())[0]).lower() if list(x.dropna()) else np.nan
            ),
            axis=1,
        )

        # drop redundant columns
        drugbank_plus_chembl_dti_df.drop(
            columns=[
                "references_x",
                "references_y",
                "source_x",
                "source_y",
                "mechanism_of_action_type_x",
                "mechanism_of_action_type_y",
            ],
            inplace=True,
        )

        t1 = time()
        logger.info(
            f"Drugbank and Chembl DTI data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        logger.debug("Started merging Drugbank+Chembl and Pharos DTI data")
        t0 = time()

        # merge drugbank+chembl and pharos dti
        drugbank_plus_chembl_plus_pharos_dti_df = (
            drugbank_plus_chembl_dti_df.merge(
                pharos_dti_df, how="outer", on=["uniprot_id", "drugbank_id"]
            )
        )

        # merge references
        drugbank_plus_chembl_plus_pharos_dti_df["references"] = (
            drugbank_plus_chembl_plus_pharos_dti_df[
                ["references_x", "references_y"]
            ].apply(self.aggregate_column_level, axis=1)
        )

        # merge mechanism of action types
        drugbank_plus_chembl_plus_pharos_dti_df[
            "mechanism_of_action_type"
        ] = drugbank_plus_chembl_plus_pharos_dti_df[
            ["mechanism_of_action_type_x", "mechanism_of_action_type_y"]
        ].apply(
            lambda x: (
                str(x.dropna().tolist()[0]).lower()
                if len(x.dropna().tolist()) > 0
                else np.nan
            ),
            axis=1,
        )

        # merge pchembl
        drugbank_plus_chembl_plus_pharos_dti_df[
            "pchembl"
        ] = drugbank_plus_chembl_plus_pharos_dti_df[
            ["pchembl_x", "pchembl_y"]
        ].apply(
            lambda x: (
                x.dropna().tolist()[0]
                if len(x.dropna().tolist()) > 0
                else np.nan
            ),
            axis=1,
        )

        # merge activity type
        drugbank_plus_chembl_plus_pharos_dti_df[
            "activity_type"
        ] = drugbank_plus_chembl_plus_pharos_dti_df[
            ["activity_type_x", "activity_type_y"]
        ].apply(
            lambda x: (
                x.dropna().tolist()[0]
                if len(x.dropna().tolist()) > 0
                else np.nan
            ),
            axis=1,
        )

        # merge sources
        drugbank_plus_chembl_plus_pharos_dti_df["source"] = (
            drugbank_plus_chembl_plus_pharos_dti_df[
                ["source_x", "source_y"]
            ].apply(self.merge_source_column, axis=1)
        )

        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_dti_df.drop(
            columns=[
                "source_x",
                "source_y",
                "references_x",
                "references_y",
                "pchembl_x",
                "pchembl_y",
                "activity_type_x",
                "activity_type_y",
                "mechanism_of_action_type_x",
                "mechanism_of_action_type_y",
            ],
            inplace=True,
        )

        t1 = time()
        logger.info(
            f"Drugbank+Chembl and Pharos DTI data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        logger.debug(
            "Started merging Drugbank+Chembl+Pharos and Dgidb DTI data"
        )
        t0 = time()

        # merge drugbank+chembl+pharos and dgidb
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df = (
            drugbank_plus_chembl_plus_pharos_dti_df.merge(
                dgidb_dti_df, how="outer", on=["uniprot_id", "drugbank_id"]
            )
        )

        # merge references
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df["references"] = (
            drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df[
                ["references_x", "references_y"]
            ].apply(self.aggregate_column_level, axis=1)
        )

        # merge mechanism of action types
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df[
            "mechanism_of_action_type"
        ] = drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df[
            ["mechanism_of_action_type_x", "mechanism_of_action_type_y"]
        ].apply(
            lambda x: (
                x.dropna().tolist()[0]
                if len(x.dropna().tolist()) > 0
                else np.nan
            ),
            axis=1,
        )

        # merge sources
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df["source"] = (
            drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df[
                ["source_x", "source_y"]
            ].apply(self.merge_source_column, axis=1)
        )

        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df.drop(
            columns=[
                "source_x",
                "source_y",
                "references_x",
                "references_y",
                "mechanism_of_action_type_x",
                "mechanism_of_action_type_y",
            ],
            inplace=True,
        )

        t1 = time()
        logger.info(
            f"Drugbank+Chembl+Pharos and Dgidb DTI data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        logger.debug(
            "Started merging Drugbank+Chembl+Pharos+Dgidb and Stitch DTI data"
        )
        t0 = time()

        # merge drugbank+chembl+pharos+dgidb and stitch
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df = (
            drugbank_plus_chembl_plus_pharos_plus_dgidb_dti_df.merge(
                stitch_dti_df, how="outer", on=["uniprot_id", "drugbank_id"]
            )
        )

        # merge sources
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df[
            "source"
        ] = drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df[
            ["source_x", "source_y"]
        ].apply(
            self.merge_source_column, axis=1
        )

        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df.drop(
            columns=["source_x", "source_y"], inplace=True
        )

        t1 = time()
        logger.info(
            f"Drugbank+Chembl+Pharos+Dgidb and Stitch DTI data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        logger.debug(
            "Started merging Drugbank+Chembl+Pharos+Dgidb+Stitch and Kegg DTI data"
        )
        t0 = time()

        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df = drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_dti_df.merge(
            kegg_dti_df, how="outer", on=["uniprot_id", "drugbank_id"]
        )

        # merge sources
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df[
            "source"
        ] = drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df[
            ["source_x", "source_y"]
        ].apply(
            self.merge_source_column, axis=1
        )

        # drop redundant columns
        drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df.drop(
            columns=["source_x", "source_y"], inplace=True
        )

        t1 = time()
        logger.info(
            f"Drugbank+Chembl+Pharos+Dgidb+Stitch and KEGG DTI data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "DTI.csv")
            else:
                full_path = os.path.join(os.getcwd(), "DTI.csv")

            drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df.to_csv(
                full_path, index=False
            )
            logger.info(f"DTI data is written: {full_path}")

        # return final dataframe
        return drugbank_plus_chembl_plus_pharos_plus_dgidb_plus_stitch_plus_kegg_dti_df

    def merge_all_ddis(self) -> pd.DataFrame:

        kegg_ddi_df = self.process_kegg_ddi_data()
        ddinter_ddi_df = self.process_ddinter_ddi_data()

        logger.debug("Started merging KEGG and DDInter DDI data")
        t0 = time()

        # merge kegg and ddinter ddi data
        kegg_plus_ddinter_ddi_df = kegg_ddi_df.merge(
            ddinter_ddi_df, how="outer", on=["drug1", "drug2"]
        )

        # merge source columns
        kegg_plus_ddinter_ddi_df["source"] = kegg_plus_ddinter_ddi_df[
            ["source_x", "source_y"]
        ].apply(self.merge_source_column, axis=1)

        # drop redundant columns
        kegg_plus_ddinter_ddi_df.drop(
            columns=["source_x", "source_y"], inplace=True
        )

        t1 = time()
        logger.info(
            f"KEGG and DDInter DDI data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "DDI.csv")
            else:
                full_path = os.path.join(os.getcwd(), "DDI.csv")

            kegg_plus_ddinter_ddi_df.to_csv(full_path, index=False)
            logger.info(f"DDI data is written: {full_path}")

        # return final dataframe
        return kegg_plus_ddinter_ddi_df

    @validate_call
    def get_drug_nodes(self, label: str = "drug") -> list[tuple]:
        """
        Merges drug node information from different sources.
        """
        drugbank_drugs = self.process_drugbank_node_data()

        node_list = []

        logger.debug("Started writing drug nodes")

        counter = 0
        for k, v in tqdm(drugbank_drugs.items()):
            drug_id = self.add_prefix_to_id(prefix="drugbank", identifier=k)

            props = {}
            for prop_key, prop_value in v.items():
                if prop_value and prop_key in self.node_fields:
                    if isinstance(prop_value, str):
                        props[prop_key.replace(" ", "_").lower()] = (
                            prop_value.replace("'", "^")
                        )
                    else:
                        props[prop_key.replace(" ", "_").lower()] = prop_value

            if DrugNodeField.SELFORMER_EMBEDDING.value in self.node_fields and self.drugbank_id_to_selformer_embedding.get(k) is not None:
                props[DrugNodeField.SELFORMER_EMBEDDING.value] = [str(emb) for emb in self.drugbank_id_to_selformer_embedding[k]]


            node_list.append((drug_id, label, props))

            counter += 1
            if self.early_stopping and counter == self.early_stopping:
                break

        return node_list

    @validate_call
    def get_edges(self, 
                  dti_label: str = "drug_targets_protein",
                  ddi_label: str = "drug_interacts_with_drug") -> list[tuple]:
        
        edge_list = []

        if DrugEdgeType.DRUG_TARGET_INTERACTION in self.edge_types:
            edge_list.extend(self.get_dti_edges(dti_label))

        if DrugEdgeType.DRUG_DRUG_INTERACTION in self.edge_types:
            edge_list.extend(self.get_ddi_edges(ddi_label))

        if DrugEdgeType.DRUG_GENE_INTERACTION in self.edge_types:
            edge_list.extend(self.get_dgi_edges())

        return edge_list

    @validate_call
    def get_dti_edges(self, label: str = "drug_targets_protein") -> list[tuple]:

        dti_df = self.merge_all_dtis()

        logger.debug("Started writing DTI edges")

        edge_list = []
        for index, row in tqdm(dti_df.iterrows(), total=dti_df.shape[0]):

            _dict = row.to_dict()
            source = self.add_prefix_to_id(
                prefix="drugbank", identifier=_dict["drugbank_id"]
            )
            target = self.add_prefix_to_id(
                prefix="uniprot", identifier=_dict["uniprot_id"]
            )

            del _dict["drugbank_id"], _dict["uniprot_id"]

            props = {}
            for k, v in _dict.items():
                if k in self.dti_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ", "_").lower()] = v.replace(
                            "'", "^"
                        ).split("|")
                    else:
                        props[str(k).replace(" ", "_").lower()] = str(
                            v
                        ).replace("'", "^")

            edge_list.append((None, source, target, label, props))

            if self.early_stopping and (index + 1) == self.early_stopping:
                break

        return edge_list

    def get_dgi_edges(self) -> list[tuple]:

        dgi_df = self.process_ctd_data()

        logger.debug("Started writing DGI edges")

        edge_list = []

        for index, row in tqdm(dgi_df.iterrows(), total=dgi_df.shape[0]):
            _dict = row.to_dict()

            source = self.add_prefix_to_id(
                prefix="drugbank", identifier=_dict["drugbank_id"]
            )
            target = self.add_prefix_to_id(
                prefix="ncbigene", identifier=_dict["entrez_id"]
            )

            if _dict["action_type"] == "decreases_expression":
                label = "_".join(["drug", "downregulates", "gene"])
            else:
                label = "_".join(["drug", "upregulates", "gene"])

            del _dict["drugbank_id"], _dict["entrez_id"], _dict["action_type"]

            props = {}
            for k, v in _dict.items():
                if k in self.dgi_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ", "_").lower()] = v.replace(
                            "'", "^"
                        ).split("|")
                    else:
                        props[str(k).replace(" ", "_").lower()] = str(
                            v
                        ).replace("'", "^")

            edge_list.append((None, source, target, label, props))

            if self.early_stopping and (index + 1) == self.early_stopping:
                break

        return edge_list

    @validate_call
    def get_ddi_edges(
        self, label: str = "drug_interacts_with_drug"
    ) -> list[tuple]:

        ddi_df = self.merge_all_ddis()

        logger.debug("Started writing DGI edges")

        edge_list = []

        for index, row in tqdm(ddi_df.iterrows(), total=ddi_df.shape[0]):
            _dict = row.to_dict()

            source = self.add_prefix_to_id(
                prefix="drugbank", identifier=_dict["drug1"]
            )
            target = self.add_prefix_to_id(
                prefix="drugbank", identifier=_dict["drug2"]
            )

            del _dict["drug1"], _dict["drug2"]

            props = {}
            for k, v in _dict.items():
                if k in self.ddi_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ", "_").lower()] = v.replace(
                            "'", "^"
                        ).split("|")
                    else:
                        props[str(k).replace(" ", "_").lower()] = str(
                            v
                        ).replace("'", "^")

            edge_list.append((None, source, target, label, props))

            if self.early_stopping and (index + 1) == self.early_stopping:
                break

        return edge_list

    def set_node_fields(self, node_fields):
        if node_fields:
            self.node_fields = [field.value for field in node_fields]
        else:
            self.node_fields = [field.value for field in DrugNodeField]

    def set_edge_fields(
        self, dti_edge_fields, ddi_edge_fields, dgi_edge_fields
    ):
        self.dti_edge_fields = (
            [field.value for field in dti_edge_fields]
            if dti_edge_fields
            else [field.value for field in DrugDTIEdgeField]
        )

        self.ddi_edge_fields = (
            [field.value for field in ddi_edge_fields]
            if ddi_edge_fields
            else [field.value for field in DrugDDIEdgeField]
        )

        self.dgi_edge_fields = (
            [field.value for field in dgi_edge_fields]
            if dgi_edge_fields
            else [field.value for field in DrugDGIEdgeField]
        )

    def set_edge_types(self, edge_types):
        self.edge_types = edge_types or list(DrugEdgeType)

    @validate_call
    def add_prefix_to_id(
        self, prefix: str, identifier: str, sep: str = ":"
    ) -> str:
        """
        Adds prefix to ids
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)

        return identifier

    def aggregate_column_level(self, element, joiner="|"):
        _set = set()
        for e in set(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _set.add(i)
            else:
                _set.add(e)

        return joiner.join(_set) if _set else np.nan

    def get_median(self, element):
        element = element.astype(np.float32)
        return round(float(element.dropna().median()), 3)

    def get_middle_row(self, element):
        if len(list(element.index)) == 1:
            return element.values[0]
        elif not list(element.dropna().index):
            return np.nan
        elif len(list(element.dropna().index)) % 2 == 1:
            middle = len(list(element.dropna().index)) // 2
            return element.dropna().values[middle]
        else:
            middle = round((len(list(element.dropna().index)) / 2 + 0.00001))
            return element.dropna().values[middle]

    def merge_source_column(self, element, joiner="|"):
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                _list.extend(iter(e.split(joiner)))
            else:
                _list.append(e)

        return joiner.join(list(dict.fromkeys(_list).keys()))
