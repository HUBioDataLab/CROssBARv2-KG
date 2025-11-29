from __future__ import annotations

import os
import collections
import h5py
import gzip

from pypath.share import curl, settings, common
from pypath.inputs import chembl, stitch, uniprot, unichem, string
from contextlib import ExitStack
from bioregistry import normalize_curie

from typing import Literal, Union, Optional

import pandas as pd
import numpy as np

from time import time

from tqdm import tqdm
from enum import Enum, EnumMeta
from pydantic import BaseModel, DirectoryPath, FilePath, validate_call

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class CompoundEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class CompoundNodeField(Enum, metaclass=CompoundEnumMeta):
    TYPE = "type"
    FULL_MWT = "full_mwt"
    SPECIES = "species"
    HEAVY_ATOMS = "heavy_atoms"
    ALOGP = "alogp"
    INCHI = "std_inchi"
    INCHIKEY = "std_inchi_key"
    QED_SCORE = "qed_weighted"
    SMILES = "canonical_smiles"

    SELFORMER_EMBEDDING = "selformer_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class CompoundCTIEdgeField(Enum, metaclass=CompoundEnumMeta):
    SOURCE = "source"
    PCHEMBL = "pchembl"
    ACTIVITY_VALUE = "activity_value"
    ACTIVITY_TYPE = "activity_type"
    ASSAY_CHEMBL = "assay_chembl"
    REFERENCES = "references"
    CONFIDENCE_SCORE = "confidence_score"
    STITCH_COMBINED_SCORE = "stitch_combined_score"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class CompoundModel(BaseModel):
    node_fields: Union[list[CompoundNodeField], None] = None
    cti_edge_fields: Union[list[CompoundCTIEdgeField], None] = None
    test_mode: bool = False
    add_prefix: bool = True
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    stitch_organism: Optional[int | Literal["*"] | None] = None,


class Compound:
    """
    Class that downloads compound data using pypath and reformats it to be ready
    for import into the BioCypher.
    """

    def __init__(
        self,
        node_fields: Optional[Union[list[CompoundNodeField], None]] = None,
        cti_edge_fields: Optional[Union[list[CompoundCTIEdgeField], None]] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
        export_csv: Optional[bool] = False,
        output_dir: Optional[DirectoryPath | None] = None,
        stitch_organism: Optional[int | Literal["*"] | None] = None,
    ):
        """
        Initialize the Compound class.

        Args:
            node_fields: Compound node fields that will be included in graph, if defined it must be values of elements from CompoundNodeField enum class (not the names)
            cti_edge_fields: Compound-target edge fields that will be included in grah, if defined it must be values of elements from CompoundCTIEdgeField enum class (not the names)
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of output data
            export_csv: if True, export data as csv
            output_dir: Location of csv export if `export_csv` is True, if not defined it will be current directory
            stitch_organism: NCBI taxonomy ID of the organism to download compound-target interaction data for, if "*", download for all organisms
        """

        model = CompoundModel(
            node_fields=node_fields,
            cti_edge_fields=cti_edge_fields,
            test_mode=test_mode,
            add_prefix=add_prefix,
            export_csv=export_csv,
            output_dir=output_dir,
            stitch_organism=stitch_organism
        ).model_dump()

        self.stitch_organism = (
            "*" if model["stitch_organism"] in ("*", None) else model["stitch_organism"]
        )

        self.add_prefix = model["add_prefix"]
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]

        # set node fields
        self.set_node_fields(node_fields=model["node_fields"])

        # set cti edge fields
        self.set_edge_fields(cti_edge_fields=model["cti_edge_fields"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_compound_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
        selformer_embedding_path: FilePath | None = None,
    ):
        """
        Wrapper function to download compound data from Chembl database using pypath.

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

            self.download_chembl_data()

            if CompoundNodeField.SELFORMER_EMBEDDING.value in self.node_fields:
                self.retrieve_selformer_embeddings(
                    selformer_embedding_path=selformer_embedding_path
                )
            self.download_stitch_cti_data()

            t1 = time()
            logger.info(
                f"All data is downloaded in {round((t1-t0) / 60, 2)} mins".upper()
            )

    def retrieve_selformer_embeddings(self, 
                                     selformer_embedding_path: FilePath | None = None) -> None:
        """
        Downloads the selformer compound embedding.

        Args:
            selformer_embedding_path: Path to the selformer compound embedding.
        """
        logger.info("Retrieving SELFormer compound embeddings")

        activities_chembl = self.get_activite_compounds()

        self.chembl_id_to_selformer_embedding = {}
        with h5py.File(selformer_embedding_path, "r") as f:
            for compound_id, embedding in tqdm(f.items(), total=len(f.keys())):
                if compound_id in activities_chembl and compound_id not in self.chembl_to_drugbank:
                    self.chembl_id_to_selformer_embedding[compound_id] = np.array(embedding).astype(np.float16)


    def process_compound_data(self):
        t0 = time()

        self.process_chembl_cti_data()
        self.process_stitch_cti_data()

        t1 = time()
        logger.info(
            f"All data is processed in {round((t1-t0) / 60, 2)} mins".upper()
        )

    def download_chembl_data(self) -> None:

        t0 = time()
        logger.debug("Started downloading Chembl data")

        self.compounds = chembl.chembl_molecules()

        self.chembl_acts = chembl.chembl_activities(
            standard_relation="=",
        )

        self.document_to_pubmed = chembl.chembl_documents()

        chembl_targets = chembl.chembl_targets()

        # filter out TrEMBL targets
        swissprots = set(uniprot._all_uniprots(organism="*", swissprot=True))
        chembl_targets = [
            i for i in chembl_targets if i.accession in swissprots
        ]
        self.target_dict = {
            i.target_chembl_id: i.accession for i in chembl_targets
        }

        chembl_assays = chembl.chembl_assays()
        self.assay_dict = {
            i.assay_chembl_id: i for i in chembl_assays if i.assay_type == "B"
        }

        # chembl to drugbank mapping
        self.chembl_to_drugbank = unichem.unichem_mapping("chembl", "drugbank")
        self.chembl_to_drugbank = {
            k: list(v)[0] for k, v in self.chembl_to_drugbank.items()
        }

        t1 = time()
        logger.info(
            f"Chembl data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_chembl_cti_data(self) -> pd.DataFrame:

        if not hasattr(self, "chembl_acts"):
            self.download_chembl_data()

        logger.debug(
            "Started Chembl processing compound-target interaction data"
        )
        t0 = time()

        # define a list for creating dataframe
        df_list = []

        # filter activities
        for act in self.chembl_acts:
            if (
                act.assay_chembl in self.assay_dict
                and act.chembl not in self.chembl_to_drugbank
                and all(
                    [
                        True if item else False
                        for item in [
                            act.standard_value,
                            act.standard_type,
                            self.target_dict.get(act.target_chembl, None),
                        ]
                    ]
                )
            ):

                df_list.append(
                    (
                        act.chembl,
                        act.pchembl,
                        act.standard_value,
                        act.standard_type,
                        act.assay_chembl,
                        self.target_dict.get(act.target_chembl, None),
                        str(self.document_to_pubmed.get(act.document, None)),
                        self.assay_dict[act.assay_chembl].confidence_score,
                    )
                )

        # create dataframe
        chembl_cti_df = pd.DataFrame(
            df_list,
            columns=[
                "chembl",
                "pchembl",
                "activity_value",
                "activity_type",
                "assay_chembl",
                "uniprot_id",
                "references",
                "confidence_score",
            ],
        )

        chembl_cti_df.fillna(value=np.nan, inplace=True)
        chembl_cti_df.replace("None", np.nan, inplace=True)

        # add source
        chembl_cti_df["source"] = "ChEMBL"

        # sort by activity value
        chembl_cti_df.sort_values(
            by="activity_value", ignore_index=True, inplace=True
        )

        # multiple processing
        chembl_cti_df = (
            chembl_cti_df.groupby(
                ["uniprot_id", "chembl"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "chembl": "first",
                    "pchembl": self.get_median,
                    "activity_value": self.get_median,
                    "activity_type": self.get_middle_row,
                    "assay_chembl": self.aggregate_column_level,
                    "uniprot_id": "first",
                    "references": self.aggregate_column_level,
                    "confidence_score": self.get_middle_row,
                    "source": "first",
                }
            )
            .replace("", np.nan)
        )

        chembl_cti_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"Chembl data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return chembl_cti_df

    @validate_call
    def download_stitch_cti_data(
        self,
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
    ) -> None:

        if self.stitch_organism is None or self.stitch_organism == "*":
            organism = string.string_species()
        else:
            organism = self.stitch_organism

        organism = common.to_list(organism)

        logger.debug("Started downloading STITCH data")
        t0 = time()

        # map string ids to swissprot ids
        uniprot_to_string = uniprot.uniprot_data("xref_string", "*", True)
        self.string_to_uniprot = collections.defaultdict(list)
        for k, v in uniprot_to_string.items():
            for string_id in list(filter(None, v.split(";"))):
                self.string_to_uniprot[string_id.split(".")[1]].append(k)

        chembl_to_pubchem = unichem.unichem_mapping("chembl", "pubchem")
        self.pubchem_to_chembl = {}
        for k, v in chembl_to_pubchem.items():
            if len(v) > 1:
                for value in list(v):
                    self.pubchem_to_chembl[value] = k
            else:
                self.pubchem_to_chembl[list(v)[0]] = k

        drugbank_to_pubchem = unichem.unichem_mapping("drugbank", "pubchem")
        self.pubchem_to_drugbank = {}
        for k, v in drugbank_to_pubchem.items():
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
                    and i.partner_a in self.pubchem_to_chembl
                    and i.partner_a not in self.pubchem_to_drugbank
                ]  # filter with swissprot ids

                logger.debug(
                    f"Downloaded STITCH data with taxonomy id {str(tax)}, interaction count is {len(organism_stitch_ints)}"
                )

                if organism_stitch_ints:
                    self.stitch_ints.extend(organism_stitch_ints)

            except (TypeError, gzip.BadGzipFile) as e:  #'NoneType' object is not an iterator

                logger.debug(
                    f"Error: {e}. Skipped tax id {tax}. This is most likely due to the empty file in database. Check the database file."
                )

        t1 = time()
        logger.info(
            f"STITCH data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_stitch_cti_data(self) -> pd.DataFrame:

        if not hasattr(self, "stitch_ints"):
            self.download_stitch_cti_data()

        logger.debug("Started processing STITCH data")
        t0 = time()

        df_list = [
            (
                self.pubchem_to_chembl[cti.partner_a],
                self.string_to_uniprot[cti.partner_b][0],
                cti.combined_score,
            )
            for cti in self.stitch_ints
        ]
        stitch_cti_df = pd.DataFrame(
            df_list, columns=["chembl", "uniprot_id", "stitch_combined_score"]
        )

        stitch_cti_df.fillna(value=np.nan, inplace=True)

        # add source
        stitch_cti_df["source"] = "STITCH"

        # sort by stitch_combined_score
        stitch_cti_df.sort_values(
            by="stitch_combined_score",
            ignore_index=True,
            inplace=True,
            ascending=False,
        )

        stitch_cti_df = (
            stitch_cti_df.groupby(
                ["chembl", "uniprot_id"], sort=False, as_index=False
            )
            .aggregate(
                {
                    "chembl": "first",
                    "uniprot_id": "first",
                    "stitch_combined_score": self.get_median,
                    "source": "first",
                }
            )
            .replace("", np.nan)
        )

        stitch_cti_df.fillna(value=np.nan, inplace=True)

        t1 = time()
        logger.info(
            f"STITCH data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return stitch_cti_df

    def merge_all_ctis(self) -> pd.DataFrame:

        stitch_cti_df = self.process_stitch_cti_data()
        chembl_cti_df = self.process_chembl_cti_data()

        logger.debug(
            "Started merging Chembl and Stitch CTI (Compound-Target Interaction) data"
        )
        t0 = time()

        # merge chembl and stitch cti data
        chembl_plus_stitch_cti_df = chembl_cti_df.merge(
            stitch_cti_df, how="outer", on=["uniprot_id", "chembl"]
        )

        # merge source column
        chembl_plus_stitch_cti_df["source"] = chembl_plus_stitch_cti_df[
            ["source_x", "source_y"]
        ].apply(self.merge_source_column, axis=1)

        # drop redundant columns
        chembl_plus_stitch_cti_df.drop(
            columns=["source_x", "source_y"], inplace=True
        )

        t1 = time()
        logger.info(
            f"Chembl and Stitch CTI data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        # write CTI data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "CTI.csv")
            else:
                full_path = os.path.join(os.getcwd(), "CTI.csv")

            chembl_plus_stitch_cti_df.to_csv(full_path, index=False)
            logger.info(f"CTI edge data is written: {full_path}")

        return chembl_plus_stitch_cti_df

    def get_activite_compounds(self) -> set[str]:
        # filter out chembl ids whether they have activity or not
        activities_chembl = set()
        for act in self.chembl_acts:
            if act.assay_chembl in self.assay_dict and all(
                [
                    True if item else False
                    for item in [
                        act.standard_value,
                        act.standard_type,
                        self.target_dict.get(act.target_chembl, None),
                    ]
                ]
            ):

                activities_chembl.add(act.chembl)

        return activities_chembl

    @validate_call
    def get_compound_nodes(
        self,
        label: str = "compound",
        rename_node_fields: dict | None = {
            "std_inchi": "inchi",
            "std_inchi_key": "inchikey",
            "canonical_smiles": "smiles",
            "qed_weighted": "qed_score",
        },
    ) -> list[tuple]:
        """
        Reformats compound node data to be ready for import into the BioCypher.
        """
        if not hasattr(self, "compounds"):
            self.download_chembl_data()

        logger.debug("Creating compound nodes.")

        # get active compounds
        activities_chembl = self.get_activite_compounds()

        compound_nodes = []
        counter = 0
        for compound in tqdm(self.compounds):

            if (
                compound.structure_type == "MOL"
                and compound.chembl not in self.chembl_to_drugbank
                and compound.chembl in activities_chembl
            ):

                compound_id = self.add_prefix_to_id("chembl", compound.chembl)

                _dict = compound._asdict()
                props = {
                    k: value
                    for k, value in _dict.items()
                    if k in self.node_fields and value
                }

                if rename_node_fields:
                    props = {
                        (
                            rename_node_fields[k]
                            if k in rename_node_fields.keys()
                            else k
                        ): v
                        for k, v in props.items()
                    }

                if CompoundNodeField.SELFORMER_EMBEDDING.value in self.node_fields and self.chembl_id_to_selformer_embedding.get(compound.chembl) is not None:
                    props[CompoundNodeField.SELFORMER_EMBEDDING.value] = [str(emb) for emb in self.chembl_id_to_selformer_embedding[compound.chembl]]

                compound_nodes.append((compound_id, label, props))

                counter += 1

                if self.early_stopping and counter >= self.early_stopping:
                    break

        # write compound node data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Compound.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Compound.csv")

            df_list = [
                {"compound_id": compound_id} | props
                for compound_id, _, props in compound_nodes
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Compound node data is written: {full_path}")

        return compound_nodes

    @validate_call
    def get_cti_edges(self, label="compound_targets_protein") -> list[tuple]:
        """
        Reformats compound-target edge data to be ready for import into the BioCypher.
        """
        chembl_cti_df = self.merge_all_ctis()

        logger.debug("Creating compound-target edges.")

        cti_edge_list = []

        for index, row in tqdm(
            chembl_cti_df.iterrows(), total=chembl_cti_df.shape[0]
        ):

            _dict = row.to_dict()
            source = self.add_prefix_to_id("chembl", _dict["chembl"])
            target = self.add_prefix_to_id("uniprot", _dict["uniprot_id"])

            del _dict["chembl"], _dict["uniprot_id"]
            props = {}
            for k, v in _dict.items():
                if k in self.cti_edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[str(k).replace(" ", "_").lower()] = v.replace(
                            "'", "^"
                        ).split("|")
                    else:
                        props[str(k).replace(" ", "_").lower()] = str(
                            v
                        ).replace("'", "^")

            cti_edge_list.append((None, source, target, label, props))

            if self.early_stopping and index == self.early_stopping:
                break

        return cti_edge_list

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

    def aggregate_column_level(self, element, joiner="|"):
        import numpy as np

        _set = set()
        for e in set(element.dropna().values):
            if joiner in e:
                for i in e.split(joiner):
                    _set.add(i)
            else:
                _set.add(e)

        return joiner.join(_set) if _set else np.nan

    def merge_source_column(self, element, joiner="|"):
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                _list.extend(iter(e.split(joiner)))
            else:
                _list.append(e)

        return joiner.join(list(dict.fromkeys(_list).keys()))

    @validate_call
    def add_prefix_to_id(
        self, prefix: str, identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to ids
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)

        return identifier

    def set_node_fields(self, node_fields):
        if node_fields:
            self.node_fields = [field.value for field in node_fields]
        else:
            self.node_fields = [field.value for field in CompoundNodeField]

    def set_edge_fields(self, cti_edge_fields):
        if cti_edge_fields:
            self.cti_edge_fields = [field.value for field in cti_edge_fields]
        else:
            self.cti_edge_fields = [field.value for field in CompoundCTIEdgeField]
