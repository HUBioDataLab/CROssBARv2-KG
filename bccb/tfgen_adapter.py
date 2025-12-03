from __future__ import annotations
from pypath.share import curl, settings

from pypath.inputs import collectri, dorothea, trrust, uniprot
import pypath.utils.mapping as mapping

from contextlib import ExitStack
from typing import Union, Literal
from enum import Enum, EnumMeta, IntEnum
from biocypher._logger import logger
from pydantic import BaseModel, DirectoryPath, validate_call

import os

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time

import pandas as pd
import numpy as np


class TFGeneEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class TFGeneEdgeField(Enum, metaclass=TFGeneEnumMeta):
    SOURCE = "source"
    PUBMED_ID = "pubmed_id"
    TF_EFFECT = "tf_effect"
    DOROTHEA_CONFIDENCE_LEVEL = "dorothea_confidence_level"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class TFGeneOrganismField(IntEnum, metaclass=TFGeneEnumMeta):
    TAX_9606 = 9606
    TAX_10090 = 10090

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class TFGeneModel(BaseModel):
    edge_fields: Union[list[TFGeneEdgeField], None] = None
    organism: Union[list[TFGeneOrganismField], None] = None
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True


class TFGene:
    def __init__(
        self,
        edge_fields: Union[list[TFGeneEdgeField], None] = None,
        organism: Union[list[TFGeneOrganismField], None] = None,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
        add_prefix: bool = True,
    ):
        """
        Args:
            edge_fields: TF-Gene edge fields that will be included in graph, if defined it must be values of elements from TFGenEdgeField enum class (not the names)
            organism: list of taxanomy ids of organisms that will be included in the graph, if defined it must be values of elements from TFGenOrganismField enum class (not the names)
            test_mode: if True, limits amount of output data
            export_csv: if True, exports the TF-Gene edge data as a csv file
            output_dir: directory path to export the csv file, if None, uses current working directory.
            add_prefix: if True, add prefix to database identifiers.
        """

        model = TFGeneModel(
            edge_fields=edge_fields,
            organism=organism,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
            add_prefix=add_prefix,
        ).model_dump()

        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]
        self.add_prefix = model["add_prefix"]

        # set edge fields
        self.set_edge_fields(edge_fields=model["edge_fields"])

        # set organism
        self.set_organism(organism=model["organism"])

        # map numbers to corresponding tf effect
        self.effect_mapping = {0: "Unknown", 1: "Activation", -1: "Repression"}

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_tfgen_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
    ) -> None:
        """
        Wrapper function to download tf-gen relation data from various databases using pypath.
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

            logger.debug("Started downloading TF-Gene data")
            t0 = time()

            self.download_dorothea_data()
            self.download_collectri_data()
            self.download_trrust_data()

            t1 = time()
            logger.info(
                f"TF-Gene data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )

    def download_dorothea_data(self) -> None:

        logger.debug("Started downloading DoRothEA data")
        t0 = time()

        if 9606 in self.organism:
            self.dorothea_interactions = list(
                dorothea.dorothea_interactions(
                    organism=9606, levels={"A", "B", "C"}
                )
            )

        t1 = time()
        logger.info(
            f"DoRothEA data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def download_collectri_data(self) -> None:
        logger.debug("Started downloading CollecTRI data")
        t0 = time()

        if 9606 in self.organism:
            self.collectri_interactions = list(
                collectri.collectri_interactions()
            )

            self.uniprot_to_entrez = {
                k: v.strip(";").split(";")[0]
                for k, v in uniprot.uniprot_data(
                    "xref_geneid", 9606, True
                ).items()
            }

        t1 = time()
        logger.info(
            f"CollecTRI data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def download_trrust_data(self) -> None:
        logger.debug("Started downloading TRRUST data")
        t0 = time()

        self.trrust_interactions = []
        self.trrust_gene_symbol_to_entrez_id = {}
        if 9606 in self.organism:
            self.trrust_gene_symbol_to_entrez_id |= {
                entry["gene_symbol"]: entry["entrez_id"]
                for entry in trrust.scrape_human()
            }
            self.trrust_interactions.extend(trrust.trrust_human())
        if 10090 in self.organism:
            self.trrust_gene_symbol_to_entrez_id |= {
                entry["gene_symbol"]: entry["entrez_id"]
                for entry in trrust.scrape_mouse()
            }
            self.trrust_interactions.extend(trrust.trrust_mouse())

        t1 = time()
        logger.info(
            f"TRRUST data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_dorothea_tf_gene(self) -> pd.DataFrame:

        if not hasattr(self, "dorothea_interactions"):
            self.download_dorothea_data()

        logger.debug("Started processing DoRothEA tf-gen data")
        t0 = time()

        df_list = []
        for interaction in self.dorothea_interactions:
            tf = self.map_gene_symbol_to_entrez_id(interaction.tf)
            target = self.map_gene_symbol_to_entrez_id(interaction.target)
            if tf and target:
                tf = list(tf)[0]
                target = list(target)[0]
                effect = self.effect_mapping[interaction.effect]

                pubmed_id = None
                if interaction.pubmed:
                    pubmed_id = interaction.pubmed

                df_list.append(
                    (tf, target, pubmed_id, effect, interaction.level)
                )

        df = pd.DataFrame(
            df_list,
            columns=[
                "tf",
                "target",
                "pubmed_id",
                "tf_effect",
                "dorothea_confidence_level",
            ],
        )

        df.drop_duplicates(
            subset=["tf", "target"], ignore_index=True, inplace=True
        )

        df.fillna(value=np.nan, inplace=True)

        df["source"] = "DoRothEA"

        t1 = time()
        logger.info(
            f"DoRothEA tf-gen data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_collectri_tf_gene(self) -> pd.DataFrame:
        if not hasattr(self, "collectri_interactions"):
            self.download_collectri_data()

        logger.debug("Started processing CollecTRI tf-gen data")
        t0 = time()

        df_list = []
        for interaction in self.collectri_interactions:
            if str(interaction.tf).startswith("COMPLEX"):
                for tf in list(interaction.tf):
                    if self.uniprot_to_entrez.get(
                        tf
                    ) and self.uniprot_to_entrez.get(interaction.target):
                        tf_entrez = self.uniprot_to_entrez[tf]
                        target_entrez = self.uniprot_to_entrez[
                            interaction.target
                        ]
                        effect = self.effect_mapping[interaction.effect]

                        pubmed_id = None
                        if interaction.pubmed:
                            pubmed_id = interaction.pubmed.replace(";", "|")

                        df_list.append(
                            (tf_entrez, target_entrez, pubmed_id, effect)
                        )

            elif self.uniprot_to_entrez.get(
                interaction.tf
            ) and self.uniprot_to_entrez.get(interaction.target):
                tf_entrez = self.uniprot_to_entrez[interaction.tf]
                target_entrez = self.uniprot_to_entrez[interaction.target]
                effect = self.effect_mapping[interaction.effect]

                pubmed_id = None
                if interaction.pubmed:
                    pubmed_id = interaction.pubmed.replace(";", "|")

                df_list.append((tf_entrez, target_entrez, pubmed_id, effect))

        df = pd.DataFrame(
            df_list, columns=["tf", "target", "pubmed_id", "tf_effect"]
        )

        df.fillna(value=np.nan, inplace=True)

        df = df.groupby(["tf", "target"], sort=False, as_index=False).aggregate(
            {
                "tf": "first",
                "target": "first",
                "pubmed_id": self.merge_source_column,
                "tf_effect": self.find_conflicting_tf_effects,
            }
        )

        df.dropna(subset="tf_effect", inplace=True)

        df["source"] = "CollecTRI"

        t1 = time()
        logger.info(
            f"CollecTRI tf-gen data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def process_trrust_tf_gene(self) -> pd.DataFrame:
        if not hasattr(self, "trrust_interactions"):
            self.download_trrust_data()

        logger.debug("Started processing TRRUST tf-gen data")
        t0 = time()

        df_list = []
        for interaction in self.trrust_interactions:
            tf = interaction.source_genesymbol
            target = interaction.target_genesymbol
            effect = interaction.effect

            if self.trrust_gene_symbol_to_entrez_id.get(tf):
                tf_entrez = self.trrust_gene_symbol_to_entrez_id[tf]

                if self.trrust_gene_symbol_to_entrez_id.get(target):
                    target_entrez = self.trrust_gene_symbol_to_entrez_id[target]
                    df_list.append(
                        (
                            tf_entrez,
                            target_entrez,
                            effect,
                        )
                    )

        df = pd.DataFrame(df_list, columns=["tf", "target", "tf_effect"])

        df = df.groupby(["tf", "target"], sort=False, as_index=False).aggregate(
            {
                "tf": "first",
                "target": "first",
                "tf_effect": self.find_conflicting_tf_effects,
            }
        )
        df.dropna(subset="tf_effect", inplace=True)

        df["source"] = "TRRUST"

        t1 = time()
        logger.info(
            f"TRRUST tf-gen data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return df

    def merge_tf_gen_data(self) -> pd.DataFrame:
        trrust_df = self.process_trrust_tf_gene()

        dorothea_df = self.process_dorothea_tf_gene()

        collectri_df = self.process_collectri_tf_gene()

        logger.debug("Started merging tf-gen edge data")
        t0 = time()

        # merge dorothea and collectri
        merged_df = pd.merge(
            dorothea_df, collectri_df, how="outer", on=["tf", "target"]
        )
        merged_df.replace("", np.nan, inplace=True)

        # merge source column
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        # merge pubmed_id column
        merged_df["pubmed_id"] = merged_df[
            ["pubmed_id_x", "pubmed_id_y"]
        ].apply(self.merge_source_column, axis=1)
        merged_df.drop(columns=["pubmed_id_x", "pubmed_id_y"], inplace=True)

        # merge tf_effect column
        merged_df["tf_effect"] = merged_df[
            ["tf_effect_x", "tf_effect_y"]
        ].apply(self.find_conflicting_tf_effects, axis=1)
        merged_df.drop(columns=["tf_effect_x", "tf_effect_y"], inplace=True)
        merged_df.dropna(subset="tf_effect", inplace=True)

        # merge dorothea+collectri and trrust
        merged_df = merged_df.merge(trrust_df, how="outer", on=["tf", "target"])
        merged_df.replace("", np.nan, inplace=True)

        # merge source column
        merged_df["source"] = merged_df[["source_x", "source_y"]].apply(
            self.merge_source_column, axis=1
        )
        merged_df.drop(columns=["source_x", "source_y"], inplace=True)

        # merge tf_effect column
        merged_df["tf_effect"] = merged_df[
            ["tf_effect_x", "tf_effect_y"]
        ].apply(self.find_conflicting_tf_effects, axis=1)
        merged_df.drop(columns=["tf_effect_x", "tf_effect_y"], inplace=True)
        merged_df.dropna(subset="tf_effect", inplace=True)

        t1 = time()
        logger.info(
            f"Tf-gene edge data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Tf_gen_edges.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Tf_gen_edges.csv")

            merged_df.to_csv(full_path, index=False)
            logger.info(f"Tf-gen edge data is written: {full_path}")

        return merged_df

    @validate_call
    def get_edges(self, label: str = "gene_regulates_gene") -> list[tuple]:

        tf_gen_edges_df = self.merge_tf_gen_data()

        logger.info("Started writing tf-gen edges")

        edge_list = []
        for index, row in tqdm(
            tf_gen_edges_df.iterrows(), total=tf_gen_edges_df.shape[0]
        ):
            _dict = row.to_dict()

            tf_id = self.add_prefix_to_id(
                prefix="ncbigene", identifier=_dict["tf"]
            )
            target_id = self.add_prefix_to_id(
                prefix="ncbigene", identifier=_dict["target"]
            )

            del _dict["tf"], _dict["target"]

            props = {}
            for k, v in _dict.items():
                if k in self.edge_fields and str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        props[k] = v.split("|")
                    else:
                        props[k] = v

            edge_list.append((None, tf_id, target_id, label, props))

            if self.early_stopping and index + 1 == self.early_stopping:
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

    @validate_call
    def map_gene_symbol_to_entrez_id(self, gene_symbol: str):
        return mapping.map_name(gene_symbol, "genesymbol", "entrez")

    def merge_source_column(self, element, joiner="|"):

        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                _list.extend(iter(e.split(joiner)))
            else:
                _list.append(e)

        return joiner.join(list(dict.fromkeys(_list).keys()))

    def find_conflicting_tf_effects(self, element):
        tf_effects = set(list(element.dropna().values))

        if len(tf_effects) > 2:
            return np.nan
        elif len(tf_effects) == 1:
            return list(tf_effects)[0]
        elif tf_effects in [
            {"Unknown", "Activation"},
            {"Unknown", "Repression"},
        ]:
            _list = list(tf_effects)
            _list.remove("Unknown")
            return _list[0]
        else:
            return np.nan

    def set_edge_fields(self, edge_fields):
        if edge_fields:
            self.edge_fields = [field.value for field in edge_fields]
        else:
            self.edge_fields = [field.value for field in TFGeneEdgeField]

    def set_organism(self, organism):
        if organism:
            self.organism = [field.value for field in organism]
        else:
            self.organism = [field.value for field in TFGeneOrganismField]
