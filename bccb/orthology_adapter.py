from __future__ import annotations

from pypath.share import curl, settings
from pypath.inputs import oma, uniprot, pharos
from pypath.utils import taxonomy

from typing import Union
from contextlib import ExitStack
from bioregistry import normalize_curie
from time import time
import os

import pandas as pd
import numpy as np

from enum import Enum, EnumMeta, IntEnum
from pydantic import BaseModel, DirectoryPath, validate_call

from biocypher._logger import logger
from tqdm import tqdm

logger.debug(f"Loading module {__name__}.")


class OrthologyEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class OrthologyEdgeField(Enum, metaclass=OrthologyEnumMeta):
    SOURCE = "source"
    RELATION_TYPE = "relation_type"
    OMA_ORTHOLOGY_SCORE = "oma_orthology_score"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class OMAOrganismField(IntEnum, metaclass=OrthologyEnumMeta):
    TAX_4932 = 4932  # s. cerevisiae
    TAX_10090 = 10090  # mouse
    TAX_3702 = 3702
    TAX_10116 = 10116  # rat
    TAX_559292 = 559292
    TAX_9913 = 9913  # cow
    TAX_1264690 = 1264690
    TAX_83333 = 83333
    TAX_6239 = 6239  # c. elegans
    TAX_1423 = 1423
    TAX_39947 = 39947
    TAX_44689 = 44689
    TAX_7227 = 7227  # drosophila
    TAX_8355 = 8355  # Xenopus laevis
    TAX_7955 = 7955  # zebrafish
    TAX_9031 = 9031  # chicken
    TAX_1773 = 1773
    TAX_9598 = 9598  # chimp - APES TOGETHER STRONG
    TAX_9544 = 9544  # Macaca - APES TOGETHER STRONG
    TAX_9595 = 9595  # GORILLA GORILLA GORILLA - APES TOGETHER STRONG
    TAX_9601 = 9601  # orangutan - APES TOGETHER STRONG

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class PharosOrganismField(Enum, metaclass=OrthologyEnumMeta):
    MOUSE = "Mouse"
    COW = "Cow"
    XENOPUS = "Xenopus"
    ZEBRAFISH = "Zebrafish"
    RAT = "Rat"
    C_ELEGANS = "C. elegans"
    S_CEREVISIA = "S.cerevisiae"
    CHICKEN = "Chicken"
    CHIMP = "Chimp"  # APES TOGETHER STRONG
    FRUITFLY = "Fruitfly"
    DOG = "Dog"
    MACAQUE = "Macaque"
    PIG = "Pig"
    HORSE = "Horse"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class OrthologyModel(BaseModel):
    edge_fields: Union[list[OrthologyEdgeField], None] = None
    oma_organisms: Union[list[OMAOrganismField], None] = None
    pharos_organisms: Union[list[PharosOrganismField], None] = None
    merge_with_pypath_taxids: bool = True
    add_prefix: bool = True
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None


class Orthology:
    """
    Class that downloads orthology data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(
        self,
        edge_fields: Union[list[OrthologyEdgeField], None] = None,
        oma_organisms: Union[list[OMAOrganismField], None] = None,
        pharos_organisms: Union[list[PharosOrganismField], None] = None,
        merge_with_pypath_taxids: bool = True,
        add_prefix: bool = True,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
    ):
        """
        Args:
            edge_fields: Gene-gene orthology edge fields that will be included in graph, if defined it must be values of elements from OrthologyEdgeField enum class (not the names)
            oma_organisms: list of taxanomy ids of organisms that will be compared against human to extract orthology relations, if defined it must be values of elements from OMAOrganismField enum class (not the names)
            pharos_organisms: list of taxanomy names of organisms that will be compared against human to extract orthology relations, if defined it must be values of elements from PharosOrganismField enum class (not the names)
            merge_with_pypath_taxids: Whether to merge `oma_organisms` list with list of pypath's tax ids
            add_prefix: if True, add prefix to database identifiers
            test_mode: if True, limits amount of output data
            export_csv: if True, export data as csv
            output_dir: Location of csv export, if not defined it will be current director. `export_csv` should be True for this arg.
        """

        model = OrthologyModel(
            edge_fields=edge_fields,
            oma_organisms=oma_organisms,
            pharos_organisms=pharos_organisms,
            merge_with_pypath_taxids=merge_with_pypath_taxids,
            add_prefix=add_prefix,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
        ).model_dump()

        self.add_prefix = model["add_prefix"]
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]

        # set edge fields
        self.set_edge_fields(edge_fields=model["edge_fields"])

        # set organism lists
        self.set_organism_lists(
            oma_organisms=model["oma_organisms"],
            pharos_organisms=model["pharos_organisms"],
            merge_with_pypath_taxids=model["merge_with_pypath_taxids"],
        )

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_orthology_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
    ):
        """
        Wrapper function to download orthology data from various databases using pypath.

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

            self.download_oma_data()
            self.download_pharos_data()

    def download_oma_data(self):
        """
        Downloads orthology data from OMA against human.
        """

        self.entry_name_to_uniprot = uniprot.uniprot_data(
            field="id", reviewed=True, organism="*"
        )
        self.entry_name_to_uniprot = {
            v: k for k, v in self.entry_name_to_uniprot.items()
        }

        uniprot_to_entrez = uniprot.uniprot_data(
            field="xref_geneid", reviewed=True, organism="*"
        )
        self.uniprot_to_entrez = {}
        for k, v in uniprot_to_entrez.items():
            self.uniprot_to_entrez[k] = v.strip(";").split(";")[0]

        logger.debug("Started downloading OMA orthology data")
        t0 = time()

        self.oma_orthology = []

        for t in tqdm(self.oma_organisms):
            tax_orthology = oma.oma_orthologs(organism_a=9606, organism_b=t)
            tax_orthology = [
                i
                for i in tax_orthology
                if i.a.id in self.entry_name_to_uniprot
                and i.b.id in self.entry_name_to_uniprot
            ]
            tax_orthology = [
                i
                for i in tax_orthology
                if self.uniprot_to_entrez.get(
                    self.entry_name_to_uniprot[i.a.id], None
                )
                and self.uniprot_to_entrez.get(
                    self.entry_name_to_uniprot[i.b.id], None
                )
            ]

            self.oma_orthology.extend(tax_orthology)
            logger.debug(f"Orthology data of tax id {t} is downloaded")

        t1 = time()
        logger.info(
            f"OMA orthology data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_oma_data(self) -> pd.DataFrame:
        """
        Processes orthology data from OMA.
        """
        if not hasattr(self, "oma_orthology"):
            self.download_oma_data()

        logger.debug("Started processing OMA orthology data")
        t0 = time()

        df_list = [
            (
                self.uniprot_to_entrez[
                    self.entry_name_to_uniprot[ortholog.a.id]
                ],
                self.uniprot_to_entrez[
                    self.entry_name_to_uniprot[ortholog.b.id]
                ],
                ortholog.rel_type,
                round(ortholog.score),
            )
            for ortholog in self.oma_orthology
        ]
        oma_orthology_df = pd.DataFrame(
            df_list,
            columns=[
                "entrez_a",
                "entrez_b",
                "relation_type",
                "oma_orthology_score",
            ],
        )

        oma_orthology_df["source"] = "OMA"

        oma_orthology_df.sort_values(
            by="oma_orthology_score", ascending=False, inplace=True
        )

        oma_orthology_df = oma_orthology_df[
            ~oma_orthology_df[["entrez_a", "entrez_b"]]
            .apply(frozenset, axis=1)
            .duplicated()
        ].reset_index(drop=True)

        t1 = time()
        logger.info(
            f"OMA orthology data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return oma_orthology_df

    def download_pharos_data(self):
        """
        Downloads orthology data from Pharos.
        """

        logger.debug("Started downloading Pharos orthology data")
        t0 = time()

        uniprot_to_entrez = uniprot.uniprot_data(
            field="xref_geneid", reviewed=True, organism="*"
        )
        self.entrez_to_uniprot = {}
        for k, v in uniprot_to_entrez.items():
            for entrez in v.strip(";").split(";"):
                if entrez:
                    self.entrez_to_uniprot[entrez] = k

        self.pharos_orthology_init = pharos.pharos_targets(orthologs=True)

        t1 = time()
        logger.info(
            f"Pharos orthology data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

    def process_pharos_data(self) -> pd.DataFrame:
        """
        Processes orthology data from Pharos.
        """
        if not hasattr(self, "pharos_orthology_init"):
            self.download_pharos_data()

        logger.debug("Started processing Pharos orthology data")
        t0 = time()

        df_list = []
        for protein in self.pharos_orthology_init:
            if protein["orthologs"]:
                for ortholog in protein["orthologs"]:
                    if (
                        ortholog["geneid"]
                        and str(ortholog["geneid"]) in self.entrez_to_uniprot
                        and str(ortholog["species"]) in self.pharos_organisms
                        and protein["uniprot"] in self.uniprot_to_entrez
                    ):

                        df_list.append(
                            (
                                self.uniprot_to_entrez[protein["uniprot"]],
                                str(ortholog["geneid"]),
                            )
                        )

        pharos_orthology_df = pd.DataFrame(
            df_list, columns=["entrez_a", "entrez_b"]
        )

        pharos_orthology_df["source"] = "Pharos"

        pharos_orthology_df = pharos_orthology_df[
            ~pharos_orthology_df[["entrez_a", "entrez_b"]]
            .apply(frozenset, axis=1)
            .duplicated()
        ].reset_index(drop=True)

        t1 = time()
        logger.info(
            f"Pharos orthology data is processed in {round((t1-t0) / 60, 2)} mins"
        )

        return pharos_orthology_df

    def merge_orthology_data(self) -> pd.DataFrame:
        """
        Merges orthology data from OMA and Pharos
        """
        logger.debug("Started merged OMA and Pharos orthology data")
        t0 = time()

        pharos_orthology_df = self.process_pharos_data()
        oma_orthology_df = self.process_oma_data()

        oma_plus_pharos_orthology_df = oma_orthology_df.merge(
            pharos_orthology_df, how="outer", on=["entrez_a", "entrez_b"]
        )

        oma_plus_pharos_orthology_df["source"] = oma_plus_pharos_orthology_df[
            ["source_x", "source_y"]
        ].apply(self.merge_source_column, axis=1)

        oma_plus_pharos_orthology_df.drop(
            columns=["source_x", "source_y"], inplace=True
        )

        t1 = time()
        logger.info(
            f"OMA and Pharos orthology data is merged in {round((t1-t0) / 60, 2)} mins"
        )

        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Orthology.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Orthology.csv")

            oma_plus_pharos_orthology_df.to_csv(full_path, index=False)
            logger.info(f"Orthology data is written: {full_path}")

        return oma_plus_pharos_orthology_df

    @validate_call
    def get_orthology_edges(
        self, label: str = "gene_is_orthologous_with_gene"
    ) -> list[tuple]:
        """
        Reformats orthology data to be ready for import into a BioCypher database.
        """
        all_orthology_df = self.merge_orthology_data()
        # define edge list
        edge_list = []

        logger.info("Preparing orthology edges.")

        for index, row in tqdm(
            all_orthology_df.iterrows(), total=all_orthology_df.shape[0]
        ):
            _dict = row.to_dict()

            source = self.add_prefix_to_id("ncbigene", _dict["entrez_a"])
            target = self.add_prefix_to_id("ncbigene", _dict["entrez_b"])

            del _dict["entrez_a"], _dict["entrez_b"]

            props = {}
            for k, v in _dict.items():
                if k in self.edge_fields and str(v) != "nan":
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

    def merge_source_column(self, element, joiner="|"):
        """
        Merges source columns' entries in selected dataframes
        """
        _list = []
        for e in list(element.dropna().values):
            if joiner in e:
                _list.extend(iter(e.split(joiner)))
            else:
                _list.append(e)

        return joiner.join(list(dict.fromkeys(_list).keys()))

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

    def set_edge_fields(self, edge_fields):
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field.value for field in OrthologyEdgeField]

    def set_organism_lists(
        self, oma_organisms, pharos_organisms, merge_with_pypath_taxids
    ):
        # define oma organisms
        if oma_organisms:
            self.oma_organisms = {field.value for field in oma_organisms}
        else:
            self.oma_organisms = {field.value for field in OMAOrganismField}

        # merge organisms with list of pypath taxids
        if merge_with_pypath_taxids:
            self.oma_organisms = self.oma_organisms.union(
                set(taxonomy.taxids.keys())
            )

        # define pharos organisms
        if pharos_organisms:
            self.pharos_organisms = [field.value for field in pharos_organisms]
        else:
            self.pharos_organisms = [field.value for field in PharosOrganismField]
