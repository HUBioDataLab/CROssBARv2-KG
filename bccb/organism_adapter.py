from __future__ import annotations
from pypath.share import curl, settings
from contextlib import ExitStack

from bioregistry import normalize_curie
from tqdm import tqdm
from time import time

import pandas as pd

import os

from enum import Enum, EnumMeta, auto
from typing import Union

from pypath.inputs import ncbi_taxonomy

from biocypher._logger import logger

from pydantic import BaseModel, DirectoryPath, validate_call

logger.debug(f"Loading module {__name__}.")


class OrganismEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class OrganismNodeField(Enum, metaclass=OrganismEnumMeta):
    NAME = "organism_name"
    COMMON_NAME = "common_name"
    RANK = "rank"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class OrganismEdgeType(Enum, metaclass=OrganismEnumMeta):
    ORGANISM_HIERARCHY = auto()


class OrganismModel(BaseModel):
    organism_node_fields: Union[list[OrganismNodeField], None] = None
    edge_types: Union[list[OrganismEdgeType], None] = None
    ranks: Union[set[str], None] = None
    include_no_rank: bool = True
    tax_ids: Union[set[str], None] = None
    test_mode: bool = False
    export_csv: bool = False
    output_dir: DirectoryPath | None = None
    add_prefix: bool = True


class Organism:
    """
    Organism taxon nodes and their hierarchy from NCBI Taxonomy.

    Args:
        organism_node_fields: node properties to include, all if None.
        edge_types: edge types to include, all if None.
        ranks: keep only taxa of these ranks (e.g. {"species", "genus"}),
            all ranks if None.
        include_no_rank: keep the taxa without a named rank ("no rank"
            and "clade"). If a taxon is excluded, its descendants are
            linked to their closest included ancestor. The root is always
            kept, and it is the parent of the taxa without any other
            included ancestor.
        tax_ids: keep only these taxa (e.g. the organisms other adapters
            link to, as "ncbitaxon:9606" or "9606") and all of their
            ancestors up to the root, regardless of rank. If given, ranks
            and include_no_rank are not used. If None, the full taxonomy
            is kept.
        test_mode: if True, only the first 100 taxa are processed.
        export_csv: if True, nodes and edges are also written to csv.
        output_dir: directory of the csv files, cwd if None.
        add_prefix: if True, "ncbitaxon:" prefix is added to the ids.
    """

    def __init__(
        self,
        organism_node_fields: Union[list[OrganismNodeField], None] = None,
        edge_types: Union[list[OrganismEdgeType], None] = None,
        ranks: Union[set[str], None] = None,
        include_no_rank: bool = True,
        tax_ids: Union[set[str], None] = None,
        test_mode: bool = False,
        export_csv: bool = False,
        output_dir: DirectoryPath | None = None,
        add_prefix: bool = True,
    ):

        model = OrganismModel(
            organism_node_fields=organism_node_fields,
            edge_types=edge_types,
            ranks=ranks,
            include_no_rank=include_no_rank,
            tax_ids=tax_ids,
            test_mode=test_mode,
            export_csv=export_csv,
            output_dir=output_dir,
            add_prefix=add_prefix,
        ).model_dump()

        self.ranks = model["ranks"]
        self.include_no_rank = model["include_no_rank"]
        self.tax_ids = model["tax_ids"]
        self.export_csv = model["export_csv"]
        self.output_dir = model["output_dir"]
        self.add_prefix = model["add_prefix"]

        # source, licence, and version fields
        # version is the release date of the taxdump, set after download
        self.data_source = "NCBI Taxonomy"
        self.data_licence = "NCBI Copyright and Data Reuse Policy"
        self.data_version = None

        # set node fields
        self.set_node_fields(organism_node_fields=model["organism_node_fields"])

        # set edge types
        self.set_edge_types(edge_types=model["edge_types"])

        # set early_stopping, if test_mode true
        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100

    @validate_call
    def download_organism_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
    ) -> None:
        """
        Wrapper function to download NCBI Taxonomy data using pypath.
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

            logger.debug("Started downloading NCBI Taxonomy data")
            t0 = time()

            tax_ids = None
            if self.tax_ids is not None:
                # "ncbitaxon:9606" or "9606"
                tax_ids = {
                    int(tax_id.rsplit(":", 1)[-1]) for tax_id in self.tax_ids
                }

            self.taxa, self.data_version = ncbi_taxonomy.ncbi_taxonomy(
                ranks=self.ranks,
                include_no_rank=self.include_no_rank,
                keep_root=True,
                tax_ids=tax_ids,
                release_date=True,
            )

            t1 = time()
            logger.info(
                f"NCBI Taxonomy data is downloaded in {round((t1-t0) / 60, 2)} mins"
            )
            logger.info(f"NCBI Taxonomy taxdump version: {self.data_version}")

            if tax_ids is not None:
                missing = sorted(tax_ids - self.taxa.keys())
                logger.info(
                    f"{len(tax_ids) - len(missing)} of the {len(tax_ids)} "
                    f"requested taxa and their ancestors give {len(self.taxa)} "
                    "organism nodes."
                )
                if missing:
                    logger.info(
                        f"{len(missing)} requested taxa are not in NCBI "
                        f"Taxonomy (e.g. merged or deleted): {missing[:10]}"
                    )

    @validate_call
    def get_nodes(self, label: str = "organism") -> list[tuple]:
        if not hasattr(self, "taxa"):
            self.download_organism_data()

        logger.info("Started writing organism nodes")

        node_list = []

        for index, taxon in tqdm(enumerate(self.taxa.values())):
            organism_id = self.add_prefix_to_id(
                prefix="ncbitaxon", identifier=str(taxon.tax_id)
            )

            props = {}
            if OrganismNodeField.NAME.value in self.organism_node_fields:
                props[OrganismNodeField.NAME.value] = self._clean(
                    taxon.scientific_name
                )

            common_name = taxon.genbank_common_name or taxon.common_name
            if (
                OrganismNodeField.COMMON_NAME.value in self.organism_node_fields
                and common_name
            ):
                props[OrganismNodeField.COMMON_NAME.value] = self._clean(
                    common_name
                )

            if OrganismNodeField.RANK.value in self.organism_node_fields:
                props[OrganismNodeField.RANK.value] = taxon.rank

            # source, licence, and version fields
            props["source"] = self.data_source
            props["licence"] = self.data_licence
            props["version"] = self.data_version

            node_list.append((organism_id, label, props))

            if self.early_stopping and index + 1 == self.early_stopping:
                break

        # write organism node data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "Organism.csv")
            else:
                full_path = os.path.join(os.getcwd(), "Organism.csv")

            df_list = [
                {"organism_id": _id} | props for _id, _, props in node_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Organism node data is written: {full_path}")

        return node_list

    def get_edges(self) -> list[tuple]:

        edge_list = []

        if OrganismEdgeType.ORGANISM_HIERARCHY in self.edge_types:
            edge_list.extend(self.get_organism_hierarchy_edges())

        return edge_list

    @validate_call
    def get_organism_hierarchy_edges(
        self, label: str = "organism_is_a_organism"
    ) -> list[tuple]:
        if not hasattr(self, "taxa"):
            self.download_organism_data()

        logger.info("Started writing organism hierarchical edges")

        edge_list = []

        for index, taxon in tqdm(enumerate(self.taxa.values())):
            # root has no parent
            if taxon.parent_tax_id is not None:
                child_id = self.add_prefix_to_id(
                    prefix="ncbitaxon", identifier=str(taxon.tax_id)
                )
                parent_id = self.add_prefix_to_id(
                    prefix="ncbitaxon", identifier=str(taxon.parent_tax_id)
                )

                edge_list.append((None, child_id, parent_id, label, {}))

            if self.early_stopping and index + 1 == self.early_stopping:
                break

        # write organism hierarchy data to csv
        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(
                    self.output_dir, "Organism_hierarchy.csv"
                )
            else:
                full_path = os.path.join(os.getcwd(), "Organism_hierarchy.csv")

            df_list = [
                {"child_id": child, "parent_id": parent, "label": label}
                for _, child, parent, label, _ in edge_list
            ]
            df = pd.DataFrame.from_records(df_list)
            df.to_csv(full_path, index=False)
            logger.info(f"Organism hierarchy edge data is written: {full_path}")

        return edge_list

    @staticmethod
    def _clean(value: str) -> str:
        return value.replace("|", ",").replace("'", "^")

    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to database id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)

        return identifier

    def set_node_fields(self, organism_node_fields):
        if organism_node_fields is not None:
            self.organism_node_fields = [
                field.value for field in organism_node_fields
            ]
        else:
            self.organism_node_fields = [
                field.value for field in OrganismNodeField
            ]

    def set_edge_types(self, edge_types):
        self.edge_types = (
            edge_types if edge_types is not None else list(OrganismEdgeType)
        )
