#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Melih Darcan
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import os
import collections
import json
from getpass import getpass

from typing import List, Union, NamedTuple, Dict, Tuple, Optional, Literal

import pypath.share.curl as curl
import pypath.share.session as session
import pypath.resources.urls as urls

_logger = session.Logger(name="disgenet_input")
_log = _logger._log


class DisgenetApi:

    _name = "DisGeNET"
    _api_url = urls.urls["disgenet"]["api_url"]
    _authenticated: bool = False
    _api_key: str = None

    def authenticate(self) -> bool:
        """
        Starts an authorization process in DisGeNET API.
        Returns a boolean which is success of authentication.
        """

        if self._authenticated and self._api_key != None:
            return True

        _log(f"Authorizing in {self._name} API...")
        # DisGeNET auth is a static per-user API key, no token exchange needed.
        # Checks DISGENET_API_KEY env var first so long batch runs aren't
        # interrupted by a prompt; falls back to getpass otherwise.
        api_key: str = os.environ.get("DISGENET_API_KEY") or getpass("API Key: ")

        if not api_key:
            _log(f"No API key provided for {self._name} API.")
            self._authenticated = False
            self._api_key = None
            return False

        self._api_key = api_key
        self._authenticated = True

        return self._authenticated and (self._api_key != None)

    def _if_authenticated(f):
        """
        Simple wrapper to get rid of the burden of
        checking authentication status and authenticating
        if already haven't.

        Wraps _retrieve_data, whose callers unpack a (payload, paging)
        tuple, so a failed authentication has to return that shape too.
        """

        def wrapper(self, *args, **kwargs):
            if self.authenticate():
                return f(self, *args, **kwargs)

            else:
                _log("DisGeNET failure in authorization, check your credentials.")
                return None, None

        return wrapper

    def _delete_cache(f):
        """
        A necessary wrapper as databases may be updated
        after the initial download of a particular data.
        """

        def wrapper(*args, **kwargs):
            with curl.cache_delete_on():
                return f(*args, **kwargs)

        return wrapper


    def get_vda_summary(
    self,
    variant: Union[str, List[str]] = None,
    gene_ncbi_id: Union[str, List[str]] = None,
    gene_symbol: Union[str, List[str]] = None,
    disease: Union[str, List[str]] = None,
    source: Union[str, List[str]] = None,
    min_score: float = None,
    max_score: float = None,
    min_ei: float = None,
    max_ei: float = None,
    min_dsi: float = None,
    max_dsi: float = None,
    min_dpi: float = None,
    max_dpi: float = None,
    dis_class_list: Union[str, List[str]] = None,
    page_number: int = None,
) -> NamedTuple(
    "VariantDiseaseAssociation",
    [
        ("variantid", str),
        ("gene_symbol", Tuple[str]),
        ("variant_dsi", float),
        ("variant_dpi", float),
        ("variant_consequence_type", str),
        ("diseaseid", str),
        ("disease_name", str),
        ("disease_classes_do", Tuple[str]),
        ("disease_classes_hpo", Tuple[str]),
        ("disease_classes_msh", Tuple[str]),
        ("disease_classes_umls_st", Tuple[str]),
        ("disease_type", str),
        ("score", float),
        ("ei", float),
        ("year_initial", int),
        ("year_final", int),
        ("source", Tuple[str]),
    ],
):
        """
        Returns Variant-Disease Associations.

        @diseaseClasses_HPO (returned field): Populated only for some
        variant/gene combinations; None otherwise.

        @source (returned field): Source names derived from
        scoreBreakdown, collected across all of its components and
        deduplicated; may contain multiple values.


        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants, up to 100.
        @gene_ncbi_id / @gene_symbol: Union[str, List[str]]
        Gene identifier, as NCBI ID or symbol (separate params).
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the VDA.
        @min_score/@max_score, @min_ei/@max_ei, @min_dsi/@max_dsi,
        @min_dpi/@max_dpi: float
            Score ranges, in [0,1]. min_score/max_score filter on
            normalized_score, matching the returned score field.
        @dis_class_list: Union[str, List[str]]
            MeSH Disease Classes.
        @page_number: int
            Page number (100 results per page; TRIAL accounts are capped
            at the top-30 results and do not support pagination).
        """

        url = f"{self._api_url}/vda/summary"
        get_params = dict()

        if variant != None:
            get_params["variant"] = self._list_to_str(variant, "Variant ID", limit=100)

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if min_ei != None:
            get_params["min_ei"] = str(min_ei)

        if max_ei != None:
            get_params["max_ei"] = str(max_ei)

        if min_dsi != None:
            get_params["min_dsi"] = str(min_dsi)

        if max_dsi != None:
            get_params["max_dsi"] = str(max_dsi)

        if min_dpi != None:
            get_params["min_dpi"] = str(min_dpi)

        if max_dpi != None:
            get_params["max_dpi"] = str(max_dpi)

        if dis_class_list != None:
            get_params["dis_class_list"] = dis_class_list

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result, paging = self._retrieve_data(url, get_params)

        if result == None:
            return None, None

        VariantDiseaseAssociation = collections.namedtuple(
            "VariantDiseaseAssociation",
            [
                "variantid",
                "gene_symbol",
                "variant_dsi",
                "variant_dpi",
                "variant_consequence_type",
                "diseaseid",
                "disease_name",
                "disease_classes_do",
                "disease_classes_hpo",
                "disease_classes_msh",
                "disease_classes_umls_st",
                "disease_type",
                "score",
                "ei",
                "year_initial",
                "year_final",
                "source",
            ],
        )

        for index, entry in enumerate(result):
            result[index] = VariantDiseaseAssociation(
                self._get_string(entry.get("variantStrID")),
                tuple(entry["geneSymbol_keyword"]) if entry.get("geneSymbol_keyword") else None,
                self._get_float(entry.get("variantDSI")),
                self._get_float(entry.get("variantDPI")),
                self._get_string(entry.get("mostSevereConsequences")),
                self._get_string(entry.get("diseaseUMLSCUI")),
                self._get_string(entry.get("diseaseName")),
                tuple(entry["diseaseClasses_DO"]) if entry.get("diseaseClasses_DO") else None,
                tuple(entry["diseaseClasses_HPO"]) if entry.get("diseaseClasses_HPO") else None,
                tuple(entry["diseaseClasses_MSH"]) if entry.get("diseaseClasses_MSH") else None,
                tuple(entry["diseaseClasses_UMLS_ST"]) if entry.get("diseaseClasses_UMLS_ST") else None,
                self._get_string(entry.get("diseaseType")).strip("[]") if entry.get("diseaseType") else None,
                self._get_float(entry.get("normalized_score")),
                self._get_float(entry.get("ei")),
                self._get_int(entry.get("yearInitial")),
                self._get_int(entry.get("yearFinal")),
                self._get_sources(entry),
            )

        return result, paging

    def get_gda_summary(
        self,
        gene_ncbi_id: Union[str, List[str]] = None,
        gene_ensembl_id: Union[str, List[str]] = None,
        gene_symbol: Union[str, List[str]] = None,
        uniprot_id: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        source: Union[str, List[str]] = None,
        min_score: float = None,
        max_score: float = None,
        min_ei: float = None,
        max_ei: float = None,
        min_dsi: float = None,
        max_dsi: float = None,
        min_dpi: float = None,
        max_dpi: float = None,
        min_pli: float = None,
        max_pli: float = None,
        type: str = None,
        dis_class_list: Union[str, List[str]] = None,
        page_number: int = None,
    ) -> NamedTuple(
        "GeneDiseaseAssociation",
        [
            ("geneid", int),
            ("gene_symbol", str),
            ("uniprotid", Tuple[str]),
            ("gene_dsi", float),
            ("gene_dpi", float),
            ("gene_pli", float),
            ("protein_class", Tuple[str]),
            ("protein_class_name", Tuple[str]),
            ("diseaseid", str),
            ("disease_name", str),
            ("disease_classes_do", Tuple[str]),
            ("disease_classes_hpo", Tuple[str]),
            ("disease_classes_msh", Tuple[str]),
            ("disease_classes_umls_st", Tuple[str]),
            ("disease_type", str),
            ("score", float),
            ("ei", float),
            ("el", str),
            ("year_initial", int),
            ("year_final", int),
            ("source", Tuple[str]),
        ],
    ):
        """
        Returns Gene-Disease Associations.

        @source (returned field): Source names derived from
        scoreBreakdown, collected across all of its components and
        deduplicated; may contain multiple values.

        @gene_ncbi_id / @gene_ensembl_id / @gene_symbol: Union[str, List[str]]
            Gene identifier(s) in the respective vocabulary, up to 100.
        @uniprot_id: Union[str, List[str]]
            Uniprot accession(s), up to 100.
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the GDA.
        @min_score/@max_score, @min_ei/@max_ei, @min_dsi/@max_dsi,
        @min_dpi/@max_dpi, @min_pli/@max_pli: float
            Score ranges, in [0,1]. min_score/max_score filter on
            normalized_score, matching the returned score field.
        @type: str
            DisGeNET Disease Type ("disease", "phenotype", "group").
        @dis_class_list: Union[str, List[str]]
            MeSH Disease Classes.
        @page_number: int
            Page number (100 results per page; TRIAL accounts are capped
            at the top-30 results and do not support pagination).
        """

        url = f"{self._api_url}/gda/summary"
        get_params = dict()

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_ensembl_id != None:
            get_params["gene_ensembl_id"] = self._list_to_str(gene_ensembl_id, "Gene Ensembl ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if uniprot_id != None:
            get_params["uniprot_id"] = self._list_to_str(uniprot_id, "Uniprot ID", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if min_ei != None:
            get_params["min_ei"] = str(min_ei)

        if max_ei != None:
            get_params["max_ei"] = str(max_ei)

        if min_dsi != None:
            get_params["min_dsi"] = str(min_dsi)

        if max_dsi != None:
            get_params["max_dsi"] = str(max_dsi)

        if min_dpi != None:
            get_params["min_dpi"] = str(min_dpi)

        if max_dpi != None:
            get_params["max_dpi"] = str(max_dpi)

        if min_pli != None:
            get_params["min_pli"] = str(min_pli)

        if max_pli != None:
            get_params["max_pli"] = str(max_pli)

        if type != None:
            get_params["type"] = type

        if dis_class_list != None:
            get_params["dis_class_list"] = dis_class_list

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result, paging = self._retrieve_data(url, get_params)

        if result == None:
            return None, None

        GeneDiseaseAssociation = collections.namedtuple(
            "GeneDiseaseAssociation",
            [
                "geneid",
                "gene_symbol",
                "uniprotid",
                "gene_dsi",
                "gene_dpi",
                "gene_pli",
                "protein_class",
                "protein_class_name",
                "diseaseid",
                "disease_name",
                "disease_classes_do",
                "disease_classes_hpo",
                "disease_classes_msh",
                "disease_classes_umls_st",
                "disease_type",
                "score",
                "ei",
                "el",
                "year_initial",
                "year_final",
                "source",
            ],
        )

        for index, entry in enumerate(result):
            result[index] = GeneDiseaseAssociation(
                self._get_int(entry.get("geneNcbiID")),
                self._get_string(entry.get("symbolOfGene")),
                tuple(entry["geneProteinStrIDs"]) if entry.get("geneProteinStrIDs") else None,
                self._get_float(entry.get("geneDSI")),
                self._get_float(entry.get("geneDPI")),
                self._get_float(entry.get("genepLI")),
                tuple(entry["geneProteinClassIDs"]) if entry.get("geneProteinClassIDs") else None,
                tuple(entry["geneProteinClassNames"]) if entry.get("geneProteinClassNames") else None,
                self._get_string(entry.get("diseaseUMLSCUI")),
                self._get_string(entry.get("diseaseName")),
                tuple(entry["diseaseClasses_DO"]) if entry.get("diseaseClasses_DO") else None,
                tuple(entry["diseaseClasses_HPO"]) if entry.get("diseaseClasses_HPO") else None,
                tuple(entry["diseaseClasses_MSH"]) if entry.get("diseaseClasses_MSH") else None,
                tuple(entry["diseaseClasses_UMLS_ST"]) if entry.get("diseaseClasses_UMLS_ST") else None,
                self._get_string(entry.get("diseaseType")).strip("[]") if entry.get("diseaseType") else None,
                self._get_float(entry.get("normalized_score")),
                self._get_float(entry.get("ei")),
                self._get_string(entry.get("el")),
                self._get_int(entry.get("yearInitial")),
                self._get_int(entry.get("yearFinal")),
                self._get_sources(entry),
            )

        return result, paging




    def get_vda_evidence(
        self,
        variant: Union[str, List[str]] = None,
        gene_ncbi_id: Union[str, List[str]] = None,
        gene_symbol: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        source: Union[str, List[str]] = None,
        min_pmYear: str = None,
        max_pmYear: str = None,
        min_score: float = None,
        max_score: float = None,
        page_number: int = None,
    ):
        """
        Returns evidences that support Variant-Disease Associations. Does
        not convert results into a namedtuple, it returns the raw list of
        record dicts as-is.

        @variant: Union[str, List[str]]
            Variant (dbSNP Identifier) or list of variants, up to 100.
        @gene_ncbi_id / @gene_symbol: Union[str, List[str]]
            Gene identifier, as NCBI ID or symbol (separate params).
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the VDA.
        @min_pmYear/@max_pmYear: str
            Publication year range (filters on pmYear, not the evidence
            timestamp).
        @min_score/@max_score: float
            Variant-disease normalized score range, in [0,1].
        @page_number: int
            Page number. Returns a single page; call again with an
            incremented page_number for more. TRIAL accounts do not
            support pagination.
        """

        url = f"{self._api_url}/vda/evidence"
        get_params = dict()

        if variant != None:
            get_params["variant"] = self._list_to_str(variant, "Variant ID", limit=100)

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_pmYear != None:
            get_params["min_pmYear"] = min_pmYear

        if max_pmYear != None:
            get_params["max_pmYear"] = max_pmYear

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result, _ = self._retrieve_data(url, get_params)

        return result



    def get_gda_evidence(
        self,
        gene_ncbi_id: Union[str, List[str]] = None,
        gene_symbol: Union[str, List[str]] = None,
        disease: Union[str, List[str]] = None,
        source: Union[str, List[str]] = None,
        min_pmYear: str = None,
        max_pmYear: str = None,
        min_score: float = None,
        max_score: float = None,
        page_number: int = None,
    ):
        """
        Returns evidences that support Gene-Disease Associations. Does not
        convert results into a namedtuple either, it returns the raw list
        of record dicts as-is.

        @gene_ncbi_id / @gene_symbol: Union[str, List[str]]
            Gene identifier, as NCBI ID or symbol (separate params).
        @disease: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @source: Union[str, List[str]]
            Source of the GDA.
        @min_pmYear/@max_pmYear: str
            Publication year range (filters on pmYear, not the evidence
            timestamp).
        @min_score/@max_score: float
            Gene-disease normalized score range, in [0,1].
        @page_number: int
            Page number. Returns a single page; call again with an
            incremented page_number for more. TRIAL accounts do not
            support pagination.
        """

        url = f"{self._api_url}/gda/evidence"
        get_params = dict()

        if gene_ncbi_id != None:
            get_params["gene_ncbi_id"] = self._list_to_str(gene_ncbi_id, "Gene NCBI ID", limit=100)

        if gene_symbol != None:
            get_params["gene_symbol"] = self._list_to_str(gene_symbol, "Gene Symbol", limit=100)

        if disease != None:
            get_params["disease"] = self._list_to_str(disease, "Disease ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_pmYear != None:
            get_params["min_pmYear"] = min_pmYear

        if max_pmYear != None:
            get_params["max_pmYear"] = max_pmYear

        if min_score != None:
            get_params["min_score"] = str(min_score)

        if max_score != None:
            get_params["max_score"] = str(max_score)

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result, _ = self._retrieve_data(url, get_params)

        return result



    def get_dda(
        self,
        disease_1: Union[str, List[str]],
        disease_2: Union[str, List[str]] = None,
        source: Union[str, List[str]] = None,
        min_jaccard_genes: float = None,
        min_jaccard_variants: float = None,
        order_by: str = None,
        order: str = None,
        page_number: int = None,
    ) -> NamedTuple(
        "DiseaseDiseaseAssociation",
        [
            ("disease1_name", str),
            ("disease2_name", str),
            ("disease1_classes_do", Tuple[str]),
            ("disease1_classes_hpo", Tuple[str]),
            ("disease1_classes_msh", Tuple[str]),
            ("disease1_classes_umls_st", Tuple[str]),
            ("disease2_classes_do", Tuple[str]),
            ("disease2_classes_hpo", Tuple[str]),
            ("disease2_classes_msh", Tuple[str]),
            ("disease2_classes_umls_st", Tuple[str]),
            ("jaccard_genes", float),
            ("pvalue_jaccard_genes", float),
            ("jaccard_variants", float),
            ("pvalue_jaccard_variants", float),
            ("source", str),
            ("ngenes1", int),
            ("ngenes2", int),
            ("nvariants1", int),
            ("nvariants2", int),
            ("shared_genes", int),
            ("shared_variants", int),
            ("diseaseid1", str),
            ("diseaseid2", str),
        ],
    ):
        """
        Returns Disease-Disease Associations for disease_1, optionally
        restricted to the diseases given in disease_2.
        The API returns gene-sharing and variant-sharing metrics together
        in a single query, but a given pair only carries the metrics it
        actually has: a pair that shares genes but no variants comes back
        with jaccard_variants (and pvalue_jaccard_variants) set to None.
        This is common rather than exceptional, so callers must guard
        against None before doing arithmetic on either jaccard field.

        @disease_1: Union[str, List[str]]
            Disease id(s) with vocabulary prefix, e.g. "UMLS_C0005745".
        @disease_2: Union[str, List[str]], optional
            Disease id(s) with vocabulary prefix. If omitted, returns all
            diseases sharing genes/variants with disease_1.
        @source: Union[str, List[str]]
            Source of the DDA.
        @min_jaccard_genes / @min_jaccard_variants: float
            Minimum jaccard index on shared genes / shared variants,
            in [0,1]. Records below the threshold are excluded
            server-side, so this also caps the number of results
            returned for a given disease_1 without needing pagination.
        @order_by: str
            Field to sort by, e.g. "jaccard_genes" or "jaccard_variants".
        @order: str
            Sort direction, "ASC" or "DESC". Only meaningful together
            with order_by.
        @page_number: int
            Page number (100 results per page; TRIAL accounts are capped
            at the top-10 results and do not support pagination).
        """

        disease_1 = self._list_to_str(disease_1, "Disease 1 ID", limit=100)

        url = f"{self._api_url}/dda"
        get_params = dict()

        get_params["disease_1"] = disease_1

        if disease_2 != None:
            get_params["disease_2"] = self._list_to_str(disease_2, "Disease 2 ID", limit=100)

        if source != None:
            get_params["source"] = source

        if min_jaccard_genes != None:
            get_params["min_jaccard_genes"] = str(min_jaccard_genes)

        if min_jaccard_variants != None:
            get_params["min_jaccard_variants"] = str(min_jaccard_variants)

        if order_by != None:
            get_params["order_by"] = order_by

        if order != None:
            get_params["order"] = order

        if page_number != None:
            get_params["page_number"] = str(page_number)

        result, paging = self._retrieve_data(url, get_params)

        if result == None:
            return None, None

        DiseaseDiseaseAssociation = collections.namedtuple(
            "DiseaseDiseaseAssociation",
            [
                "disease1_name",
                "disease2_name",
                "disease1_classes_do",
                "disease1_classes_hpo",
                "disease1_classes_msh",
                "disease1_classes_umls_st",
                "disease2_classes_do",
                "disease2_classes_hpo",
                "disease2_classes_msh",
                "disease2_classes_umls_st",
                "jaccard_genes",
                "pvalue_jaccard_genes",
                "jaccard_variants",
                "pvalue_jaccard_variants",
                "source",
                "ngenes1",
                "ngenes2",
                "nvariants1",
                "nvariants2",
                "shared_genes",
                "shared_variants",
                "diseaseid1",
                "diseaseid2",
            ],
        )

        for index, entry in enumerate(result):
            result[index] = DiseaseDiseaseAssociation(
                self._get_string(entry.get("disease1_Name")),
                self._get_string(entry.get("disease2_Name")),
                tuple(entry["disease1_Classes_DO"]) if entry.get("disease1_Classes_DO") else None,
                tuple(entry["disease1_Classes_HPO"]) if entry.get("disease1_Classes_HPO") else None,
                tuple(entry["disease1_Classes_MSH"]) if entry.get("disease1_Classes_MSH") else None,
                tuple(entry["disease1_Classes_UMLS_ST"]) if entry.get("disease1_Classes_UMLS_ST") else None,
                tuple(entry["disease2_Classes_DO"]) if entry.get("disease2_Classes_DO") else None,
                tuple(entry["disease2_Classes_HPO"]) if entry.get("disease2_Classes_HPO") else None,
                tuple(entry["disease2_Classes_MSH"]) if entry.get("disease2_Classes_MSH") else None,
                tuple(entry["disease2_Classes_UMLS_ST"]) if entry.get("disease2_Classes_UMLS_ST") else None,

                self._get_float(entry.get("jaccard_genes")),
                self._get_float(entry.get("pvalue_jaccard_genes")),
                self._get_float(entry.get("jaccard_variants")),
                self._get_float(entry.get("pvalue_jaccard_variants")),
                self._get_string(entry.get("source")),
                self._get_int(entry.get("ngenes_diseaseID_1")),
                self._get_int(entry.get("ngenes_diseaseID_2")),
                self._get_int(entry.get("nvariants_diseaseID_1")),
                self._get_int(entry.get("nvariants_diseaseID_2")),
                self._get_int(entry.get("shared_genes")),
                self._get_int(entry.get("shared_variants")),
                self._get_string(entry.get("disease1_UMLSCUI")),
                self._get_string(entry.get("disease2_UMLSCUI")),
            )

        return result, paging



    def _list_to_str(
        self, list_obj: List[str], name: str, limit: Optional[int] = None
    ) -> List[str]:
        """
        Joins the list object and returns a string

        @list_obj : List[str]
            List object to be processed
        @name: str
            Name of items, like Gene ID or something
        @limit: Optional[int]
            Maximum number of list items
        """

        if isinstance(list_obj, list):
            if limit != None and len(list_obj) > limit:
                _log(
                    f"DisGeNET maximum length of {name}'s are {limit}, "
                    f"first {limit} {name}'s will be used."
                )

                return ",".join(list_obj[:limit])

            return ",".join(list_obj)

        return list_obj

    # Guards against calling the API without a valid key.
    @_if_authenticated
    def _retrieve_data(
        self, url: str, get_params: Union[List[str], Dict[str, str]]
    ) -> Tuple[Optional[List[Dict[str, str]]], Optional[Dict[str, str]]]:
        """
        Retrieves the data with given request body

        @url : str
            Query url
        @get_params: Union[List[str], Dict[str, str]]
            GET request parameters

        Returns a (payload, paging) tuple; both are None on failure.
        """

        headers = ["accept: */*", f"Authorization: Bearer {self._api_key}"]

        # List-valued params (e.g. source, dis_class_list) are expanded into
        # repeated "key=value" entries, one per item, since the API doesn't
        # parse a literal Python list string like "source=['CURATED']".
        get_params_list = []
        for key, value in get_params.items():
            if isinstance(value, list):
                get_params_list.extend([f"{key}={item}" for item in value])
            else:
                get_params_list.append(f"{key}={value}")

        get_params = get_params_list

        c = curl.Curl(url=url, get=get_params, req_headers=headers)

        if c.status == 0 or c.status == 200:
            raw_result = c.result

            # The body can contain more than one concatenated JSON document:
            # pypath's retry loop appends a successful retry after the failed
            # attempt (e.g. a quota-exceeded error followed by a valid
            # payload). The valid response is therefore the LAST document,
            # so all documents are parsed and the last one is used.
            decoder = json.JSONDecoder()
            documents = []
            index = 0

            try:
                while index < len(raw_result):
                    # raw_decode() does not skip leading whitespace
                    while index < len(raw_result) and raw_result[index].isspace():
                        index += 1

                    if index >= len(raw_result):
                        break

                    document, index = decoder.raw_decode(raw_result, index)
                    documents.append(document)
            except json.JSONDecodeError as e:
                _log(f"DisGeNET: could not parse response as JSON: {e}")
                _log(f"DisGeNET response body: {repr(raw_result)[:200]}")
                return None, None

            if not documents:
                _log("DisGeNET: empty response body.")
                return None, None

            if len(documents) > 1:
                _log(
                    f"DisGeNET: response contained {len(documents)} JSON "
                    "documents; the last one is used, since a retried "
                    "request is appended after the failed attempt."
                )

            result = documents[-1]

            if not isinstance(result, dict) or result.get("status") != "OK":
                _log(
                    f"DisGeNET: API returned a non-OK status: "
                    f"{repr(raw_result)[:300]}"
                )
                return None, None

            # Records live under "payload"; "paging" carries pageSize,
            # totalElements, and currentPageNumber, which callers doing
            # their own pagination need to detect a short final page,
            # a mismatched page count, or a truncated result reliably —
            # rather than inferring completion from len(result) alone.
            return result.get("payload"), result.get("paging")

        _log(f"DisGeNET: an error occurred with the code {c.status}")
        _log(f"DisGeNET response body: {repr(c.result)[:200]}")
        return None, None # keep return shape consistent for (result, paging) unpacking
    def _get_int(self, str_obj) -> int:
        """
        Returns an int if the object is not None

        @str_obj : str
            String to be processed
        """

        if str_obj != None and not isinstance(str_obj, int):
            return int(str_obj)

        return str_obj

    def _get_float(self, str_obj) -> float:
        """
        Returns a float if the object is not None

        @str_obj : str
            String to be processed
        """

        if str_obj != None and not isinstance(str_obj, float):
            return float(str_obj)

        return str_obj

    def _get_string(self, obj) -> str:
        """
        Returns a string if the object is not None

        @obj : object
            Object to be processed
        """

        if obj != None and not isinstance(obj, str):
            return str(obj).strip()

        return obj

    def _get_sources(self, entry) -> Optional[Tuple[str]]:
        """
        Returns source names derived from scoreBreakdown, if present.

        Collects every component (curated, clinical, inferred, models,
        literature, biobank); callers wanting a single component can
        filter downstream.

        A component may be present but null (biobank commonly is), and
        the same source can appear under more than one component
        (TEXTMINING_MODELS shows up under both models and literature), so
        the union is deduplicated. It is also sorted to keep the value
        stable across runs, since these tuples end up in namedtuples that
        disgenet_annotations() collects into a set.

        @entry : dict
            A single record from a /summary response.
        """
        breakdown = entry.get("scoreBreakdown")
        if not breakdown:
            return None

        components = breakdown[0].get("components") or {}
        sources = set()

        for component in components.values():
            if not component:
                continue

            sources.update(component.get("sources") or ())

        return tuple(sorted(sources)) if sources else None


@DisgenetApi._delete_cache
def variant_gene_mappings(
    api: "DisgenetApi",
    gene_ncbi_ids: List[str],
    batch_size: int = 10,
) -> Tuple[Dict[str, List["VariantGeneMapping"]], List[List[str]]]:
    """
    Builds a {snpId: [VariantGeneMapping(geneId, geneSymbol, sourceIds),
    ...]} mapping by querying the DisGeNET API.

    Returns a (mapping, failed_batches) tuple. A batch whose pagination
    is cut short by a retrieval error - an expired quota mid-run being
    the common case - is appended to failed_batches, so a caller can
    tell a partial result from a complete one instead of silently
    building on truncated input.

    Requires a list of known NCBI Gene IDs to query against - e.g. from
    uniprot_adapter.py's xref_geneid field (same source used for
    disgenet_annotations()). Queries /entity/variant using gene_ncbi_id,
    which returns each matching variant's own variantToGenes field - a
    list of {geneNcbiID, symbolOfGene, geneEnsemblID, sources} grouped
    per variant.

    @api: DisgenetApi
        An already-authenticated DisgenetApi instance.
    @gene_ncbi_ids: List[str]
        Known NCBI Gene (Entrez) IDs to query, e.g. ["672", "675", ...].
    @batch_size: int
        Number of gene IDs to send per request.
    """

    VariantGeneMapping = collections.namedtuple(
        "VariantGeneMapping",
        [
            "geneId",
            "geneSymbol",
            "sourceIds",
        ],
    )

    mapping = dict()
    failed_batches = []

    for i in range(0, len(gene_ncbi_ids), batch_size):
        batch = gene_ncbi_ids[i : i + batch_size]
        page_number = 0

        while True:
            url = f"{api._api_url}/entity/variant"
            get_params = {
                "gene_ncbi_id": batch,
                "page_number": str(page_number),
            }

            result, paging = api._retrieve_data(url, get_params)

            if paging is None:
                # A real error occurred (auth failure, bad response, etc.),
                # not just an empty final page — stop this batch rather than
                # silently treating a failure as "pagination complete".
                _log(
                    f"DisGeNET: stopping pagination for batch {batch} "
                    f"at page {page_number} due to a retrieval error; "
                    f"results for these ids are incomplete."
                )
                failed_batches.append(batch)
                break

            if not result:
                break

            for entry in result:
                snp_id = entry.get("strID")
                variant_to_genes = entry.get("variantToGenes")

                if snp_id == None or not variant_to_genes:
                    continue

                if snp_id not in mapping:
                    mapping[snp_id] = []

                existing_pairs = {(m.geneId, m.geneSymbol) for m in mapping[snp_id]}

                for vtg in variant_to_genes:
                    gene_id = vtg.get("geneNcbiID")
                    gene_symbol = vtg.get("symbolOfGene")
                    sources = vtg.get("sources")

                    gene_id_str = api._get_string(gene_id) if gene_id != None else None
                    pair = (gene_id_str, gene_symbol)

                    if pair in existing_pairs:
                        # variantToGenes is variant-level data: the same
                        # variant reached through a different gene query
                        # comes back with an identical list, sources
                        # included, so the skipped record has nothing to
                        # merge in.
                        continue

                    mapping[snp_id].append(
                        VariantGeneMapping(
                            gene_id_str,
                            gene_symbol,
                            tuple(sources) if sources else None,
                        )
                    )
                    existing_pairs.add(pair)

            if paging.get("totalElementsInPage", 0) < paging.get("pageSize", 100):
                break

            page_number += 1

    return mapping, failed_batches


@DisgenetApi._delete_cache
def disease_id_mappings(
    api: "DisgenetApi",
    mondo_ids: List[str],
    batch_size: int = 10,
) -> Tuple[Dict[str, "DiseaseIdMapping"], List[List[str]]]:
    """
    Builds a {diseaseId: DiseaseIdMapping(name, vocabularies)} mapping by
    querying the DisGeNET API.

    Returns a (mapping, failed_batches) tuple. A batch whose pagination
    is cut short by a retrieval error - an expired quota mid-run being
    the common case - is appended to failed_batches, so a caller can
    tell a partial result from a complete one instead of silently
    building on truncated input.

    Requires a list of known MONDO disease IDs to query against - e.g.
    from disease_adapter.py's MONDO ontology download. Each ID must be in
    the API's expected format, e.g. "MONDO_0007254".

    This is a standalone function (not a DisgenetApi method); it takes an
    already-authenticated DisgenetApi instance as a parameter and calls
    _retrieve_data directly against the /entity/disease endpoint, which
    returns disease-level records (name + diseaseCodes) rather than the
    per-gene-association rows get_gda_summary() returns.

    Note on data shape: a single queried MONDO id can match more than
    one disease entity in the response - e.g. an umbrella disease and
    a related subtype can both list the same MONDO code among their
    diseaseCodes - so the same disease_id can appear more than once
    across a batch's paginated results. We only need each disease's
    name + vocabularies once, so once a disease_id has been recorded
    we skip it on subsequent rows - this does not lose any information,
    since this function only returns name + vocabularies per disease.

    @api: DisgenetApi
        An already-authenticated DisgenetApi instance.
    @mondo_ids: List[str]
        Known MONDO disease IDs to query, e.g. ["MONDO_0007254", ...].
    @batch_size: int
        Number of disease IDs to send per request.
    """
    Vocabulary = collections.namedtuple(
        "Vocabulary",
        [
            "vocabulary",
            "code",
            "vocabularyName",
        ],
    )
    DiseaseIdMapping = collections.namedtuple(
        "DiseaseIdMapping",
        [
            "name",
            "vocabularies",
        ],
    )
    mapping = dict()
    failed_batches = []
    for i in range(0, len(mondo_ids), batch_size):
        batch = mondo_ids[i : i + batch_size]
        page_number = 0
        while True:
            url = f"{api._api_url}/entity/disease"
            get_params = {
                "disease": batch,
                "page_number": str(page_number),
            }
            result, paging = api._retrieve_data(url, get_params)
            if paging is None:
                # A real error occurred (auth failure, bad response, etc.),
                # not just an empty final page — stop this batch rather than
                # silently treating a failure as "pagination complete".
                _log(
                    f"DisGeNET: stopping pagination for batch {batch} "
                    f"at page {page_number} due to a retrieval error; "
                    f"results for these ids are incomplete."
                )
                failed_batches.append(batch)
                break
            if not result:
                break
            for entry in result:
                disease_id = entry.get("diseaseUMLSCUI")
                raw_codes = entry.get("diseaseCodes")
                if disease_id is None or not raw_codes:
                    continue
                # The same disease_id can recur across a batch's results
                # (see docstring note above); only the first occurrence
                # needs to be processed.
                if disease_id in mapping:
                    continue
                name = entry.get("name")
                vocab_list = []
                for raw in raw_codes:
                    vocabulary = raw.get("vocabulary")
                    code = raw.get("code")
                    if not vocabulary or not code:
                        continue
                    vocab_list.append(
                        Vocabulary(
                            vocabulary,
                            code,
                            # diseaseCodes entries carry only the short
                            # vocabulary prefix (e.g. "MSH", "MONDO"),
                            # no full vocabulary name.
                            None,
                        )
                    )
                mapping[disease_id] = DiseaseIdMapping(name, tuple(vocab_list))
            if paging.get("totalElementsInPage", 0) < paging.get("pageSize", 100):
                break
            page_number += 1
    return mapping, failed_batches


@DisgenetApi._delete_cache
def disgenet_annotations(
    api: "DisgenetApi",
    gene_ncbi_ids: List[str],
    dataset: str = "curated",
    batch_size: int = 10,
) -> Tuple[Dict[str, set], List[List[str]]]:
    """
    Builds a {uniprot_id: {DisGeNetAnnotation(...), ...}} mapping by
    querying the DisGeNET API.

    Returns a (mapping, failed_batches) tuple. A batch whose pagination
    is cut short by a retrieval error - an expired quota mid-run being
    the common case - is appended to failed_batches, so a caller can
    tell a partial result from a complete one instead of silently
    building on truncated input.

    Requires a list of known NCBI Gene IDs to query against - e.g. from
    uniprot_adapter.py's xref_geneid field. The API returns UniProt IDs
    directly in geneProteinStrIDs, so those are read straight from the
    response.

    This calls _retrieve_data directly rather than get_gda_summary(),
    because get_gda_summary() does not return numPMIDs or
    numDBSNPsupportingAssociation, which map to this function's
    nof_pmids/nof_snps fields.

    @api: DisgenetApi
        An already-authenticated DisgenetApi instance.
    @gene_ncbi_ids: List[str]
        Known NCBI Gene (Entrez) IDs to query, e.g. ["672", "675", ...].
    @dataset: str
        Only "curated" and "all" are supported (mapped to
        source=["CURATED"] and no source filter, respectively).
    @batch_size: int
        Number of gene IDs to send per request.
    """

    dataset_source_map = {
        "curated": ["CURATED"],
        "all": None,
    }

    if dataset not in dataset_source_map:
        _log(
            f"DisGeNET dataset='{dataset}' is not supported. Only 'curated' and "
            "'all' are mapped to the API's source filter; 'literature' "
            "and 'befree' have no clear equivalent in the source list."
        )
        return {}, []

    source = dataset_source_map[dataset]

    DisGeNetAnnotation = collections.namedtuple(
        "DisGeNetAnnotation",
        [
            "disease",
            "type",
            "score",
            "dsi",
            "dpi",
            "nof_pmids",
            "nof_snps",
            "source",
        ],
    )

    data = collections.defaultdict(set)
    failed_batches = []

    for i in range(0, len(gene_ncbi_ids), batch_size):
        batch = gene_ncbi_ids[i : i + batch_size]
        page_number = 0

        while True:
            url = f"{api._api_url}/gda/summary"
            get_params = {
                "gene_ncbi_id": batch,
                "page_number": str(page_number),
            }

            if source != None:
                get_params["source"] = source

            result, paging = api._retrieve_data(url, get_params)

            if paging is None:
                # A real error occurred (auth failure, bad response, etc.),
                # not just an empty final page — stop this batch rather than
                # silently treating a failure as "pagination complete".
                _log(
                    f"DisGeNET: stopping pagination for batch {batch} "
                    f"at page {page_number} due to a retrieval error; "
                    f"results for these ids are incomplete."
                )
                failed_batches.append(batch)
                break

            if not result:
                break

            for entry in result:
                uniprot_ids = entry.get("geneProteinStrIDs")

                if not uniprot_ids:
                    continue

                disease = entry.get("diseaseName")
                disease_type = entry.get("diseaseType")
                disease_type = disease_type.strip("[]") if disease_type else None
                score = entry.get("normalized_score")
                dsi = entry.get("geneDSI")
                dpi = entry.get("geneDPI")
                nof_pmids = entry.get("numPMIDs")
                nof_snps = entry.get("numDBSNPsupportingAssociation")


                record_source = api._get_sources(entry)

                annotation = DisGeNetAnnotation(
                    disease,
                    disease_type,
                    api._get_float(score) if score != None else None,
                    api._get_float(dsi) if dsi != None else None,
                    api._get_float(dpi) if dpi != None else None,
                    api._get_int(nof_pmids) if nof_pmids != None else None,
                    api._get_int(nof_snps) if nof_snps != None else None,
                    record_source,
                )

                for uniprot in uniprot_ids:
                    data[uniprot].add(annotation)

            if paging.get("totalElementsInPage", 0) < paging.get("pageSize", 100):
                break

            page_number += 1

    return dict(data), failed_batches
