import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
# Add parent directory to Python path
sys.path.insert(0, str(project_root))

from bccb.uniprot_adapter import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)

from bccb.organism_adapter import (
    Organism
)

from bccb.ppi_adapter import (
    PPI
)

from bccb.interpro_adapter import (
    InterPro
)

from bccb.go_adapter import (
    GO
)

from bccb.drug_adapter import (
    Drug
)

from bccb.compound_adapter import (
    Compound
)

from bccb.orthology_adapter import (
    Orthology
)

from bccb.disease_adapter import (
    Disease
)

from bccb.phenotype_adapter import (
    HPO
)

from bccb.pathway_adapter import (
    Pathway
)

from bccb.side_effect_adapter import (
    SideEffect
)

from bccb.ec_adapter import (
    EC
)

from bccb.tfgen_adapter import (
    TFGene
)

from biocypher import BioCypher

bc = BioCypher(biocypher_config_path= project_root / "config/biocypher_config.yaml",
               schema_config_path= project_root / "config/schema_config.yaml",
)

# Whether to cache data by pypath for future usage
CACHE = True

# Flag for exporting node and edge files as csv format
export_as_csv = True

# Flag for test mode
TEST_MODE = True

# dirs
output_dir_path = "YOUR_PATH"

# user and passwd
drugbank_user = "YOUR_DRUGBANK_USER"
drugbank_passwd = "YOUR_DRUGBANK_PASSWD"

# embedding files, download them from the link in the README and unzip into the project root
embeddings_dir = project_root / "embeddings"

# uniprot configuration
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    UniprotNodeType.GENE,
]

uniprot_node_fields = [
    UniprotNodeField.PRIMARY_GENE_NAME,
    UniprotNodeField.LENGTH,
    UniprotNodeField.MASS,
    UniprotNodeField.ORGANISM,
    UniprotNodeField.ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.ENSEMBL_GENE_IDS,
    UniprotNodeField.ENTREZ_GENE_IDS,
    UniprotNodeField.KEGG_IDS,
    UniprotNodeField.PROTEOME,
    UniprotNodeField.SEQUENCE,
    UniprotNodeField.PROTT5_EMBEDDING,
    UniprotNodeField.ESM2_EMBEDDING,
    UniprotNodeField.NT_EMBEDDING,
]

uniprot_edge_types = [
     UniprotEdgeType.PROTEIN_TO_ORGANISM,
     UniprotEdgeType.GENE_TO_PROTEIN,
]

uniprot_id_type = [
     UniprotIDField.GENE_ENTREZ_ID,
]


uniprot_adapter = Uniprot(
        organism="*",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        id_fields=uniprot_id_type,
        test_mode=TEST_MODE,
    )

uniprot_adapter.download_uniprot_data(
    cache=CACHE,
    retries=6,
    prott5_embedding_output_path=embeddings_dir / "prott5_protein_embeddings.h5",
    esm2_embedding_path=embeddings_dir / "esm2_t33_650M_UR50D_protein_embedding.h5",
    nucleotide_transformer_embedding_path=embeddings_dir / "nucleotide_transformerv2_2.5b_multispecies_embeddings.h5",
)

uniprot_edges = uniprot_adapter.get_edges()


bc.write_nodes(uniprot_adapter.get_nodes())
bc.write_edges(uniprot_edges)


if export_as_csv:
    uniprot_adapter.export_data_to_csv(path=output_dir_path,
                                    node_data=uniprot_adapter.get_nodes(),
                                    edge_data=uniprot_edges)

# PPI
ppi_adapter = PPI(organism=None, 
                  output_dir=output_dir_path,
                  export_csv=export_as_csv,
                  test_mode=TEST_MODE)

ppi_adapter.download_ppi_data(cache=CACHE)

ppi_adapter.process_ppi_data()

bc.write_edges(ppi_adapter.get_ppi_edges())

# protein domain
interpro_adapter = InterPro(
    organism="*", 
    test_mode=TEST_MODE
)

interpro_adapter.download_interpro_data(
    cache=CACHE,
    dom2vec_embedding_path=embeddings_dir / "dom2vec_domain_embedding.h5",
)


if export_as_csv:
    interpro_adapter.export_as_csv(path=output_dir_path)

bc.write_nodes(interpro_adapter.get_interpro_nodes())
bc.write_edges(interpro_adapter.get_interpro_edges())

# gene ontology

go_adapter = GO(
    organism="*", 
    test_mode=TEST_MODE
)
go_adapter.download_go_data(
    cache=CACHE,
    anc2vec_embedding_path=embeddings_dir / "anc2vec_go_term_embedding.h5",
)
bc.write_nodes(go_adapter.get_go_nodes())
bc.write_edges(go_adapter.get_go_edges())
if export_as_csv:
    go_adapter.export_as_csv(path=output_dir_path)

# drug

drug_adapter = Drug(
    drugbank_user=drugbank_user, 
    drugbank_passwd=drugbank_passwd,    
    export_csv=export_as_csv, 
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)
drug_adapter.download_drug_data(
    cache=CACHE,
    selformer_embedding_path=embeddings_dir / "selformer_drug_embeddings.h5",
)
drug_adapter.process_drug_data()
bc.write_nodes(drug_adapter.get_drug_nodes())
bc.write_edges(drug_adapter.get_edges())

# compound
compound_adapter = Compound(
    stitch_organism="*",
    export_csv=export_as_csv, 
    output_dir=output_dir_path
)
compound_adapter.download_compound_data(
    cache=CACHE,
    selformer_embedding_path=embeddings_dir / "selformer_compound_embeddings.h5",
)
compound_adapter.process_compound_data()
bc.write_nodes(compound_adapter.get_compound_nodes())
bc.write_edges(compound_adapter.get_cti_edges())

# orthology
orthology_adapter = Orthology(
    export_csv=export_as_csv, 
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)
orthology_adapter.download_orthology_data(cache=CACHE)
bc.write_edges(orthology_adapter.get_orthology_edges())

# disease
disease_adapter = Disease(
    drugbank_user=drugbank_user, 
    drugbank_passwd=drugbank_passwd,
    export_csv=export_as_csv, 
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)
disease_adapter.download_disease_data(
    cache=CACHE,
    doc2vec_embedding_path=embeddings_dir / "doc2vec_disease_embedding.h5",
)
bc.write_nodes(disease_adapter.get_nodes())
disease_edges = disease_adapter.get_edges()
bc.write_edges(disease_edges)

# organism
# only the organisms the uniprot and disease edges link to, and their
# ancestors up to the root
organism_tax_ids = {
    node_id
    for edges in (uniprot_edges, disease_edges)
    for _, source_id, target_id, _, _ in edges
    for node_id in (source_id, target_id)
    if isinstance(node_id, str) and node_id.startswith("ncbitaxon:")
}
organism_adapter = Organism(
    tax_ids=organism_tax_ids,
    export_csv=export_as_csv,
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)
organism_adapter.download_organism_data(cache=CACHE)
bc.write_nodes(organism_adapter.get_nodes())
bc.write_edges(organism_adapter.get_edges())

# phenotype
phenotype_adapter = HPO(
    export_csv=export_as_csv, 
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)
phenotype_adapter.download_hpo_data(
    cache=CACHE,
    cada_embedding_path=embeddings_dir / "cada_phenotype_embedding.h5",
)
bc.write_nodes(phenotype_adapter.get_nodes())
bc.write_edges(phenotype_adapter.get_edges())

# pathway
pathway_adapter = Pathway(
    drugbank_user=drugbank_user, 
    drugbank_passwd=drugbank_passwd,
    export_csv=export_as_csv,
    output_dir=output_dir_path,
    test_mode=TEST_MODE,
)
pathway_adapter.download_pathway_data(
    cache=CACHE,
    biokeen_embedding_path=embeddings_dir / "biokeen_pathway_embedding.h5",
)
bc.write_nodes(pathway_adapter.get_nodes())
bc.write_edges(pathway_adapter.get_edges())


# side effect
side_effect_adapter = SideEffect(
    drugbank_user=drugbank_user,
    drugbank_passwd=drugbank_passwd,
    export_csv=export_as_csv,
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)

side_effect_adapter.download_side_effect_data(cache=CACHE)
bc.write_nodes(side_effect_adapter.get_nodes())
bc.write_edges(side_effect_adapter.get_edges())

# ec numbers
ec_adapter = EC(
    export_csv=export_as_csv,
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)
ec_adapter.download_ec_data(
    cache=CACHE,
    rxnfp_embedding_path=embeddings_dir / "rxnfp_ec_number_embedding.h5",
)
bc.write_nodes(ec_adapter.get_nodes())
bc.write_edges(ec_adapter.get_edges())

# tf-gen
tfgene_adapter = TFGene(
    export_csv=export_as_csv,
    output_dir=output_dir_path,
    test_mode=TEST_MODE
)
tfgene_adapter.download_tfgen_data(cache=CACHE)
bc.write_edges(tfgene_adapter.get_edges())


# Write import call and other post-processing
bc.write_import_call()
bc.summary()
