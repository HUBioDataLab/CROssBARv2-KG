
- [CROssBAR-KG](#crossbar-kg)
  - [About the KG](#about-the-kg)
  - [Adapter Usage](#adapter-usage)
    - [Enum Classes in Adapters](#enum-classes-in-adapters)
  - [Installation](#installation)
  - [Data](#data)

# CROssBAR-KG
This is the repo for CROssBARv2 Knowledge Graph (KG) data. CROssBARv2 is a heterogeneous general-purpose biomedical KG-based system. 

This repository provides a collection of configurable and reusable **adapter scripts** designed to collect, process, and harmonize biological data. Each adapter is responsible for handling a specific entity type or data source, making the system flexible and easy to extend.  

You can explore the available adapters in the [bccb directory](/bccb).

These adapters employ the [pypath](https://github.com/saezlab/pypath) for data retrieval and [BioCypher](https://github.com/biocypher/biocypher) for KG creation.


## About the KG

The CROssBARv2 KG comprises approximately 2.7 million nodes spanning 14 distinct node types and around 12.6 million edges representing 51 different edge types, all integrated from 34 biological data sources. We also incorporated several ontologies (e.g., Gene Ontology, Mondo Disease Ontology) along with rich metadata captured as node and edge properties.

Building upon this foundation, we further enhanced the semantic depth of CROssBARv2. This was achieved by generating and storing embeddings for key biological entities, such as proteins, drugs, and Gene Ontology terms. These embeddings are managed using the native [vector index](https://neo4j.com/developer/genai-ecosystem/vector-search/) feature in Neo4j, enabling powerful semantic similarity searches.

![CROssBARv2 KG](https://crossbarv2.hubiodatalab.com/static/images/crossbar-schema.png)

## Adapter Usage

Building on the modular design of the framework, we devised adapters object-oriented manner. We showed how to create nodes and interaction types from adapters with enum classes, what information provided in the data sources can be used as node/edge attributes, and paved the way for highly configurable usage. 

To achieve this, we use **enum classes** to standardize the creation of node and edge types, properties, and labels. With these enums, each adapter can be flexibly configured to determine exactly what data it retrieves and how it is represented. This design offers a structured, modular, and consistent approach to managing different data sources within the KG.

---
### Adapter Workflow

You can import any adapter using the following pattern: 
```
from bccb.adapter_name import ADAPTER_CLASS
```

Example:
```
from bccb.uniprot_adapter import Uniprot
```

All adapters follow a standardized lifecycle consisting of four main steps: Initialization, Data Retrieval, Processing, and Extraction.
```
from bccb.adapter_name import ADAPTER_CLASS

# 1. Initialization: Configure the adapter with necessary arguments
adapter = ADAPTER_CLASS(...)

# 2. Data Retrieval: Download raw data from external sources
adapter.download_data()

# 3. Processing: Harmonize, and filter the data
adapter.process_data()

# 4. Extraction: Retrieve nodes and edges formatted for BioCypher
nodes = adapter.get_nodes()
edges = adapter.get_edges()
```

Adapters are designed to integrate seamlessly with BioCypher. You can configure exactly which node properties, edge properties, and edge types the adapter outputs by passing Enum classes during initialization.

Below is a detailed example using the `Drug` adapter. This demonstrates how to inspect available fields, configure specific outputs, and write the results to Neo4j-importable files using BioCypher.
```
from bccb.drug_adapter import (
    Drug,
    DrugNodeField,
    DrugDTIEdgeField,
    DrugEdgeType
)

from biocypher import BioCypher

# Initialize BioCypher instance
bc = BioCypher(
    biocypher_config_path=...,
    schema_config_path=...
)

print([edge.name for edge in DrugEdgeType])
# Output: ['DRUG_DRUG_INTERACTION', 'DRUG_TARGET_INTERACTION', 'DRUG_GENE_INTERACTION']

print([field.name for field in DrugNodeField])
# Output: ['SMILES', 'INCHI', 'INCHIKEY', 'CAS', 'NAME', 'GROUPS', 'GENERAL_REFERENCES', 
# 'ATC_CODES', 'ZINC', 'CHEMBL', 'BINDINGDB', 'CLINICALTRIALS', 'CHEBI', 'PUBCHEM', 'KEGG_DRUG', 
# 'RXCUI', 'PHARMGKB', 'PDB', 'DRUGCENTRAL', 'SELFORMER_EMBEDDING']

print([field.name for field in DrugDTIEdgeField])
# Output: ['SOURCE', 'MECHANISM_OF_ACTION_TYPE', 'MECHANISM_OF_ACTION', 'REFERENCES', 'KNOWN_ACTION', 
# 'DGIDB_SCORE', 'ACTIVITY_VALUE', 'ACTIVITY_TYPE', 'PCHEMBL', 'CONFIDENCE_SCORE', 'DISEASE_EFFICACY', 
# 'DIRECT_INTERACTION', 'STITCH_COMBINED_SCORE']

# define node fields
node_fields = [DrugNodeField.NAME]

# define dti edge fields
dti_edge_fields = [DrugDTIEdgeField.SOURCE, DrugDTIEdgeField.PCHEMBL, DrugDTIEdgeField.REFERENCES]

# define edge types
edge_types = [DrugEdgeType.DRUG_TARGET_INTERACTION, DrugEdgeType.DRUG_DRUG_INTERACTION]

# initialize Drug adapter
drug_adapter = Drug(
    drugbank_user=...,
    drugbank_password=...,
    edge_types=edge_types,
    node_fields=node_fields,
    dti_edge_fields=dti_edge_fields,
    export_csv=...,
    output_dir=...,
    test_mode=...
)

# download drug data
drug_adapter.download_drug_data(cache=...)

# process drug data
drug_adapter.process_drug_data()

# write Neo4j-importable CSV files via BioCypher
bc.write_nodes(drug_adapter.get_drug_nodes())
bc.write_edges(drug_adapter.get_edges())
bc.write_import_call()
```

Some adapters also provide additional enums for configuration (for example, you can select organisms in the [orthology adapter](/bccb/orthology_adapter.py)). If you do not pass any enum values, the adapter will use all available fields in the relevant enum class by default.

Certain adapters require extra information for data retrieval, such as `drugbank_user` and `drugbank_password` for the `Drug` adapter. Besides enum-related arguments, most adapters share a common set of options:

- `export_csv`: whether to export the processed data as structured CSV files (separate from Neo4j importable CSVs).
- `output_dir`: directory where these exported CSV files will be written.
- `test_mode`: limits the amount of output data for testing or debugging.
- `add_prefix`: whether to add a prefix to node identifiers used in the adapter (e.g., adding `uniprot:` to protein IDs).






### Enum Classes in Adapters

Depending on the adapter, you can encounter following enum classes:

Enums ending in `..NodeField` control which properties are created for a specific node type. They map the fields from the source data to the node property names in the KG.

Example:
```
class SideEffectNodeField(Enum):
    NAME = "name"
    SYNONYMS = "synonyms"
```

By selecting keys from this enum class (explained further below), you can determine which properties the side effect nodes will have in the KG.

In many cases, you can also change how properties are represented by editing the values of the enum class:
```
# Default representation
SYNONYMS = "synonyms"

# Custom representation  
SYNONYMS = "alternative names"
```

>:warning:**Important Exception:** While most adapters support this flexibility, the `uniprot_adapter` is a notable exception. You can use `uniprot_adapter.gene_property_name_mappings` or `uniprot_adapter.protein_property_name_mappings` dictionaries for this purpose.

Similarly, Enums ending in `..EdgeField` control which properties are created for a specific edge type. They map the fields from the source data to the edge property names in the KG.

Example:
```
class TFGenEdgeField(Enum):
    SOURCE = "source"
    PUBMED_ID = "pubmed_id"
    TF_EFFECT = "tf_effect"
```

By selecting keys from this enum class, you can determine which properties the gene regulation edges will have in the KG.

Just like with node fields, you can customize how these edge properties represented in the KG by modifying the enum value.

Enums ending in `..EdgeType` control which edges will be generated as output in particular adapter.


## Installation
The project uses [uv](https://docs.astral.sh/uv/). You can install like this:

Create a virtual environment:
```
uv venv crossbarv2 --python=3.10.8
```
Activate the environment:
```
source crossbarv2/bin/activate
```
Install dependencies:
```
uv pip install -r requirements.txt
```

## Data

Neo4j-importable CSV files required to reconstruct the KG are publicly available at [here](https://drive.google.com/file/d/1KoMAxlvy_4IOo8MPi4TrSbMlQtBf8Pch/view?usp=sharing).
