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

On top of modular design of the framework, we devised adapters object-oriented manner. We showed how to create nodes and interaction types from adapters with enum classes, what information provided in the data sources can be used as node/edge attributes, and paved the way for using them in a configurable way. 

To achieve this, we use **enum classes** to standardize the creation of node and edge types, properties, and labels. With these enums, each adapter can be flexibly configured to determine exactly what data it retrieves and how it is represented. This design offers a structured, modular, and consistent approach to managing different data sources within the KG.

---

### Enum Classes in Adapters

Enums ending in `..NodeField` control which properties are created for a specific node type. They map the fields from the source data to the property names in the KG.

Example:
```
class SideEffectNodeField(Enum, metaclass=SideEffectEnumMeta):
    NAME = "name"
    SYNONYMS = "synonyms"
```

By selecting keys from this enum class (explained further below), you can determine which properties the side effect nodes will have in the KG.

Depending on the adapter you can encounter following enum classes:

Enums ending in `..NodeField` control which properties are created for a specific node type. They map the fields from the source data to the property names in the KG.

In many cases, you can also change how properties are represented by editing the values of the enum class:
```
# Default representation
SYNONYMS = "synonyms"

# Custom representation  
SYNONYMS = "alternative names"
```

>:warning:**Important Exception:** While most adapters support this flexibility, the `uniprot_adapter` is a notable exception. You can use `uniprot_adapter.gene_property_name_mappings` or `uniprot_adapter.protein_property_name_mappings` dictionaries for this purpose.



## Installation
The project uses [Poetry](https://python-poetry.org). You can install like this:

```
git clone https://github.com/HUBioDataLab/CROssBAR-BioCypher-Migration.git
cd CROssBAR-BioCypher-Migration
poetry install
```

Poetry will create a virtual environment according to your configuration (either centrally or in the project folder). You can activate it by running `poetry shell` inside the project directory.

### Note about pycurl
You may encounter an error when executing the UniProt adapter about the SSL
backend in pycurl: `ImportError: pycurl: libcurl link-time ssl backend (openssl)
is different from compile-time ssl backend (none/other)`

Should this happen, it can be fixed as described here:
https://stackoverflow.com/questions/68167426/how-to-install-a-package-with-poetry-that-requires-cli-args
by running `poetry shell` followed by `pip list`, noting the version of pycurl,
and then running `pip install --compile --install-option="--with-openssl"
--upgrade --force-reinstall pycurl==<version>` to provide the correct SSL
backend.
