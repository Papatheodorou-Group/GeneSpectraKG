from biocypher import BioCypher
from template_package.adapters.genespectra_adapter import (
    GeneSpectraAdapter,
    GeneSpectraAdapterNodeType,
    GeneSpectraAdapterEdgeType,
    GeneSpectraAdapterCellTypeField,
    GeneSpectraAdapterGeneField,
    GeneSpectraAdapterOrthologousGroupField,
    GeneSpectraAdapterSpeciesField,
)

# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher(
    biocypher_config_path="config/biocypher_config.yaml",
)

# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_types = [
    GeneSpectraAdapterNodeType.CELL_TYPE,
    GeneSpectraAdapterNodeType.GENE,
    GeneSpectraAdapterNodeType.SPECIES,
    GeneSpectraAdapterNodeType.ORTHOLOGOUS_GROUP,
]

# Choose protein adapter fields to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_fields = [
    # Proteins
    GeneSpectraAdapterCellTypeField.CELL_TYPE_ID,
    GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME,
    GeneSpectraAdapterCellTypeField.TISSUE_ID,
    GeneSpectraAdapterCellTypeField.TISSUE_NAME,
    GeneSpectraAdapterGeneField.GENE_ID,
    GeneSpectraAdapterGeneField.GENE_NAME,
    GeneSpectraAdapterSpeciesField.SPECIES_ID,
    GeneSpectraAdapterSpeciesField.SPECIES_NAME,
    GeneSpectraAdapterOrthologousGroupField.DATASET,
    GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID,
]

edge_types = [
    GeneSpectraAdapterEdgeType.CELL_TYPE_FROM_SPECIES,
    GeneSpectraAdapterEdgeType.GENE_ENHANCED_IN_CELL_TYPE,
    GeneSpectraAdapterEdgeType.GENE_FROM_SPECIES,
    GeneSpectraAdapterEdgeType.GENE_IN_ORTHOLOGOUS_GROUP,
]


# Create a protein adapter instance
adapter = GeneSpectraAdapter(
    node_types=node_types,
    node_fields=node_fields,
    edge_types=edge_types,
    # we can leave edge fields empty, defaulting to all fields in the adapter
)

adapter.load_genespectra_data()


# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Write admin import statement
bc.write_import_call()

# Print summary
bc.summary()