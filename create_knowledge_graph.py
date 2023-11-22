from biocypher import BioCypher
from genespectrakg.adapters.genespectra_adapter_individual import (
    GeneSpectraAdapter,
    GeneSpectraAdapterNodeType,
    GeneSpectraAdapterEdgeType,
    GeneSpectraAdapterCellTypeField,
    GeneSpectraAdapterGeneField,
    GeneSpectraAdapterOrthologousGroupField,
    GeneSpectraAdapterSpeciesField,
    GeneSpectraAdapterEdgeField,
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
    GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_NAME,
    GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID,
    GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_ID,
]

edge_types = [
    GeneSpectraAdapterEdgeType.CELL_TYPE_FROM_SPECIES,
    GeneSpectraAdapterEdgeType.GENE_ENHANCED_IN_CELL_TYPE,
    GeneSpectraAdapterEdgeType.GENE_ENRICHED_IN_CELL_TYPE,
    GeneSpectraAdapterEdgeType.GENE_FROM_SPECIES,
    GeneSpectraAdapterEdgeType.GENE_IN_ORTHOLOGOUS_GROUP,
]

edge_fields = [
    GeneSpectraAdapterEdgeField.SPECIFICITY_CATEGORY,
    GeneSpectraAdapterEdgeField.DISTRIBUTION_CATEGORY,
    GeneSpectraAdapterEdgeField.SPECIFICITY_CATEGORY_TYPE,
    GeneSpectraAdapterEdgeField.FRACTION_EXPRESSED,
    GeneSpectraAdapterEdgeField.MAX_EXPRESSION,
    GeneSpectraAdapterEdgeField.MEAN_EXPRESSION,
    GeneSpectraAdapterEdgeField.NUMBER_EXPRESSED,
    GeneSpectraAdapterEdgeField.SPECIFICITY_SCORE,
]


# Create a protein adapter instance
adapter = GeneSpectraAdapter(
    node_types=node_types,
    node_fields=node_fields,
    edge_types=edge_types,
    edge_fields=edge_fields,
)

adapter.load_genespectra_data(eggnog_file='data/MTG_eggnog_ensembl_mapped_mammalia.csv',
                              cell_ontology_file='data/MTG_cell_type_to_ontology.csv', 
                              genespectra_file='data/human_classes_subclass_processed.csv')


# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Write admin import statement
bc.write_import_call()

# Print summary
bc.summary()