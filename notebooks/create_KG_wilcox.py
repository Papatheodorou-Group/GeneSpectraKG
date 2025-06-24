from biocypher import BioCypher
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger
import pandas as pd
import hashlib

import os

os.chdir("/Users/ysong/SOFTWARE/GeneSpectraKG")

from genespectrakg.adapters.genespectra_adapter_wilcox_logit import *

bc = BioCypher(
    biocypher_config_path="config/biocypher_config.yaml",
)
# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_types = [
    GeneSpectraAdapterNodeType.CELL_TYPE,
    GeneSpectraAdapterNodeType.GENE,
    GeneSpectraAdapterNodeType.ORTHOLOGOUS_GROUP,
]

# Choose protein adapter fields to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_fields = [
    # Proteins
    GeneSpectraAdapterCellTypeField.SPECIES_OF_ORIGIN,
    GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME,
    GeneSpectraAdapterCellTypeField.CELL_TYPE_ONTOLOGY_NAME,
    GeneSpectraAdapterCellTypeField.CELL_TYPE_ID,
    GeneSpectraAdapterCellTypeField.BROAD_TYPE,
    GeneSpectraAdapterCellTypeField.BROAD_TYPE_2,
    GeneSpectraAdapterCellTypeField.BROAD_TYPE_3,
    GeneSpectraAdapterCellTypeField.BROAD_TAXO_CS,
    GeneSpectraAdapterGeneField.PEPTIDE_ID,
    GeneSpectraAdapterGeneField.PREFERRED_NAME_AND_SPECIES,
    GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX,
    GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN,
    GeneSpectraAdapterGeneField.IS_A_GO_TF,
    GeneSpectraAdapterGeneField.DESCRIPTION,
    GeneSpectraAdapterGeneField.PREFERRED_NAME,
    GeneSpectraAdapterGeneField.PFAMS,
    GeneSpectraAdapterGeneField.GOS,
    GeneSpectraAdapterGeneField.KEGG_KO,
    GeneSpectraAdapterGeneField.KEGG_PATHWAY,
    GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_NAME,
    GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID,
    GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_ID,
]

edge_types = [
    GeneSpectraAdapterEdgeType.GENE_WILCOX_MARKER_IN_CELL_TYPE,
    GeneSpectraAdapterEdgeType.GENE_IN_ORTHOLOGOUS_GROUP,
]

edge_fields = [
    GeneSpectraAdapterEdgeField.AVG_LOG2FC,
    GeneSpectraAdapterEdgeField.P_VAL,
    GeneSpectraAdapterEdgeField.P_VAL_ADJ,
    GeneSpectraAdapterEdgeField.P_VAL_ADJ_RANKING,
    GeneSpectraAdapterEdgeField.AVG_LOG2FC_RANKING,
]


# Create a protein adapter instance
adapter = GeneSpectraAdapter(
    node_types=node_types,
    node_fields=node_fields,
    edge_types=edge_types,
    edge_fields=edge_fields,
)

adapter.load_genespectra_data(
    eggnog_file="data/wilcox/emapper_mammalia_all_species.csv",
    cell_ontology_file="data/wilcox/MTG_cell_type_to_ontology_broad_with_taxo_cs_wilcox.csv",
    wilcox_marker_file="data/wilcox/all_species_wilcox_marker_processed_ensembl_name_use.csv",
)


# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())

bc.write_edges(adapter.get_edges())
# Write admin import statement
bc.write_import_call()
bc.summary()
