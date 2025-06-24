from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger
import pandas as pd
import hashlib
import numpy as np

logger.debug(f"Loading module {__name__}.")

# reference:
# https://github.com/biocypher/decider-genetics/blob/main/decider_genetics/adapters/cn_genes_adapter.py


class GeneSpectraAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    GENE = auto()
    CELL_TYPE = auto()
    ORTHOLOGOUS_GROUP = auto()


class GeneSpectraAdapterGeneField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    The names correspond to the column names in the genespectra input file
    Some columns are present in the emapper file
    """

    PEPTIDE_ID = "peptide"
    PREFERRED_NAME_AND_SPECIES = "preferred_name_and_species"
    SPECIES_OF_ORIGIN = "species_of_origin"
    IS_A_GO_TF = "is_a_GO_tf"
    DESCRIPTION = "Description"
    PREFERRED_NAME = "Preferred_name"
    PREFERRED_NAME_WILCOX = "ensembl_gene_name_use"  # specifically for the MTG data, the data i have is in ncbi gene names but eggnog is ensembl gene names, so i corrected for differences
    PFAMS = "PFAMs"
    GOS = "GOs"
    KEGG_KO = "KEGG_ko"
    KEGG_PATHWAY = "KEGG_Pathway"


class GeneSpectraAdapterCellTypeField(Enum):
    """
    Define possible fields the adapter can provide for cell types.
    These are in the cell clade annotation file
    """

    SPECIES_OF_ORIGIN = "species_scientific_name"
    CELL_TYPE_ONTOLOGY_NAME = "cell_type_name"
    CELL_TYPE_ID = "cell_ontology_id"
    CELL_TYPE_NAME = "cell_type"
    CELL_TYPE_NAME_AND_SPECIES = "cell_type_name_and_species"
    BROAD_TYPE = "broad_type"  ## for MTG dataset cell broad types at various levels
    BROAD_TYPE_2 = "broad_type_2"
    BROAD_TYPE_3 = "broad_type_3"
    BROAD_TAXO_CS = "broad_taxo_cs"


class GeneSpectraAdapterOrthologousGroupField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    """

    ORTHOLOGOUS_GROUP_ID = "og_id"
    EGGNOG_DATASET_NAME = "eggnog_dataset_name"
    EGGNOG_DATASET_ID = "eggnog_dataset_id"
    ORTHOLOGOUS_GROUP_ID_AND_DATASET = "og_id_and_dataset"


class GeneSpectraAdapterEdgeType(Enum):
    """
    Enum for the types of the edges.
    """

    GENE_IN_ORTHOLOGOUS_GROUP = "gene_in_orthologous_group"
    GENE_WILCOX_MARKER_IN_CELL_TYPE = "gene_wilcox_marker_in_cell_type"


## fields relevant to the edges


class GeneSpectraAdapterEdgeField(Enum):
    """
    Enum for the fields of the edges.
    Thee are in the genespectra results

    """

    P_VAL = "pvals"
    AVG_LOG2FC = "logfoldchanges"  # log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
    P_VAL_ADJ = "pvals_adj"  # Adjusted p-value, based on bonferroni correction using all features in the dataset.
    P_VAL_ADJ_RANKING = "pvals_adj_rank"
    AVG_LOG2FC_RANKING = "logfc_rank"


class GeneSpectraAdapter:
    """
    GeneSpectra BioCypher adapter. Generates nodes and edges for creating a
    knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

        ## load cells and ids
        ## load eggnog ogs and genes memberships
        ## load genespectra data
        ## so far it is an all-in-one file

    def load_genespectra_data(
        self,
        eggnog_file=None,
        cell_ontology_file=None,
        wilcox_marker_file=None,
    ):
        """
        Read genespectra output csv and get dataframe for relevant fields.
        """
        logger.info("Loading data.")

        # read data

        self.eggnog = pd.read_csv(
            eggnog_file,
            delimiter=",",
            header=0,
            dtype=str,
        )

        self.eggnog = self.eggnog[
            [
                field.value
                for field in chain(
                    self.node_fields,
                    self.edge_fields,
                )
                if field.value in self.eggnog.columns
            ]
        ]

        self.eggnog = self.eggnog.loc[:, ~self.eggnog.columns.duplicated()]

        # often there are duplicated entrieds in eggnog as the same gene can have several different peptide products
        # in most cases the info we need from these dups are similar
        # so i just leep one for now... keep in mind that dedup should happen per-species

        self.eggnog = self.eggnog.drop_duplicates(
            subset=[
                GeneSpectraAdapterGeneField.PREFERRED_NAME.value,
                GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value,
            ]
        )

        self.wilcox_markers = pd.read_csv(
            wilcox_marker_file,
            sep=",",
            header=0,
            dtype=str,  # read all properties as str, nxbi_txid can get error for being int
        )

        self.wilcox_markers.replace("\\/", "-", regex=True, inplace=True)

        # filter only relevant fields
        self.wilcox_markers = self.wilcox_markers[
            [
                field.value
                for field in chain(
                    self.node_fields,
                    self.edge_fields,
                )
                if field.value in self.wilcox_markers.columns
            ]
        ]

        self.wilcox_markers = self.wilcox_markers.loc[
            :, ~self.wilcox_markers.columns.duplicated()
        ]

        self.cell_ontology = pd.read_csv(
            cell_ontology_file,
            sep=",",
            header=0,
            dtype=str,  # read all properties as str, nxbi_txid can get error for being int
        )

        # filter only relevant fields
        self.cell_ontology = self.cell_ontology[
            [
                field.value
                for field in chain(
                    self.node_fields,
                    self.edge_fields,
                )
                if field.value in self.cell_ontology.columns
            ]
        ]

        self.cell_ontology = self.cell_ontology.loc[
            :, ~self.cell_ontology.columns.duplicated()
        ]

        # replace the slashes in the cell ontology names to avoid potential neo4j prolblems
        self.cell_ontology.replace("\\/", "-", regex=True, inplace=True)

        # get relevant records
        # to include both genes in eggnog but not specific
        # and genes in genespectra but not in eggnog
        # we need to merge both dataframes

        self.gene = pd.merge(
            left=self.eggnog[
                [
                    field.value
                    for field in GeneSpectraAdapterGeneField
                    if field.value in self.eggnog.columns
                ]
            ].drop_duplicates(),
            right=self.wilcox_markers[
                [
                    GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX.value,
                    GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value,
                ]
            ].drop_duplicates(),
            left_on=[
                GeneSpectraAdapterGeneField.PREFERRED_NAME.value,
                GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value,
            ],
            right_on=[
                GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX.value,
                GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value,
            ],
            how="outer",
        ).drop_duplicates()

        # Give a unique idnetifier for gene nodes

        self.gene["gene_name_use"] = np.select(
            [
                self.gene[
                    GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX.value
                ].isna(),
                ~self.gene[
                    GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX.value
                ].isna(),
            ],
            [
                self.gene[GeneSpectraAdapterGeneField.PREFERRED_NAME.value],
                self.gene[GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX.value],
            ],
            default=None,
        )  # hard lessions learned from a long debug: the outer join can result in na values in gene names for those non-overlapping
        # and if i need the gene name for indexing, i need to combine the two gene names so there is no NA
        # if the node index col is na, biocypher will return a error that it is " a float value" and not a string

        self.gene[GeneSpectraAdapterGeneField.PREFERRED_NAME_AND_SPECIES.value] = (
            self.gene["gene_name_use"]
            + "_"
            + self.gene[GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value]
        )

        self.orthologous_group = self.eggnog[
            [
                field.value
                for field in GeneSpectraAdapterOrthologousGroupField
                if field.value in self.eggnog.columns
            ]
        ].drop_duplicates()

        # Give a unique identifier for orthologous group nodes
        self.orthologous_group[
            GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID_AND_DATASET.value
        ] = (
            self.orthologous_group[
                GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID.value
            ]
            + "_"
            + self.orthologous_group[
                GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_NAME.value
            ]
        )

        self.cell_type = self.cell_ontology[
            [
                field.value
                for field in GeneSpectraAdapterCellTypeField
                if field.value in self.cell_ontology.columns
            ]
        ].drop_duplicates()

        self.cell_type[
            GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME_AND_SPECIES.value
        ] = (
            self.cell_type[GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME.value]
            + "_"
            + self.cell_type[GeneSpectraAdapterCellTypeField.SPECIES_OF_ORIGIN.value]
        )

    def get_node_count(self):
        """
        Returns the number of nodes generated by the adapter.
        """
        return len(list(self.get_nodes()))

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in GeneSpectraAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    GeneSpectraAdapterCellTypeField,
                    GeneSpectraAdapterGeneField,
                    GeneSpectraAdapterOrthologousGroupField,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type_now for type_now in GeneSpectraAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field_now for field_now in GeneSpectraAdapterEdgeField]

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        logger.info("Generating nodes.")

        # GENES: for each node (row), yield a 3-tuple of node id (the 'NAME'
        # column with 'hgnc:' prefix), node label (hardcode to 'gene' for now),
        # and node properties (dict of column names and values, except the
        # 'NAME')

        self.nodes = []

        print("get gene nodes from genespectra and EggNOG")

        for _, node in self.gene.iterrows():
            yield (
                node[
                    GeneSpectraAdapterGeneField.PREFERRED_NAME_AND_SPECIES.value
                ],  # for this specifric dataset i need to use gene name and species name as encoding id
                "gene",
                {
                    "species_of_origin": node[
                        GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value
                    ],
                    "is_a_go_tf": node[GeneSpectraAdapterGeneField.IS_A_GO_TF.value],
                    "description": node[GeneSpectraAdapterGeneField.DESCRIPTION.value],
                    "preferred_name": node[
                        GeneSpectraAdapterGeneField.PREFERRED_NAME.value
                    ],
                    "perferred_name_wilcox": node[
                        GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX.value
                    ],
                    "pfams": node[GeneSpectraAdapterGeneField.PFAMS.value],
                    "gos": node[GeneSpectraAdapterGeneField.GOS.value],
                    "kegg_ko": node[GeneSpectraAdapterGeneField.KEGG_KO.value],
                    "peptide_id": node[GeneSpectraAdapterGeneField.PEPTIDE_ID.value],
                    "kegg_pathway": node[
                        GeneSpectraAdapterGeneField.KEGG_PATHWAY.value
                    ],
                },
            )
        print("finish writing gene nodes")
        print("get OG nodes from EggNOG")

        for _, node in self.orthologous_group.iterrows():
            yield (
                node[
                    GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID_AND_DATASET.value
                ],
                "orthologous_group",
                {
                    "eggnog_dataset_name": node[
                        GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_NAME.value
                    ],
                    "eggnog_dataset_id": node[
                        GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_ID.value
                    ],
                    "orthologous_group_id": node[
                        GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID.value
                    ],
                },
            )
        print("finish writing OG nodes")
        print("get cell type nodes from cell ontology info")

        for _, node in self.cell_type.iterrows():
            yield (
                node[GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME_AND_SPECIES.value],
                "cell_type",
                {
                    "cell_type_name": node[
                        GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME.value
                    ],
                    "species_of_origin": node[
                        GeneSpectraAdapterCellTypeField.SPECIES_OF_ORIGIN.value
                    ],
                    "cell_ontology_id": node[
                        GeneSpectraAdapterCellTypeField.CELL_TYPE_ID.value
                    ],
                    "broad_type": node[
                        GeneSpectraAdapterCellTypeField.BROAD_TYPE.value
                    ],
                    "broad_type_2": node[
                        GeneSpectraAdapterCellTypeField.BROAD_TYPE_2.value
                    ],
                    "broad_type_3": node[
                        GeneSpectraAdapterCellTypeField.BROAD_TYPE_3.value
                    ],
                    "broad_taxo_cs": node[
                        GeneSpectraAdapterCellTypeField.BROAD_TAXO_CS.value
                    ],
                    "cell_type_ontology_name": node[
                        GeneSpectraAdapterCellTypeField.CELL_TYPE_ONTOLOGY_NAME.value
                    ],
                },
            )

    def get_edges(self):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        """

        logger.info("Generating edges.")

        # if not self.nodes:
        #     raise ValueError("No nodes found. Please run get_nodes() first.")

        # yield 5-tuple of edge id, source node id, target node id, edge label
        # (hardcode to 'variant_in_gene' for now), and edge properties (empty
        # dict for now)

        print("Get gene belongs to OG edges")
        gene_from_orthologous_group_df = self.eggnog[
            [
                field.value
                for field in GeneSpectraAdapterGeneField
                if field.value in self.eggnog.columns
            ]
            + [
                field.value
                for field in GeneSpectraAdapterOrthologousGroupField
                if field.value in self.eggnog.columns
            ]
        ].drop_duplicates()

        gene_from_orthologous_group_df[
            GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID_AND_DATASET.value
        ] = (
            gene_from_orthologous_group_df[
                GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID.value
            ]
            + "_"
            + gene_from_orthologous_group_df[
                GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_NAME.value
            ]
        )

        gene_from_orthologous_group_df[
            GeneSpectraAdapterGeneField.PREFERRED_NAME_AND_SPECIES.value
        ] = (
            gene_from_orthologous_group_df[
                GeneSpectraAdapterGeneField.PREFERRED_NAME.value
            ]
            + "_"
            + gene_from_orthologous_group_df[
                GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value
            ]
        )

        for _, row in gene_from_orthologous_group_df.iterrows():
            gene_id = row[GeneSpectraAdapterGeneField.PREFERRED_NAME_AND_SPECIES.value]
            og_id = row[
                GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID_AND_DATASET.value
            ]
            _id = hashlib.md5((gene_id + og_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                gene_id,
                og_id,
                "gene_in_orthologous_group",
                {},
            )

        print("Get gene is a wilcox marker in cell type edges")
        gene_wilcox_marker_in_cell_type_df = self.wilcox_markers[
            [
                field.value
                for field in GeneSpectraAdapterGeneField
                if field.value in self.wilcox_markers.columns
            ]
            + [
                field.value
                for field in GeneSpectraAdapterCellTypeField
                if field.value in self.wilcox_markers.columns
            ]
            + [
                field.value
                for field in GeneSpectraAdapterEdgeField
                if field.value in self.wilcox_markers.columns
            ]
        ].drop_duplicates()

        gene_wilcox_marker_in_cell_type_df = gene_wilcox_marker_in_cell_type_df.loc[
            :, ~gene_wilcox_marker_in_cell_type_df.columns.duplicated()
        ]

        gene_wilcox_marker_in_cell_type_df[
            GeneSpectraAdapterGeneField.PREFERRED_NAME_AND_SPECIES.value
        ] = (
            gene_wilcox_marker_in_cell_type_df[
                GeneSpectraAdapterGeneField.PREFERRED_NAME_WILCOX.value
            ]
            + "_"
            + gene_wilcox_marker_in_cell_type_df[
                GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value
            ]
        )

        gene_wilcox_marker_in_cell_type_df[
            GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME_AND_SPECIES.value
        ] = (
            gene_wilcox_marker_in_cell_type_df[
                GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME.value
            ]
            + "_"
            + gene_wilcox_marker_in_cell_type_df[
                GeneSpectraAdapterGeneField.SPECIES_OF_ORIGIN.value
            ]
        )

        print("Yield gene is a wilcox marker in cell type edges")
        print(gene_wilcox_marker_in_cell_type_df.columns)
        print(gene_wilcox_marker_in_cell_type_df.head())
        print(gene_wilcox_marker_in_cell_type_df.shape)

        for _, row in gene_wilcox_marker_in_cell_type_df.iterrows():
            gene_id = row[GeneSpectraAdapterGeneField.PREFERRED_NAME_AND_SPECIES.value]
            ct_id = row[
                GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME_AND_SPECIES.value
            ]
            _id = hashlib.md5((gene_id + ct_id).encode("utf-8")).hexdigest()

            yield (
                _id,
                gene_id,
                ct_id,
                "gene_wilcox_marker_in_cell_type",
                {
                    "p_val": row[GeneSpectraAdapterEdgeField.P_VAL.value],
                    "avg_log2fc": row[GeneSpectraAdapterEdgeField.AVG_LOG2FC.value],
                    "p_val_adj": row[GeneSpectraAdapterEdgeField.P_VAL_ADJ.value],
                    "p_val_adj_ranking": row[
                        GeneSpectraAdapterEdgeField.P_VAL_ADJ_RANKING.value
                    ],
                    "avg_log2fc_ranking": row[
                        GeneSpectraAdapterEdgeField.AVG_LOG2FC_RANKING.value
                    ],
                },
            )


class Node:
    """
    Base class for nodes.
    """

    def __init__(self):
        self.id = None
        self.label = None
        self.properties = {}

    def get_id(self):
        """
        Returns the node id.
        """
        return self.id

    def get_label(self):
        """
        Returns the node label.
        """
        return self.label

    def get_properties(self):
        """
        Returns the node properties.
        """
        return self.properties
