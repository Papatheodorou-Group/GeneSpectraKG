import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger
import pandas as pd
import hashlib

logger.debug(f"Loading module {__name__}.")

# reference:
#https://github.com/biocypher/decider-genetics/blob/main/decider_genetics/adapters/cn_genes_adapter.py

class GeneSpectraAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    GENE = auto()
    CELL_TYPE = auto()
    SPECIES = auto()
    ORTHOLOGOUS_GROUP = auto()


class GeneSpectraAdapterGeneField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    The names correspond to the column names in the genespectra input file
    """

    GENE_ID = "ensembl_gene_id"
    # GENE_NAME = "ncbi_gene_name"
    GENE_NAME = "external_gene_name"


class GeneSpectraAdapterCellTypeField(Enum):
    """
    Define possible fields the adapter can provide for cell types.
    """

    CELL_TYPE_ID = "cell_ontology_id"
    CELL_TYPE_NAME = "cell_type_name"
    TISSUE_ID = "uberon_tissue_id"
    TISSUE_NAME = "tissue_name"

class GeneSpectraAdapterSpeciesField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    """

    SPECIES_ID = "ncbi_txid"
    SPECIES_NAME = "species_scientific_name"


class GeneSpectraAdapterOrthologousGroupField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    """

    ORTHOLOGOUS_GROUP_ID = "eggnog_og_id"
    EGGNOG_DATASET_NAME = "eggnog_dataset_name"
    EGGNOG_DATASET_ID = 'eggnog_dataset_id'


class GeneSpectraAdapterEdgeType(Enum):
    """
    Enum for the types of the edges.
    """

    GENE_IN_ORTHOLOGOUS_GROUP = "gene_in_orthologous_group"
    GENE_FROM_SPECIES = "gene_from_species"
    CELL_TYPE_FROM_SPECIES = "cell_type_from_species"
    GENE_ENRICHED_IN_CELL_TYPE = "gene_enriched_in_cell_type"
    GENE_ENHANCED_IN_CELL_TYPE = "gene_enhanced_in_cell_type"


## fields relevant to the edges

class GeneSpectraAdapterEdgeField(Enum):
    """
    Enum for the fields of the edges.
        
    """

    SPECIFICITY_CATEGORY_TYPE = 'specificity_category_type'
    SPECIFICITY_CATEGORY = 'specificity_category'
    SPECIFICITY_SCORE = 'specificity_score'
    DISTRIBUTION_CATEGORY = 'distribution_catehory'
    NUMBER_EXPRESSED = 'n_expressed'
    MEAN_EXPRESSION = 'mean_expression'
    MAX_EXPRESSION = 'max_expression'
    FRACTION_EXPRESSED = 'fraction_expressed'



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

    def load_genespectra_data(self, eggnog_file=None, cell_ontology_file=None, genespectra_file=None):

        """
        Read genespectra output csv and get dataframe for relevant fields.
        """
        logger.info("Loading data.")

        # read data

        self.eggnog = pd.read_csv(
            eggnog_file, delimiter=',', header=0, dtype=str,
        
        )
        
        self.eggnog = self.eggnog[
            [
                field.value
                for field in chain(
                    self.node_fields,
                )
                if field.value in self.eggnog.columns
            ]
        ]



        self.genespectra = pd.read_csv(
            genespectra_file, sep=",", header=0, dtype=str,# read all properties as str, nxbi_txid can get error for being int
        )

        # filter only relevant fields
        self.genespectra = self.genespectra[
            [
                field.value
                for field in chain(
                    self.node_fields,
                    self.edge_fields,
                )
                if field.value in self.genespectra.columns
            ]
        ]

        self.genespectra_enriched = self.genespectra.loc[self.genespectra.specificity_category_type == 'enriched'].drop_duplicates()
        self.genespectra_enhanced = self.genespectra.loc[self.genespectra.specificity_category_type  == 'enhanced'].drop_duplicates()

        self.cell_ontology = pd.read_csv(
            cell_ontology_file, sep=",", header=0, dtype=str, # read all properties as str, nxbi_txid can get error for being int
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

        # get relevant records

        self.gene = self.eggnog[
            [
                field.value
                for field in GeneSpectraAdapterGeneField
                if field.value in self.eggnog.columns
            ]
        ].drop_duplicates()

        self.orthologous_group = self.eggnog[
            [
                field.value
                for field in GeneSpectraAdapterOrthologousGroupField
                if field.value in self.eggnog.columns
            ]
        ].drop_duplicates()

        self.species = self.eggnog[
            [
                field.value
                for field in GeneSpectraAdapterSpeciesField
                if field.value in self.eggnog.columns
            ]
        ].drop_duplicates()

        self.cell_type = self.cell_ontology[
            [
                field.value
                for field in GeneSpectraAdapterCellTypeField
                if field.value in self.cell_ontology.columns
            ]
        ].drop_duplicates()



    # def get_nodes_eggnog(self):
    #     """generate gene, orthologous group and species nodes from eggnog data

    #     :yield: _description_
    #     :rtype: _type_
    #     """

        


    # def get_nodes_cell_ontology(self):
    #     """generate cell types from cell ontology data

    #     :yield: _description_
    #     :rtype: _type_
    #     """
    #     print('get nodes cell ontology')

    #     for _, node in self.cell_type.iterrows():
    #         yield (
    #             node[GeneSpectraAdapterCellTypeField.CELL_TYPE_ID.value],
    #             "cell_type",
    #             {"cell_type_name": node[GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME.value],
    #             "uberon_tissue_id": node[GeneSpectraAdapterCellTypeField.TISSUE_ID.value],
    #             "tissue_name": node[GeneSpectraAdapterCellTypeField.TISSUE_NAME.value],},
    #         )




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

        print('get gene nodes from EggNOG')

        for _, node in self.gene.iterrows():
            yield (
                node[GeneSpectraAdapterGeneField.GENE_ID.value],
                "gene",
                {"external_gene_name": node[GeneSpectraAdapterGeneField.GENE_NAME.value],},
            )

        print('get OG nodes from EggNOG')

        for _, node in self.orthologous_group.iterrows():
                yield (
                    node[GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID.value],
                    "orthologous_group",
                    {"eggnog_dataset_name": node[GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_NAME.value],
                    "eggnog_dataset_id": node[GeneSpectraAdapterOrthologousGroupField.EGGNOG_DATASET_ID.value],},
                )

        print('get species nodes from EggNOG')

        for _, node in self.species.iterrows():
            yield (
                node[GeneSpectraAdapterSpeciesField.SPECIES_ID.value],
                "species",
                {"species_scientific_name": node[GeneSpectraAdapterSpeciesField.SPECIES_NAME.value],},
            )

        print('get cell type nodes from cell ontology info')

        for _, node in self.cell_type.iterrows():
            yield (
                node[GeneSpectraAdapterCellTypeField.CELL_TYPE_ID.value],
                "cell_type",
                {"cell_type_name": node[GeneSpectraAdapterCellTypeField.CELL_TYPE_NAME.value],
                "uberon_tissue_id": node[GeneSpectraAdapterCellTypeField.TISSUE_ID.value],
                "tissue_name": node[GeneSpectraAdapterCellTypeField.TISSUE_NAME.value],},
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



        print('Get cell type from species edges')

        ct_from_species_df = self.cell_ontology[
            [
                field.value
                for field in GeneSpectraAdapterCellTypeField
                if field.value in self.cell_ontology.columns
            ] +
            [
                field.value
                for field in GeneSpectraAdapterSpeciesField
                if field.value in self.cell_ontology.columns
            ]
        ].drop_duplicates()

        for _, row in ct_from_species_df.iterrows():
            ct_id = row[GeneSpectraAdapterCellTypeField.CELL_TYPE_ID.value]
            sp_id = row[GeneSpectraAdapterSpeciesField.SPECIES_ID.value]
            _id = hashlib.md5((ct_id + sp_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                ct_id,
                sp_id,
                "cell_type_from_species",
                {},
            )
    
    
        print('Get gene from species edges')
        gene_from_species_df = self.eggnog[
            [
                field.value
                for field in GeneSpectraAdapterGeneField
                if field.value in self.eggnog.columns
            ] +
            [
                field.value
                for field in GeneSpectraAdapterSpeciesField
                if field.value in self.eggnog.columns
            ]
        ].drop_duplicates()


        for _, row in gene_from_species_df.iterrows():
            gene_id = row[GeneSpectraAdapterGeneField.GENE_ID.value]
            sp_id = row[GeneSpectraAdapterSpeciesField.SPECIES_ID.value]
            _id = hashlib.md5((gene_id + sp_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                gene_id,
                sp_id,
                "gene_from_species",
                {},
            )


    
        print('Get gene from OG edges')
        gene_from_orthologous_group_df = self.eggnog[
            [
                field.value
                for field in GeneSpectraAdapterGeneField
                if field.value in self.eggnog.columns
            ] +
            [
                field.value
                for field in GeneSpectraAdapterOrthologousGroupField
                if field.value in self.eggnog.columns
            ]
        ].drop_duplicates()


        for _, row in gene_from_orthologous_group_df.iterrows():
            gene_id = row[GeneSpectraAdapterGeneField.GENE_ID.value]
            og_id = row[GeneSpectraAdapterOrthologousGroupField.ORTHOLOGOUS_GROUP_ID.value]
            _id = hashlib.md5((gene_id + og_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                gene_id,
                og_id,
                "gene_in_orthologous_group",
                {},
            )


    
        print('Get gene enriched in cell type edges')
        gene_enriched_in_cell_type_df = self.genespectra_enriched[
            [
                field.value
                for field in GeneSpectraAdapterGeneField
                if field.value in self.genespectra_enriched.columns
            ] +
            [
                field.value
                for field in GeneSpectraAdapterCellTypeField
                if field.value in self.genespectra_enriched.columns
            ] +
            [
                field.value
                for field in GeneSpectraAdapterEdgeField
                if field.value in self.genespectra_enriched.columns
            ]
        ].drop_duplicates()

        for _, row in gene_enriched_in_cell_type_df.iterrows():
            gene_id = row[GeneSpectraAdapterGeneField.GENE_ID.value]
            ct_id = row[GeneSpectraAdapterCellTypeField.CELL_TYPE_ID.value]
            _id = hashlib.md5((gene_id + ct_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                gene_id,
                ct_id,
                "gene_enriched_in_cell_type",
                {"specificity_category": row[GeneSpectraAdapterEdgeField.SPECIFICITY_CATEGORY.value],
                "specificity_score": row[GeneSpectraAdapterEdgeField.SPECIFICITY_SCORE.value],
                "distribution_category": row[GeneSpectraAdapterEdgeField.DISTRIBUTION_CATEGORY.value],
                "max_expression" : row[GeneSpectraAdapterEdgeField.MAX_EXPRESSION.value],
                "mean_expression": row[GeneSpectraAdapterEdgeField.MEAN_EXPRESSION.value],
                "number_expressed": row[GeneSpectraAdapterEdgeField.NUMBER_EXPRESSED.value],
                "fraction_expressed": row[GeneSpectraAdapterEdgeField.FRACTION_EXPRESSED.value],},
            )

    
        print('Get gene enhanced in cell type edges')
        gene_enhanced_in_cell_type_df = self.genespectra_enhanced[
            [
                field.value
                for field in GeneSpectraAdapterGeneField
                if field.value in self.genespectra_enhanced.columns
            ] +
            [
                field.value
                for field in GeneSpectraAdapterCellTypeField
                if field.value in self.genespectra_enhanced.columns
            ]+
            [
                field.value
                for field in GeneSpectraAdapterEdgeField
                if field.value in self.genespectra_enriched.columns
            ]
        ].drop_duplicates()

        for _, row in gene_enhanced_in_cell_type_df.iterrows():
            gene_id = row[GeneSpectraAdapterGeneField.GENE_ID.value]
            ct_id = row[GeneSpectraAdapterCellTypeField.CELL_TYPE_ID.value]
            _id = hashlib.md5((gene_id + ct_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                gene_id,
                ct_id,
                "gene_enhanced_in_cell_type",
                {"specificity_category": row[GeneSpectraAdapterEdgeField.SPECIFICITY_CATEGORY.value],
                "specificity_score": row[GeneSpectraAdapterEdgeField.SPECIFICITY_SCORE.value],
                "distribution_category": row[GeneSpectraAdapterEdgeField.DISTRIBUTION_CATEGORY.value],
                "max_expression" : row[GeneSpectraAdapterEdgeField.MAX_EXPRESSION.value],
                "mean_expression": row[GeneSpectraAdapterEdgeField.MEAN_EXPRESSION.value],
                "number_expressed": row[GeneSpectraAdapterEdgeField.NUMBER_EXPRESSED.value],
                "fraction_expressed": row[GeneSpectraAdapterEdgeField.FRACTION_EXPRESSED.value],},
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
                    GeneSpectraAdapterSpeciesField
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


# class Protein(Node):
#     """
#     Generates instances of proteins.
#     """

#     def __init__(self, fields: Optional[list] = None):
#         self.fields = fields
#         self.id = self._generate_id()
#         self.label = "uniprot_protein"
#         self.properties = self._generate_properties()

#     def _generate_id(self):
#         """
#         Generate a random UniProt-style id.
#         """
#         lets = [random.choice(string.ascii_uppercase) for _ in range(3)]
#         nums = [random.choice(string.digits) for _ in range(3)]

#         # join alternating between lets and nums
#         return "".join([x for y in zip(lets, nums) for x in y])

#     def _generate_properties(self):
#         properties = {}

#         ## random amino acid sequence
#         if (
#             self.fields is not None
#             and GeneSpectraAdapterProteinField.SEQUENCE in self.fields
#         ):
#             # random int between 50 and 250
#             l = random.randint(50, 250)

#             properties["sequence"] = "".join(
#                 [random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(l)],
#             )

#         ## random description
#         if (
#             self.fields is not None
#             and GeneSpectraAdapterProteinField.DESCRIPTION in self.fields
#         ):
#             properties["description"] = " ".join(
#                 [random.choice(string.ascii_lowercase) for _ in range(10)],
#             )

#         ## taxon
#         if self.fields is not None and GeneSpectraAdapterProteinField.TAXON in self.fields:
#             properties["taxon"] = "9606"

#         return properties


# class Disease(Node):
#     """
#     Generates instances of diseases.
#     """

#     def __init__(self, fields: Optional[list] = None):
#         self.fields = fields
#         self.id = self._generate_id()
#         self.label = "do_disease"
#         self.properties = self._generate_properties()

#     def _generate_id(self):
#         """
#         Generate a random disease id.
#         """
#         nums = [random.choice(string.digits) for _ in range(8)]

#         return f"DOID:{''.join(nums)}"

#     def _generate_properties(self):
#         properties = {}

#         ## random name
#         if self.fields is not None and GeneSpectraAdapterDiseaseField.NAME in self.fields:
#             properties["name"] = " ".join(
#                 [random.choice(string.ascii_lowercase) for _ in range(10)],
#             )

#         ## random description
#         if (
#             self.fields is not None
#             and GeneSpectraAdapterDiseaseField.DESCRIPTION in self.fields
#         ):
#             properties["description"] = " ".join(
#                 [random.choice(string.ascii_lowercase) for _ in range(10)],
#             )

#         return properties
