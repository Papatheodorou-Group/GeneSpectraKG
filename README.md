# GeneSpectraKG: a knowledge graph to describe the cell type pattern of gene expression across species

Yuyao Song <ysong@ebi.ac.uk>

**The GeneSpectraKG is built using the BioCypther project template**


# Description

The GeneSpectraKG is a method to generate a knowledge graph that integrates [GeneSpectra](https://github.com/Papatheodorou-Group/GeneSpectra) results of gene classification with orthologous groups, species, and cell types. This knowledge graph can be used to explore gene expression patterns across species and cell types, and to quantify expression similarity between cell types leveraging all types of orthologous genes (one-to-one orthologs and non-one-to-one orthologs).

## Adapters

In `genespectrakg/adapters`, there are adapters for generating a BioCypher knowledge graph from GeneSpectra results. The adapters are designed to convert the output of GeneSpectra into nodes and edges that can be made into a BioCypher knowledge graph.

 - Nodes:
    - Genes (from ENSEMBL or EggNOG)
    - Orthologs Groups (from EggNOG)
    - Species (from single cell data)
    - Cell Types (from single cell data)

- Edges:
    - Gene to Ortholog Group (from EggNOG)
    - Gene to Species (from ENSEMBL or EggNOG)
    - Cell type to Species (from single cell data)
    - Gene to to Cell Type (various methods for expression specificity)

`genespectra_adapter_individual.py`: generate node and edges from [GeneSpectra](https://github.com/Papatheodorou-Group/GeneSpectra) gene classification results
`genespectra_adapter_wilcox_logit.py`: generate node and edges from marker calling using Wilcoxon test or the Logistic regression method.


## Config

Configs are stored in `genespectrakg/config`. The config files define the structure of the knowledge graph, including the graph schema, and BioCypher settings.

## Notebooks

Practical examples of how to use the adapters and configs to generate a BioCypher knowledge graph from GeneSpectra results are provided in `genespectrakg/notebooks`.

 - Prepare the data as input for the adapters (process results from GeneSpectra, emapper, marker calling etc.)
 - Create the KG using the adapters and configs
 - Query the Local KG with the neo4j python driver


## BioCypher project template

A quick way to set up a BioCypher-driven knowledge graph pipeline.

Please refer to [BioCypher project template](https://github.com/biocypher/project-template) for more information.