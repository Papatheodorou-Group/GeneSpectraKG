# The GeneSpectraKG is built using the BioCypther project template
Yuyao Song <ysong@ebi.ac.uk>

# Description

The GeneSpectraKG is a method to generate a knowledge graph that integrates GeneSpectra results of gene classification with orthologous groups, species, and cell types. The knowledge graph describes the cell type pattern of gene expression across species. It is built using the BioCypher project template, which provides a framework for creating and managing knowledge graphs in the life sciences.

## Adapters

In `genespectrakg/adapters`, there are adapters for generating a BioCypher knowledge graph from GeneSpectra results. The adapters are designed to convert the output of GeneSpectra into nodes and edges that can be made into a BioCypher knowledge graph.

`genespectra_adapter_individual.py`: generate node and edges from [GeneSpectra](https://github.com/Papatheodorou-Group/GeneSpectra) gene classification results
`genespectra_adapter_wilcox_logit.py`: generate node and edges from marker calling using Wilcoxon test or the Logistic regression method 

 - Nodes:
    - Genes
    - Orthologs Groups
    - Species
    - Cell Types

- Edges:
    - Gene to Ortholog Group
    - Gene to Species
    - Cell type to Species
    - Gene to to Cell Type

## Config

Configs are stored in `genespectrakg/config`. The config files define the structure of the knowledge graph, including the graph schema, and BioCypher settings.

## Notebooks

Practical examples of how to use the adapters and configs to generate a BioCypher knowledge graph from GeneSpectra results are provided in `genespectrakg/notebooks`.

 - Prepare the data as input for the adapters (process results from GeneSpectra, emapper, marker calling etc.)
 - Create the KG
 - Query the Local KG with the neo4j python driver


# BioCypher project template

A quick way to set up a BioCypher-driven knowledge graph pipeline.

Please refer to [BioCypher project template](https://github.com/biocypher/project-template) for more information.