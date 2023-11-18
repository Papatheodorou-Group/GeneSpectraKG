import pandas as pd
import re
import biomart

def filter_and_extract(row, species_ids: list):
    """filter and extract the records relevant to the species of interest

    :param row: rows in the input table
    :type row: _type_
    :param species_ids: _description_
    :type species_ids: list
    :return: _description_
    :rtype: _type_
    """
    second_last_col = row['species_id'].split(",")  # Split the second-to-last column by comma
    last_col = row['sequences_id'].split(",")  # Split the last column by comma

    if all(species_id in second_last_col for species_id in species_ids):
        second_last_col = [col for col in second_last_col if re.match(fr'^({"|".join(species_ids)})', col)]
        last_col = [col for col in last_col if re.match(fr'^({"|".join(species_ids)})\.', col)]
        row['species_id_filtered'] = second_last_col
        row['sequence_id_filtered'] = last_col
        return row
    else:
        row['species_id_filtered'] = None
        row['sequence_id_filtered'] = None
        return row
    
def read_eggnog6_raw(file_path: str):

    return pd.read_csv(file_path, delimiter='\t', header=None, names=['taxa_id', 'og_id', 'num_species', 'num_sequences', 'species_id', 'sequences_id'])

    
def filter_eggnog_ogs(og_members: pd.DataFrame, species_ids_keep=None):
    """filter eggnog ogs table by species ids

    :param og_members: _description_
    :type og_members: _type_
    :param species_ids_keep: _description_, defaults to ["9544", "9606"]
    :type species_ids_keep: list, optional
    :return: _description_
    :rtype: _type_
    """

    res = pd.DataFrame(og_members.apply(filter_and_extract, species_ids = species_ids_keep, axis=1))
    eggnog_ogs = res.dropna().reset_index(drop=True).drop(['species_id', 'sequences_id', 'species_id_filtered'], axis=1).explode('sequence_id_filtered').reset_index(drop=True)
    eggnog_ogs[['species_id', 'sequence_id']] = eggnog_ogs['sequence_id_filtered'].str.split('.',expand=True)
    return eggnog_ogs


def eggnog_to_ensembl_gene(eggnog_ogs: pd.DataFrame, species_name_id: dict):
    species_id_value = str(list(species_name_id.values())[0])
    ## remember that after strsplit the species_id is a str
    species_name = str(list(species_name_id.keys())[0])
    print(f"processing species {species_name}, eggNOG id {species_id_value}")
    server = biomart.BiomartServer("http://www.ensembl.org/biomart")
    dataset = server.datasets[f'{species_name}_gene_ensembl']
    
    attributes = ['ensembl_gene_id', 'external_gene_name', 'ensembl_peptide_id']
    filters = {'ensembl_peptide_id': eggnog_ogs.loc[eggnog_ogs.species_id == species_id_value]['sequence_id'].values}  
    response = dataset.search({'query': filters,
            'attributes': attributes,
            'mart_instance': 'ensembl'})
    
    gene_to_peptide = pd.read_csv(response.url, sep='\t', header=None, names=attributes)
    gene_to_peptide = gene_to_peptide.dropna()
    # re-filter because sometimes biomart query doesnt work
    gene_to_peptide = gene_to_peptide.loc[gene_to_peptide.ensembl_peptide_id.isin(eggnog_ogs.loc[eggnog_ogs.species_id == species_id_value]['sequence_id'].values)]
    
    eggnog_and_ensembl = gene_to_peptide.merge(eggnog_ogs.loc[eggnog_ogs.species_id == species_id_value], left_on='ensembl_peptide_id', right_on='sequence_id')
    
    return eggnog_and_ensembl


def all_species_eggnog_to_ensembl(eggnog_ogs: pd.DataFrame, species_names_ids: dict):
    all_species_eggnog_and_ensembl = pd.DataFrame()
    for species_name_id_now in species_names_ids.items():
        print(species_name_id_now)
        species_eggnog_and_ensembl = eggnog_to_ensembl_gene(eggnog_ogs, dict({species_name_id_now[0]: species_name_id_now[1]}))
        all_species_eggnog_and_ensembl = all_species_eggnog_and_ensembl.append(species_eggnog_and_ensembl, ignore_index=True)
    return all_species_eggnog_and_ensembl




