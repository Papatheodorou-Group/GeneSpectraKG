import pandas as pd
import re
import biomart
from typing import Literal


def filter_and_extract(row, species_ids: list):
    """filter and extract the records relevant to the species of interest from eggnog6 raw data

    :param row: rows in the input table
    :type row:
    :param species_ids: a list of species ids (NCBI txid) to include, e.g. ['9606', '9544']
    :type species_ids: list
    :return: a row with the filtered species ids and sequence ids, or none if there is no relevant data in that row
    :rtype: _type_
    """
    second_last_col = row["species_id"].split(",")  # split by comma
    last_col = row["sequences_id"].split(",")  # split by comma

    if any(species_id in second_last_col for species_id in species_ids):
        second_last_col = [
            col
            for col in second_last_col
            if re.match(rf'^({"|".join(species_ids)})', col)
        ]  # select relevant records
        last_col = [
            col for col in last_col if re.match(rf'^({"|".join(species_ids)})\.', col)
        ]  # select relevant records
        row["species_id_filtered"] = second_last_col
        row["sequence_id_filtered"] = last_col
        return row
    else:
        row["species_id_filtered"] = None
        row["sequence_id_filtered"] = None
        return row


def read_eggnog_raw(
    file_path: str, eggnog_dataset_name: str, eggnog_version: Literal["v5", "v6"]
):
    """read eggnog raw data downloaded from the database

    :param file_path: path to the file, should be named e6.og2seqs_and_species.tsv and have 6 columns
    :param eggnog_dataset_name: name of the taxa that eggnog is calculated on, such as 'mammalia'
    :param eggnog_version: version of eggnog, either 'v5' or 'v6'
    :type file_path: str
    :return: pd.DataFrame
    :rtype: pd.DataFrame
    """
    if eggnog_version == "v6":
        df = pd.read_csv(
            file_path,
            delimiter="\t",
            header=None,
            names=[
                "eggnog_dataset_id",
                "eggnog_og_id",
                "num_species",
                "num_sequences",
                "species_id",
                "sequences_id",
            ],
            dtype=str,
        )
        df["eggnog_dataset_name"] = eggnog_dataset_name
        return df

    elif eggnog_version == "v5":
        df = pd.read_csv(
            file_path,
            delimiter="\t",
            header=None,
            names=[
                "eggnog_dataset_id",
                "eggnog_og_id",
                "num_sequences",
                "num_species",
                "sequences_id",
                "species_id",
            ],
            dtype=str,
        )
        df["eggnog_dataset_name"] = eggnog_dataset_name
        return df
    else:
        raise ValueError("eggnog_version should be either v5 or v6")


def filter_eggnog_ogs(og_members: pd.DataFrame, species_ids_keep=None):
    """filter eggnog ogs table by species ids

    :param og_members: the read-in raw eggnog6 file e6.og2seqs_and_species.tsv from read_eggnog6_raw
    :type og_members: pd.DataFrame
    :param species_ids_keep: the species ids to extract that are relevant to the analysis
    :type species_ids_keep: list
    :return:a filtered eggNOG ogs table
    :rtype: pd.DataFrame
    """

    res = pd.DataFrame(
        og_members.apply(filter_and_extract, species_ids=species_ids_keep, axis=1)
    )
    if res.shape[0] == 0:
        raise ValueError(
            "not all species are in EggNOG, this is quite unlikely so please check the species NCBI txids are correct"
        )
    eggnog_ogs = (
        res.dropna()
        .reset_index(drop=True)
        .drop(["species_id", "sequences_id", "species_id_filtered"], axis=1)
        .explode("sequence_id_filtered")
        .reset_index(drop=True)
    )
    eggnog_ogs[["species_id", "sequence_id"]] = eggnog_ogs[
        "sequence_id_filtered"
    ].str.split(".", expand=True)
    return eggnog_ogs


def eggnog_to_ensembl_gene(eggnog_ogs: pd.DataFrame, species_name_id: dict):
    """convert eggnog ogs to ensembl gene ids
    :param eggnog_ogs: filtered eggnog og table
    :type eggnog_ogs: pd.DataFrame
    :param species_name_id: a distionary mapping of species scientific name (for ensembl) and species id (for eggnog), should be only 1 species
    :type species_name_id: dict
    :return: a dataframe mapping eggnog records to ensembl gene ids
    :rtype: pd.DataFrame
    """
    species_id_value = str(list(species_name_id.values())[0])
    ## remember that after strsplit the species_id is a str

    species_name = str(list(species_name_id.keys())[0])

    print("Connecting to ensembl server")

    server = biomart.BiomartServer("http://www.ensembl.org/biomart")
    dataset = server.datasets[f"{species_name}_gene_ensembl"]

    attributes = ["ensembl_gene_id", "external_gene_name", "ensembl_peptide_id"]
    filters = {
        "ensembl_peptide_id": eggnog_ogs.loc[eggnog_ogs.species_id == species_id_value][
            "sequence_id"
        ].values
    }

    print("Downloading ensembl response")

    response = dataset.search(
        {"query": filters, "attributes": attributes, "mart_instance": "ensembl"}
    )

    gene_to_peptide = pd.read_csv(response.url, sep="\t", header=None, names=attributes)
    gene_to_peptide = gene_to_peptide.dropna()
    # re-filter because sometimes biomart query doesnt work
    gene_to_peptide = gene_to_peptide.loc[
        gene_to_peptide.ensembl_peptide_id.isin(
            eggnog_ogs.loc[eggnog_ogs.species_id == species_id_value][
                "sequence_id"
            ].values
        )
    ]

    print("Mapping peptides in EggNOG OGs to ensembl genes")

    eggnog_and_ensembl = gene_to_peptide.merge(
        eggnog_ogs.loc[eggnog_ogs.species_id == species_id_value],
        left_on="ensembl_peptide_id",
        right_on="sequence_id",
    )

    eggnog_and_ensembl["ncbi_txid"] = species_id_value
    eggnog_and_ensembl["species_scientific_name"] = species_name

    return eggnog_and_ensembl


def all_species_eggnog_to_ensembl(eggnog_ogs: pd.DataFrame, species_names_ids: dict):
    """convert all species eggnog ogs to ensembl gene ids

    :param eggnog_ogs: filtered eggnog og table
    :type eggnog_ogs: pd.DataFrame
    :param species_names_id: a distionary mapping of species scientific name (for ensembl) and species id (for eggnog), can have several species
    :type species_names_id: dict
    :return: a dataframe mapping eggnog records to ensembl gene ids
    :rtype: pd.DataFrame
    """
    all_species_eggnog_and_ensembl = pd.DataFrame()
    for species_name_id_now in species_names_ids.items():
        dict_now = dict({species_name_id_now[0]: species_name_id_now[1]})
        print(
            f"Processing species name {str(list(dict_now.keys())[0])}, EggNOG id {str(list(dict_now.values())[0])}"
        )
        species_eggnog_and_ensembl = eggnog_to_ensembl_gene(
            eggnog_ogs, species_name_id=dict_now
        )
        all_species_eggnog_and_ensembl = pd.concat(
            [all_species_eggnog_and_ensembl, species_eggnog_and_ensembl],
            ignore_index=True,
        )
        print(
            f"Finish species name {str(list(dict_now.keys())[0])}, EggNOG id {str(list(dict_now.keys())[0])}"
        )
    return all_species_eggnog_and_ensembl


def read_ncbi_gene2ensembl_raw(file_path):

    df = pd.read_csv(
        file_path,
        delimiter="\t",
        header=None,
        names=[
            "species_id",
            "ncbi_gene_id",
            "ensembl_gene_id",
            "ensembl_refseq_rna",
            "ensembl_rna_id",
            "ensembl_refseq_peptide",
            "ensembl_peptide_id",
        ],
        dtype=str,
    )
    df.ensembl_peptide_id = df.ensembl_peptide_id.apply(
        lambda x: re.sub("\..*$", "", x)
    )
    return df


def read_ncbi_gene2name_raw(file_path):
    # gene2accession file processed to spmplify using
    # awk '$1==9595 {print$1"\t"$2"\t"$16}' gene2accession.tsv | sort | uniq | grep -v Assembly | grep -v - > gorilla_gene_id_to_symbol.tsv
    return pd.read_csv(
        file_path,
        delimiter="\t",
        header=None,
        names=["species_id", "ncbi_gene_id", "ncbi_gene_name"],
        dtype=str,
    )


def eggnog_to_ncbi_gene(
    eggnog_ogs: pd.DataFrame,
    species_name_id: dict,
    ncbi_gene2ensembl: pd.DataFrame,
    ncbi_gene2accession: pd.DataFrame,
):
    species_id_value = str(list(species_name_id.values())[0])
    ## remember that after strsplit the species_id is a str

    species_name = str(list(species_name_id.keys())[0])

    ncbi_full = ncbi_gene2accession.merge(
        ncbi_gene2ensembl, left_on="ncbi_gene_id", right_on="ncbi_gene_id", how="right"
    )

    eggnog_and_ncbi = ncbi_full.merge(
        eggnog_ogs.loc[eggnog_ogs.species_id == species_id_value],
        left_on="ensembl_peptide_id",
        right_on="sequence_id",
    )

    eggnog_and_ncbi["ncbi_txid"] = species_id_value
    eggnog_and_ncbi["species_scientific_name"] = species_name

    return eggnog_and_ncbi


def all_species_eggnog_to_ncbi_ensembl(
    eggnog_ogs: pd.DataFrame, species_names_ids: dict, ncbi_files_dir: str
):
    all_species_eggnog_and_ncbi_ensembl = pd.DataFrame()
    for species_name_id_now in species_names_ids.items():
        dict_now = dict({species_name_id_now[0]: species_name_id_now[1]})
        species_name_now = str(list(dict_now.keys())[0])
        species_id_now = str(list(dict_now.values())[0])
        print(f"Processing species name {species_name_now}, EggNOG id {species_id_now}")
        gene2ensembl_now = read_ncbi_gene2ensembl_raw(
            file_path=f"{ncbi_files_dir}/{species_name_now}_gene2ensembl.tsv"
        )
        gene2name_now = read_ncbi_gene2name_raw(
            file_path=f"{ncbi_files_dir}/{species_name_now}_gene_id_to_symbol.tsv"
        )
        mapped_now = eggnog_to_ncbi_gene(
            eggnog_ogs=eggnog_ogs,
            species_name_id=dict_now,
            ncbi_gene2accession=gene2name_now,
            ncbi_gene2ensembl=gene2ensembl_now,
        )

        all_species_eggnog_and_ncbi_ensembl = pd.concat(
            [all_species_eggnog_and_ncbi_ensembl, mapped_now], ignore_index=True
        )

        print(f"Finishing species name {species_name_now}, EggNOG id {species_id_now}")

    return all_species_eggnog_and_ncbi_ensembl


def process_emapper_results(emapper_file_path, species_name, species_scientific_name):
    """
    This is for processing emapper results to get gene functional annotations
    Contains all OG annotations of all genes, and their GO, KEGG, PFAM annotations etc.
    Processed from emapper out.emapper.annotations

    :param emapper_file_path: path to briefly processed emapper file
    :type emapper_file_path: str
    :param species_name: convenient name used for the species file name output
    :type species_name: str
    :param species_scientific_name: species scientific name for the KG
    :type species_scientific_name: str
    """
    eggnog = pd.read_csv(emapper_file_path)

    eggnog = eggnog.assign(ogs=eggnog.eggNOG_OGs.str.split(",")).explode("ogs")
    eggnog[["og_name", "eggnog_dataset_name"]] = eggnog["ogs"].str.split(
        "|", n=1, expand=True
    )
    eggnog[["og_id", "eggnog_dataset_id"]] = eggnog["og_name"].str.split(
        "@", n=1, expand=True
    )
    eggnog.rename(columns={"query": "peptide"}, inplace=True)

    eggnog["species_of_origin"] = species_scientific_name

    eggnog = eggnog.loc[
        ~eggnog.gene.str.contains(
            "|".join(["queries scanned", "## Rate:", "## Total time \(seconds\):"])
        ),
    ]

    eggnog.replace("3'", "3 prime", regex=True, inplace=True)
    eggnog.replace("2'", "2 prime", regex=True, inplace=True)
    eggnog.replace("5'", "5 prime", regex=True, inplace=True)
    eggnog.replace("'", "prime", regex=True, inplace=True)

    eggnog.to_csv(
        f"../data/wilcox_{species_name}_annot_skim_OGs_processed.csv", index=False
    )
