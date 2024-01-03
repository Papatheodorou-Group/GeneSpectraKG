def split_by_comma(row):
    row['enriched_group_split'] = row['enriched_groups'].split(";")
    return row

def determine_enriched_enhanced(row):
    if row['spec_category'] in ['cell type enriched', 'group enriched']:
        row['specificity_category'] = 'enriched'
    if row['spec_category'] in ['cell type enhanced', 'group enhanced']:
        row['specificity_category'] = 'enhanced'
    return row

def preprocess_genespectra_output(gene_classes):
    result_specific = gene_classes.loc[gene_classes.spec_category.isin(['cell type enriched', 'cell type enhanced', 'group enriched', 'group enhanced'])]
    result_specific = result_specific.apply(split_by_comma, axis=1)
    result_specific = result_specific.apply(determine_enriched_enhanced, axis=1)
    results_new = result_specific.explode(['enriched_group_split']).reset_index(drop=True)
    results_few = results_new[['gene', 'specificity_category', 'spec_category', 'spec_score', 'dist_category', 'n_exp', 'mean_exp', 'max_exp', 'frac_exp', 'enriched_group_split']]
    results_few.columns = ['external_gene_name', 'specificity_category_type', 'specificity_category', 'specificity_score', 'distribution_catehory', 'n_expressed', 'mean_expression', 'max_expression', 'fraction_expressed', 'cell_type_name']
    return results_few

