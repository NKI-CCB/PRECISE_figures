"""
@author: Soufiane Mourragui
"""

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()



def read_cell_line_data(cell_line_file,\
                        gene_lookup_file=None,\
                        cell_line_lookup_file=None,\
                        tumor_type=None):
    
    #Read data
    dataRDS = robjects.r['readRDS']
    df_data = dataRDS(cell_line_file)

    #read the gene expression, the samples' and genes' name
    numpy_formatted_data = pandas2ri.ri2py(df_data)
    samples_name = np.array(df_data.names[0])
    genes_name = np.array(df_data.names[1])

    if cell_line_lookup_file is not None and tumor_type is not None:
        #Load the data types
        cell_lines_cancer_types = pd.read_csv(cell_line_lookup_file, delimiter='\t')
        relevant_samples = cell_lines_cancer_types[cell_lines_cancer_types.tcga_type == tumor_type]['sample_name'].values

        #Filter on type
        relevant_samples_index = np.where(np.isin(samples_name, relevant_samples))[0]
        numpy_formatted_data = numpy_formatted_data[relevant_samples_index,:]
        samples_name = samples_name[relevant_samples_index]

    if gene_lookup_file is not None:
        #Remove non-protein coding genes
        genes_lookup_table = pd.read_csv(gene_lookup_file)
        genes_lookup_table = genes_lookup_table[genes_lookup_table.status == 'protein_coding']
        protein_coding_genes = np.array(genes_lookup_table.TCGA_name.values)

        #Filter and retrieve data
        protein_coding_genes_index = np.where(np.isin(genes_name, protein_coding_genes))[0]
        filtered_formatted_data = numpy_formatted_data[:,protein_coding_genes_index]
        filtered_genes = genes_name[protein_coding_genes_index]

        return filtered_formatted_data, filtered_genes, samples_name
    
    return numpy_formatted_data, genes_name, samples_name


def read_cell_line_data_cell_passport(cell_line_file,\
                                gene_lookup_file=None,\
                                cell_line_lookup_file=None,\
                                tumor_type=None):
    # Read Cell Model Passport
    cell_model_df = pd.read_csv(cell_line_file, sep=',')
    cell_model_gene_df = pd.read_csv(gene_lookup_file, sep=',')
    cell_model_sample_df = pd.read_csv(cell_line_lookup_file, sep=',')

    if tumor_type is not None:
        cell_model_sample_df = cell_model_sample_df[cell_model_sample_df.tissue == tumor_type]

    cell_model_df = cell_model_df.merge(cell_model_sample_df[['model_id', 'model_name']], on='model_id', how='inner')
    cell_model_df = cell_model_df.merge(cell_model_gene_df[['gene_id', 'ensembl_gene_id']], on='gene_id', how='inner')
    cell_model_df = cell_model_df[['model_name', 'ensembl_gene_id', 'read_count', 'fpkm']]

    gene_caract_file = '/DATA/s.mourragui/data/2019_04_cell_line_data/pybiomart_gene_status.csv'
    gene_lookup_df = pd.read_csv(gene_caract_file, sep='\t')[['ENSEMBL','Hugo', 'status']]
    gene_lookup_df = gene_lookup_df.drop_duplicates(subset=['Hugo'], keep=False)
    cell_model_df = cell_model_df.merge(gene_lookup_df, left_on='ensembl_gene_id', right_on='ENSEMBL')
    cell_model_df = cell_model_df[cell_model_df.status == 'protein_coding']
    cell_model_df.head()

    # Put data in data matrix (samples per genes)
    cell_model_df = cell_model_df[['model_name', 'ensembl_gene_id', 'read_count']]
    cell_model_data_df = cell_model_df.pivot_table(columns=['ensembl_gene_id'], index='model_name', fill_value=0)

    source_data = cell_model_data_df.values
    source_gene_names = np.array(cell_model_data_df.columns.droplevel()).astype(str)
    source_samples = np.array(cell_model_data_df.index).astype(str)

    return source_data, source_gene_names, source_samples


def read_cell_line_data_fpkm(cell_line_file,\
                        gene_lookup_file=None,\
                        cell_line_lookup_file=None,\
                        tumor_type=None):

    # Read data
    df_data = pd.read_csv(cell_line_file, sep='\t', index_col=0)

    # Read the gene expression, the samples' and genes' name
    numpy_formatted_data = np.array(df_data)
    samples_name = np.array(df_data.index).astype(str)
    genes_name = np.array(df_data.columns).astype(str)

    if cell_line_lookup_file is not None and tumor_type is not None:
        #Load the data types
        cell_lines_cancer_types = pd.read_csv(cell_line_lookup_file, delimiter='\t')
        
        relevant_samples = cell_lines_cancer_types[cell_lines_cancer_types.tcga_type == tumor_type]['sample_name'].values

        #Filter on type
        relevant_samples_index = np.where(np.isin(samples_name, relevant_samples))[0]
        numpy_formatted_data = numpy_formatted_data[relevant_samples_index,:]
        samples_name = samples_name[relevant_samples_index]

    if gene_lookup_file is not None:
        #Remove non-protein coding genes
        genes_lookup_table = pd.read_csv(gene_lookup_file)
        genes_lookup_table = genes_lookup_table[genes_lookup_table.status == 'protein_coding']
        protein_coding_genes = np.array(genes_lookup_table.TCGA_name.values)

        #Filter and retrieve data
        protein_coding_genes_index = np.where(np.isin(genes_name, protein_coding_genes))[0]
        filtered_formatted_data = numpy_formatted_data[:,protein_coding_genes_index]
        filtered_genes = genes_name[protein_coding_genes_index]

        return filtered_formatted_data, filtered_genes, samples_name

    return numpy_formatted_data, genes_name, samples_name
