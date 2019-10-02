# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

READ_PDX_DATA
"""

import numpy as np
import pandas as pd


def read_pdx_data(pdx_file,
                gene_lookup_file=None,
                cell_line_lookup_file=None,
                tumor_type=None):
    
    #Read data
    pdx_data = pd.read_csv(pdx_file)
    pdx_data = pdx_data.sort_values(by='TCGA_gene_name',axis=0)
    pdx_data = pdx_data.pivot_table(columns=['TCGA_gene_name'], index=['sample'], values='counts', aggfunc='sum')

    #Retrieve meta data
    genes_name = np.array(list(pdx_data.columns))
    samples_name = list(pdx_data.index)
    
    #Put the pandas in a numpy array format
    numpy_formatted_data = np.array(pdx_data)

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