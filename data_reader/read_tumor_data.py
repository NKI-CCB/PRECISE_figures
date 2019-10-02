# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

READ_TUMOR_DATA
"""

import xarray as xr
import numpy as np
import pandas as pd



def read_tumor_data(tumor_file,\
                    gene_lookup_file=None,\
                    biospecimen_file=None):
    
    #Read data
    tumors_data = xr.open_dataset(tumor_file)
    tumors_data = tumors_data.sortby('ensemble_gene')
    
    #Transform data
    formatted_data = np.array(tumors_data['counts'])
    
    #Pick gene names
    genes_names = np.array(tumors_data['ensemble_gene'])
    genes_names = np.array([str(name.split('.')[0]) for name in genes_names])
    
    #Remove non-protein coding genes
    if gene_lookup_file is not None:
        genes_lookup_table = pd.read_csv(gene_lookup_file)
        genes_lookup_table = genes_lookup_table[genes_lookup_table.status == 'protein_coding']
        protein_coding_genes = np.array(genes_lookup_table.TCGA_name.values)
    
        #Filter and retrieve data
        protein_coding_genes_index = np.where(np.isin(genes_names, protein_coding_genes))[0]
        formatted_data = formatted_data[:,protein_coding_genes_index]
        genes_names = genes_names[protein_coding_genes_index]
    
    #Reads biospec and return the barcode of each TCGA sample
    if biospecimen_file is not None:
        biospec_data = xr.open_dataset(biospecimen_file)[['barcode']]
        biospec_data = biospec_data.to_dataframe()
        del biospec_data.index.name #for convenience
        aliquot_values = np.array(tumors_data['aliquot']).astype(str)
        barcode_value = [biospec_data[biospec_data.index == al].values[0][0] for al in aliquot_values]
        return formatted_data, genes_names, barcode_value
    
    return formatted_data, genes_names, []