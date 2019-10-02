# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

READ_CNA_TUMORS

Reads the Copy Number Alterations and associate it to the tumor data using
barcode identification.
"""

import pandas as pd
import numpy as np

default_location = '../../data/2018_03_29_cna/tcga_breast/data.txt'

def read_cna_tumors(gene_name,\
                    tumor_barcodes,\
                    cna_data_location=default_location):
    """
    For a given set of tumor barcode and a gene, finds with a lookup the cna for this
    particular gene on the TCGA dataset.
    The code takes as input a csv table separated with '\t' with genes in column and 
    barcodes in columns
    
    INPUT:
        - gene_name (str): name of the gene considered
        - tumor_barcodes (array): array of tumor barcodes
        - cna_data_location (str, optional, default to breast for tcga): where data is
        located
    OUTPUT:
        - list of the cna with nan if tumor has not been found
    """
    
    #Read CNA
    cna_df = pd.read_csv(cna_data_location, sep='\t')
    
    #Select the relevant row and transpose for later filtering
    cna_df_filtered = cna_df[cna_df['Hugo_Symbol'] == gene_name]
    cna_df_filtered = cna_df_filtered.transpose()
    
    #Shorten the barcode for later usage
    length_barcode = np.min([len(e) for e in cna_df_filtered.index.values][2:])
    shortened_barcode = [barcode[:length_barcode] for barcode in tumor_barcodes]
    cna_df_filtered.index = [e[:length_barcode] for e in cna_df_filtered.index.values]
    
    #Retrieve cna
    cna_tumors = [cna_df_filtered[cna_df_filtered.index == barcode] for barcode in shortened_barcode]
    cna_tumors = [cna.values[0,0] if cna.shape[0] >= 1 else np.nan for cna in cna_tumors]
    
    return np.array(cna_tumors)



if __name__ == '__main__':
    from read_data_one_type import read_data_one_type
    
    #X_tumors, X_cell_lines, gene_names, source_sample_names, barcodes = read_data_one_type('Lung',\
    #                                                        None, True)
    
    cna_tumors = read_cna_tumors('EGFR', barcodes, '../../data/2018_03_29_cna/tcga_lung/data.txt')