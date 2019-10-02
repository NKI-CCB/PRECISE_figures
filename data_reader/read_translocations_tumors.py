# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

READ_TRANSLOCATIONS_TUMORS

Reads the translocation data from the TCGA data downloaded from http://www.tumorfusions.org
"""

import pandas as pd
import numpy as np

default_location = '../../data/2018_07_20_biomarkers/translocations/pancanfus.txt'

def read_translocations_tumors(gene_A, gene_B,\
                    tumor_barcodes,\
                    data_location=default_location):
    """
    For a given set of tumor barcode and a gene, finds with a lookup the mutation for this
    particular gene on the TCGA dataset.
    
    INPUT:
        - gene_A (str): first gene of translocation
        - gene_B (str): second gene of translocation
        - tumor_barcodes (list): list of tumor barcodes
        - data_location (str, optional): where data is located
    OUTPUT:
        - indicator list with 1 on tumor barcodes with a translocation
    """

    translocated_genes = [gene_A, gene_B]

    # Read data and filter
    df = pd.read_csv(data_location, sep='\t')
    df = df[np.isin(df.Gene_A, translocated_genes)]
    df = df[np.isin(df.Gene_B, translocated_genes)]

    # Common barcode length
    barcode_length = np.unique([len(e) for e in df['sampleId'].values])
    if barcode_length.shape[0] > 1:
        raise ValueError('File does not the same barcoding length')
    barcode_length = barcode_length[0]
    print(barcode_length)

    # Map translocated tumors
    translocated_barcodes = df['sampleId'].values.astype(str)
    translocated_barcodes = [e.replace('.', '-') for e in translocated_barcodes]
    print(translocated_barcodes)
    translocated_tumors = np.where(np.isin([e[5:5+barcode_length] for e in tumor_barcodes], translocated_barcodes))
    print(translocated_barcodes)
    is_translocated = np.zeros(len(tumor_barcodes))
    is_translocated[translocated_tumors] = 1

    return is_translocated



