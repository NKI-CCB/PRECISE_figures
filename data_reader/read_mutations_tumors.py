# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

READ_MUTATIONS_TUMORS

Reads the mutations data from the TCGA data.
Data should be in the format of the cBIO portal. Tumor barcodes are shortened to
comply with the input data.
"""

import pandas as pd
import numpy as np

default_location = '../../data/2018_07_20_biomarkers/tcga_skin/BRAF_mutation_status.csv'
default_detail_location = '../../data/2018_07_20_biomarkers/tcga_skin/BRAF_mutation_detailed.csv'

def read_mutations_tumors(gene_name,\
                    tumor_barcodes,\
                    mutation_data_location=default_location,\
                    mutation_detailed_data_location=default_detail_location):
    """
    For a given set of tumor barcode and a gene, finds with a lookup the mutation for this
    particular gene on the TCGA dataset.
    
    INPUT:
        - gene_name (str): name of the gene considered
        - tumor_barcodes (list): list of tumor barcodes
        - mutation_data_location (str, optional, default to skin for tcga): where data is
        located
    OUTPUT:
        - list of mutation status with -1 if tumor has not been found
    """

    # Read mutation data
    mutation_data = pd.read_csv(mutation_data_location, sep ='\t').transpose()
    mutation_data = mutation_data[[0,2]]
    mutation_data.columns = ['CLINICAL', 'MUTATIONS']
    mutation_data = mutation_data.iloc[2:]
    mutation_data = mutation_data[mutation_data.CLINICAL == 'Yes']

    # Filter the samples according to their mutation status
    samples = mutation_data.index
    non_mutated_samples = mutation_data[pd.isna(mutation_data.MUTATIONS)].index.astype(np.str)
    mutated_samples = mutation_data[~pd.isna(mutation_data.MUTATIONS)].index.astype(np.str)

    # Compute mutation status
    mutation_status = - np.ones(len(tumor_barcodes))

    length_barcode = min([len(e) for e in samples])
    truncated_tumor_names = [e[:length_barcode] for e in tumor_barcodes]
    non_mutated_samples = [e[:length_barcode] for e in non_mutated_samples]
    mutated_samples = [e[:length_barcode] for e in mutated_samples]

    mutation_status[np.where(np.isin(truncated_tumor_names, non_mutated_samples))] = 0
    mutation_status[np.where(np.isin(truncated_tumor_names, mutated_samples))] = 1

    if mutation_detailed_data_location is None:
        return mutation_status

    # Read more details
    detail_data = pd.read_csv(mutation_detailed_data_location, sep='\t')
    detail_data = detail_data[['Sample ID', 'Protein Change']]

    mutation_status = mutation_status.astype(np.str)
    for tumor_name, protein_change in detail_data.values:
        mutation_status[np.isin(truncated_tumor_names, [tumor_name[:length_barcode]])] = protein_change

    return mutation_status