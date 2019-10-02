# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

Harmonize the features between the target and the source data so that:
- same feature space is considered between the source and the target.
- features are odered in the same way, avoiding permutation issue.
"""

import numpy as np
import pandas as pd


def harmonize_feature_naming(target_data, 
                            source_data,
                            target_gene_names,
                            source_gene_names,
                            remove_mytochondria=False,
                            gene_lookup_file=None):
    
    #Find common genes
    common_genes = np.intersect1d(target_gene_names, source_gene_names)
    
    #Find common gene location
    target_common_gene_index = np.where(np.isin(target_gene_names, common_genes))[0]
    source_common_gene_index = np.where(np.isin(source_gene_names, common_genes))[0]
    
    #Stack data
    target_data = target_data[:,target_common_gene_index]
    source_data = source_data[:,source_common_gene_index]

    if remove_mytochondria and gene_lookup_file is not None:
        #Load table
        genes_lookup_table = pd.read_csv(gene_lookup_file, delimiter=',')
        print(genes_lookup_table.shape)
        genes_lookup_table = genes_lookup_table.drop_duplicates(subset=['ENSEMBL'], keep=False)
        print(genes_lookup_table.shape)

        #Transform gene names
        df_genes_name = pd.DataFrame(common_genes, columns=['ENSEMBL'])
        df_genes_name = df_genes_name.merge(genes_lookup_table, on='ENSEMBL', how='left')
        chromosome_name = df_genes_name['chromosome_name'].values

        #Filter data on mitochondria
        accepted_chromosomes = np.array(list(range(1,23)) + ['X','Y'])
        filtered_genes_index = np.where(np.isin(chromosome_name, accepted_chromosomes))[0]
        target_data = target_data[:,filtered_genes_index]
        source_data = source_data[:,filtered_genes_index]
        common_genes = common_genes[filtered_genes_index]

    
    return target_data, source_data, common_genes