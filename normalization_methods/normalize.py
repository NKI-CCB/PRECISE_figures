# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 18:21:03 2018

@author: soufiane

NORMALIZE

Centralize normalization methods. All the code has been picked from different notebooks
used in /results/. The idea is just to get a central code that has not to be called
everytime.

Requirement:
    R 3.14
    numpy
    pandas
"""

import numpy as np

from normalization_methods.TMM_normalization import TMM_normalization
from normalization_methods.DESeq_normalization import DESeq_normalization
from normalization_methods.total_count_normalization import total_count_normalization
from normalization_methods.voom_normalization import voom_normalization
from normalization_methods.normalization_parameters import NormalizationParameter
from normalization_methods.quantile_normalization import quantile_normalization

def normalize_data(count_data, normalization_method, return_instance=False, coef=True):
    
    normalization_method = normalization_method or ''

    if normalization_method.lower() == 'tmm':
        return TMM_normalization(count_data, return_instance, coef)
    
    elif normalization_method.lower() == 'deseq':
        return DESeq_normalization(count_data, return_instance, coef)
    
    elif normalization_method.lower() == 'total_count':
        return total_count_normalization(count_data, return_instance, coef)

    elif normalization_method.lower() == 'upper_quartile':
        upper_quartile_value = [np.percentile(x[np.where(x != 0)],75) for x in count_data]
        mean_count = np.mean(upper_quartile_value)
        normalized_counts =  (1. * count_data.transpose() / upper_quartile_value).transpose() * mean_count

    elif normalization_method.lower() == 'median':
        median_counts = [np.median(x[np.where(x!=0)]) for x in count_data]
        mean_count = np.mean(median_counts)
        normalized_counts = (1. * count_data.transpose() / median_counts).transpose() * mean_count

    elif normalization_method.lower() == 'quantile':
        return quantile_normalization(count_data, return_instance, coef)

    elif normalization_method.lower() == 'voom':
        return voom_normalization(count_data, return_instance, coef)
    
    else:
        print('ERROR: not an available normalization method')
        if return_instance:
            return count_data, NormalizationParameter()
        else:
            return count_data

    if return_instance:
        return normalized_counts, NormalizationParameter()
    return normalized_counts


if __name__ == '__main__':
    from time import time
    import sys
    sys.path.insert(0, '../../src/')
    from data_reader.read_data_one_type import read_data_one_type
    
    X_T, X_CL, gene_names, CL_sample_names = read_data_one_type('Breast', 'BRCA')
    
    start_time = time()
    X_CL_DESeq = normalize_data(X_CL, 'DESeq')
    print('DESeq normalization. Computed in %s s'%(time() - start_time))
    
    start_time = time()
    X_CL_TMM = normalize_data(X_CL, 'TMM')
    print('TMM normalization. Computed in %s s'%(time() - start_time))
    
    start_time = time()
    X_CL_TC = normalize_data(X_CL, 'total_count')
    print('TC normalization. Computed in %s s'%(time() - start_time))
    
    start_time = time()
    X_CL_voom = normalize_data(X_CL, 'voom')
    print('Voom normalization. Computed in %s s'%(time() - start_time))