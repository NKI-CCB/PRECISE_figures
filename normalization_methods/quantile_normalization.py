# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 18:26:46 2018

@author: soufiane
"""
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
from normalization_methods.normalization_parameters import NormalizationParameter


def quantile_normalization(count_data, return_instance=False, coef=None):
    
    coef = coef or NormalizationParameter()

    #Custom import
    rpy2.robjects.numpy2ri.activate()
    importr('limma')
    
    #Transform the input count data and feed it to R
    n_samples, n_genes = count_data.shape
    count_cell_lines_R = robjects.r.matrix(count_data.transpose(), nrow=n_genes, ncol=n_samples)
    robjects.r.assign("count_data", count_cell_lines_R)

    X_quantile = robjects.r('''
        quantile_normalized_data <- normalizeQuantiles(count_data)
        ''')

    if coef.parameters is None:
        coef.parameters = True

    if not return_instance:
        #Transpose it back to have it in a scikit-learn format
        return X_quantile.transpose()
    return X_quantile.transpose(), coef
