# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 18:31:41 2018

@author: soufiane
"""

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
from normalization_methods.normalization_parameters import NormalizationParameter


def DESeq_normalization(count_data, return_instance=False, coef=None):
    
    coef = coef or NormalizationParameter()

    #Custom import
    rpy2.robjects.numpy2ri.activate()
    importr('DESeq')
    
    #Transform the input count data and feed it to R
    n_samples, n_genes = count_data.shape
    count_cell_lines_R = robjects.r.matrix(count_data.astype(int).transpose(), nrow=n_genes, ncol=n_samples)
    robjects.r.assign("count_data", count_cell_lines_R)

    # Run DESeq if asked
    if coef.parameters is None:
        robjects.r('''
            size_factor <- estimateSizeFactorsForMatrix(count_data)
            ''')
        coef.parameters = robjects.r["size_factor"]
    else:
        robjects.r.assign('size_factor', coef.parameters)

    #Run DESeq normalization with API
    X_DESeq = robjects.r('''
    DESeq_normalized_data <- sweep(count_data, 2, size_factor,'/')
    ''')

    #Transpose it back
    if not return_instance:
        return X_DESeq.transpose()
    return X_DESeq.transpose(), coef
