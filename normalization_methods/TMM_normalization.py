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


def TMM_normalization(count_data, return_instance=False, coef=None):
    
    coef = coef or NormalizationParameter()

    #Custom import
    rpy2.robjects.numpy2ri.activate()
    importr('edgeR')
    
    #Transform the input count data and feed it to R
    n_samples, n_genes = count_data.shape
    count_cell_lines_R = robjects.r.matrix(count_data.transpose(), nrow=n_genes, ncol=n_samples)
    robjects.r.assign("count_data", count_cell_lines_R)

    # Recompute the coefficients is asked, use precomputed otherwise
    robjects.r('''
        D <- DGEList(counts=count_data)
        ''')
    if coef.parameters is None:
        robjects.r('''
            #TMM normalization
            Dnorm <- calcNormFactors(D)
            ''')
        coef.parameters = robjects.r["Dnorm"]
    else:
        robjects.r.assign('Dnorm', coef.parameters)
    
    #Run TMM normalization with API
    X_TMM = robjects.r('''
        rellibsize <- colSums(D$counts)/exp(mean(log(colSums(D$counts))))
        nf = Dnorm$samples[,3]*rellibsize
        TMM_normalized_data = round(sweep(D$counts, 2, nf, "/"))
        ''')

    if not return_instance:
        #Transpose it back to have it in a scikit-learn format
        return X_TMM.transpose()
    return X_TMM.transpose(), coef
