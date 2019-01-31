"""
Created on Mon Mar  5 18:21:03 2018

@author: soufiane

NORMALIZE

Centralize Transformation methods. All the code has been picked from different notebooks
used in /results/. The idea is just to get a central code that has not to be called
everytime.

Requirement:
    numpy
    pandas
"""

import numpy as np
from sklearn.preprocessing import StandardScaler, quantile_transform
from normalization_methods.normalization_parameters import NormalizationParameter


def transform_data(count_data, transformation_method,\
                   mean_center=False,\
                   std_unit=False,\
                   return_instance=False,\
                   coef=None):
    
    transformation_method = transformation_method or ''
    
    if transformation_method.lower() == 'log':
        transformed_data = np.log(count_data + 1)
    
    elif transformation_method.lower() == 'voom':
        transformed_data = np.log(count_data)
    
    elif transformation_method.lower() == 'anscombe':
        transformed_data = np.sqrt(count_data + 3./8.)
        
    elif transformation_method.lower() == 'quantile':
        transformed_data = quantile_transform(count_data, output_distribution='normal')
    
    else:
        print('ERROR: not an available transformation method')
        transformed_data = count_data
        
        
    # Create NormalizationParameter instance for further processing.
    coef = coef or NormalizationParameter()
    if coef.scaler is None:
        coef.scaler = StandardScaler(with_mean=mean_center, with_std=std_unit)
        coef.scaler.fit(transformed_data)

    if return_instance:
        return coef.scaler.transform(transformed_data), coef
    else:
        return coef.scaler.transform(transformed_data)
    



if __name__ == '__main__':
    from time import time
    import sys
    sys.path.insert(0, '../../src/')
    from data_reader.read_data_one_type import read_data_one_type
    
    X_T, X_CL, gene_names, CL_sample_names = read_data_one_type('Breast', 'BRCA')
    
    start_time = time()
    X_CL_log = transform_data(X_CL, 'log')
    print('LOG transformation. Computed in %s s'%(time() - start_time))
    
    start_time = time()
    from normalize import normalize_data
    X_CL_voom = transform_data(normalize_data(X_CL,'voom'), 'voom')
    print('VOOM transformation. Computed in %s s'%(time() - start_time))
    
    start_time = time()
    X_CL_log_z_score = transform_data(X_CL, 'log', True, True)
    print('LOG_Z-SCORE transformation. Computed in %s s'%(time() - start_time))
    
    start_time = time()
    from normalize import normalize_data
    X_CL_voom_z_score = transform_data(normalize_data(X_CL,'voom'), 'voom', True, True)
    print('VOOM_Z-SCORE transformation. Computed in %s s'%(time() - start_time))
    
    start_time = time()
    from normalize import normalize_data
    X_CL_anscombe = transform_data(X_CL, 'anscombe')
    print('ANSCOMBE transformation. Computed in %s s'%(time() - start_time))