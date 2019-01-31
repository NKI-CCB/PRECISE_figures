# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 11:27:09 2018

@author: soufiane
"""

from normalization_methods.transform import transform_data
from normalization_methods.normalize import normalize_data
import numpy as np


def feature_engineering(data,\
                        normalization_method = None,\
                        transformation_method = None,\
                        mean_center=False,\
                        std_unit=False,\
                        return_instance=False,\
                        coef=None):
    
    if normalization_method == 'voom':
        voom_transformed_data = normalize_data(data, 'voom', False, coef)
        if transformation_method =='voom':
            voom_transformed_data = np.log(normalize_data(data, 'voom', False, coef))
        else:
            voom_transformed_data = transform_data(voom_transformed_data, transformation_method, False, False, return_instance, coef)
            
        return transform_data(voom_transformed_data, None, mean_center, std_unit, return_instance, coef)
    
    if return_instance:
        normalized_data, coef = normalize_data(data, normalization_method, return_instance, coef)
        return transform_data(normalized_data, transformation_method, mean_center, std_unit, return_instance, coef)
    else:
        normalized_data = normalize_data(data, normalization_method, return_instance, coef)
        return transform_data(normalized_data, transformation_method, mean_center, std_unit, return_instance, coef)
    