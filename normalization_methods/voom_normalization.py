# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 18:36:18 2018

@author: soufiane
"""

import numpy as np
from normalization_methods.normalization_parameters import NormalizationParameter

def voom_normalization(count_data, return_instance=False, coef=None):
	
	coef = coef or NormalizationParameter()
	if coef.parameters is None:
		coef.parameters = (np.sum(count_data,1) + 1.)

	if not return_instance:
		return (1.*(count_data.transpose() + 0.5) / coef.parameters).transpose() * 10**6
	else:
		return (1.*(count_data.transpose() + 0.5) / coef.parameters).transpose() * 10**6, coef