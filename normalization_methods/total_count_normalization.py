# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 18:34:07 2018

@author: soufiane
"""

import numpy as np
from normalization_methods.normalization_parameters import NormalizationParameter

def total_count_normalization(count_data, return_instance=False, coef=None):

	coef = coef or NormalizationParameter()
	if coef.parameters is None:
		coef.parameters = (np.sum(count_data,1))

	if not return_instance:
		return (1.*count_data.transpose() / coef.parameters).transpose() * np.mean(coef.parameters)
	return (1.*count_data.transpose() / coef.parameters).transpose() * np.mean(coef.parameters), coef