"""
@author: soufiane

NORMALIZATION_PARAMETERS

Keeps a trace of the normalization parameters in an attempt to use them afterwards.

Requirement:
    numpy
    pandas
"""

import numpy as np
from sklearn.preprocessing import StandardScaler

class NormalizationParameter():

	def __init__(self, parameters=None, scaler=None):
		self.parameters = parameters
		self.scaler = scaler
		pass