# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

READ_DRUG_RESPONSE.PY

Reads the drug response from the given files and returns the IC50 for cell lines.
Comprehensive reader with lookup function.
"""

import numpy as np
import pandas as pd


def read_drug_response(drug_id,
                       source_data,
                       source_sample_names,
                       type_data='rnaseq'):
    """
    Reads data for drug response.
    
    INPUT:
        - drug_id (int): ID of the drug
        - source_data (np.ndarray): data corresponding to the source, 
        in the form (n_samples, n_genes).
        - source_sample_names (np.ndarray): array of sample names, in the same order 
        as source_data rows.
        - type_data (optional, default to count): type of data considered, count for cell lines.
    OUTPUT:
        - data (2d array) : data filtered to give only the filtered source.
        - drug_response (array): drug response in the same order as data.
        - sample_names : names of the source samples with known drug response. Same ordering
        as in data and drug_response.
        - drug_name (str) : name of the drug.
    """
    
    if type_data == 'count':
        #Files of data
        cell_lines_data_folder = './data/cell_line/'
        cell_lines_drug_response_file = 'ic50.csv'
        drugs_specifications_file = 'drugs_specifications.csv'
        
        return read_drug_response_cell_lines(drug_id,
                                            cell_lines_data_folder + cell_lines_drug_response_file,
                                            source_data,
                                            source_sample_names,
                                            cell_lines_data_folder + drugs_specifications_file)

    else:
        raise ValueError('%s if not an available type. Should be \'count\' or \'fpkm\'')
    
    
    
def read_drug_response_cell_lines(drug_id,
                       drug_response_file,
                       cell_lines_data,
                       cell_lines_sample_names,
                       drug_specification_file):
    
    """
    Reads data for drug response and filter cell lines that have an output
    for the drug response problem.
    
    INPUT:
        - drug_id (int): ID of the drug
        - drug_response_file (str): where the drug responses are stored
        - drug_specification_file (str): where specification for drugs are stored
        - cell_lines_data (2d array): cell lines data previously found
        - cell_lines_sample_names (array): names of samples
    OUTPUT:
        - cell_lines_data (2d array) : data filtered to give only the filtered cell lines
        - drug_response (array): drug response in the same order as cell_lines_data
        - cell_lines_sample_names : names of the cell lines with known drug response
        - drug_name (str) : name of the drug
    """
        
    #Read drug response and process data
    drug_response = pd.read_csv(drug_response_file)
    drug_response.columns = ['sample_name'] + list(drug_response.columns[1:])
    drug_response = drug_response[['sample_name', str(drug_id)]]
    drug_response = drug_response.dropna()
    
    #Take common sample
    common_samples_names = np.intersect1d(cell_lines_sample_names, drug_response.sample_name.values)
    common_samples_names_index = np.where(np.isin(cell_lines_sample_names, common_samples_names))[0]
    cell_lines_sample_names = cell_lines_sample_names[common_samples_names_index]
    
    #Retrieve values
    cell_lines_data = cell_lines_data[common_samples_names_index,:]
    drug_response = np.array([drug_response[drug_response.sample_name == e][str(drug_id)].values[0] \
                         for e in cell_lines_sample_names])
    
    #Read name
    drug_name = pd.read_csv(drug_specification_file, sep=',')
    drug_name = drug_name[drug_name.Identifier == drug_id]['Name'].values[0]
    
    return cell_lines_data, drug_response, cell_lines_sample_names, drug_name