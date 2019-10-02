# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

READ DATA

General function to read data from. Given a type (RNA-Seq or FPKM), a source type
(cell-line, PDX) and a target type (PDX or tumor), retrieve the data from the downloaded
files and convert them into a (samples, genes) format. Further processing is done so that
source and target have the same genes, in the same order.

Source and target sample names (barcodes in the case of tumors) are also retrieved.
"""

from data_reader.read_tumor_data import read_tumor_data
from data_reader.read_cell_line_data import read_cell_line_data, read_cell_line_data_fpkm, read_cell_line_data_cell_passport
from data_reader.read_pdx_data import read_pdx_data
from data_reader.harmonize_feature_naming import harmonize_feature_naming

# For PRECISE-like files, i.e. in RDS
precise_cl_folder = '../data/cell_line/'
# For integrating the new dataset from Cell Passport
cell_passport_cl_folder = '../data/cell_line/'
# For tumors
tumor_folder = '../data/tumor/'
# For PDX
pdx_folder = '../data/pdx/'
# For lookup data
lookup_folder = '../data/lookup/'

def read_data(source_type,
            target_type,
            data_type,
            source_tissue=None,
            target_tissue=None,
            remove_mytochondria=False):
    """
    Read data located in different files with one source and one target type and 
    homogenize the features to keep the overlapping ones in the same order.
    Sample names are also retrieved.
    
    INPUT:
        - source_type (str): type of the source, i.e. cell_line or pdx.
        - target_type (str): type of the target, i.e. pdx or tumor.
        - data_type (str): type of data, i.e. fpkm, count or count_passport.
        - source_tissue (str, optional, default to None): tissue of origin for the source
        - target_tissue (str, optional, default to None): tissue of origin for the target
        - remove_mytochondria (bool, optional, default to False): whether mythocondrial 
        genes should be removed.
    OUTPUT:
        - target_data (np.ndarray): array with target data in the form (n_samples, n_genes).
        - source_data (np.ndarray): array with source data in the form (n_samples, n_genes).
        - gene_names (np.ndarray): array with the names of the genes, in the same order than 
        the features in target_data and source_data.
        - source_samples (np.ndarray): array with the source sample names, in the same order
        than in source_data.
        - target_samples (np.ndarray): array with the target sample names, in the same order
        than in target_data.
    """

    # Read source data
    source_data, source_gene_names, source_samples = read_one_data_source(source_type,
                                                                          data_type,
                                                                          source_tissue)

    # Read target data
    target_data, target_gene_names, target_samples = read_one_data_source(target_type,
                                                                          data_type,
                                                                          target_tissue)

    # Homegenize to get same genes, i.e. features
    gene_lookup_file = '%sgene_status.csv'%(lookup_folder)
    target_data, source_data, gene_names = harmonize_feature_naming(target_data,
                                                                  source_data,
                                                                  target_gene_names,
                                                                  source_gene_names,
                                                                  remove_mytochondria,
                                                                  gene_lookup_file)

    return target_data, source_data, gene_names, source_samples, target_samples


def read_one_data_source(model_type, data_type, tissue):
    """
    Read data corresponding to one source and one data type.
    
    INPUT:
        - model_type (str): type of the source, i.e. cell_line pdx or tumor.
        - data_type (str): type of data, i.e. fpkm, count or count_passport.
        - tissue (str, optional, default to None): tissue of origin
    OUTPUT:
        - data (np.ndarray): array with target data in the form (n_samples, n_genes).
        - gene_names (np.ndarray): array with the names of the genes, in the same order than 
        the features in data.
        - samples (np.ndarray): array with the target sample names, in the same order
        than in data.
    """
    gene_lookup_file = '%sgene_status.csv'%(lookup_folder)
    biospec_file = '%sbiospec_%s'%(tumor_folder, tissue)

    if model_type.lower() == 'tumor':
        tumor_file = '%s%s_%s_netcdf'%(tumor_folder,
                                      'count' if data_type == 'count_passport' else data_type,
                                      tissue)

        return read_tumor_data(tumor_file,
                              gene_lookup_file,
                              biospec_file)

    elif model_type.lower() == 'cell_line':
        cell_line_lookup_file = '%scancer_type.csv'%(precise_cl_folder)

        if data_type.lower() == 'fpkm':
            cell_line_file = '%srnaseq_fpkm_protein_coding.csv'%(precise_cl_folder)
            return read_cell_line_data_fpkm(cell_line_file,
                                          gene_lookup_file,
                                          cell_line_lookup_file,
                                          tissue)

        elif data_type.lower() == 'count':
            cell_line_file = '%srnaseq_readcounts_TCGA.RDS'%(precise_cl_folder)
            return read_cell_line_data(cell_line_file,
                                      gene_lookup_file,
                                      cell_line_lookup_file,
                                      tissue)

        elif data_type.lower() == 'count_passport':
            cell_line_file = '%srnaseq_latest.csv'%(cell_passport_cl_folder)
            gene_lookup_file = '%sgene_identifiers_latest.csv'%(cell_passport_cl_folder)
            cell_line_lookup_file = '%smodel_list_latest.csv'%(cell_passport_cl_folder)
            return read_cell_line_data_cell_passport(cell_line_file,
                                                  gene_lookup_file,
                                                  cell_line_lookup_file,
                                                  tissue)

    elif model_type.lower() == 'pdx':
        if data_type.lower() != 'fpkm':
            raise ValueError('FPKM not available for PDX.')

        pdx_file = '%spdx_%s_TCGA_index_fpkm.csv'%(pdx_folder, tissue)
        return read_pdx_data(pdx_file,
                            gene_lookup_file,
                            None)