#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 15:41:39 2022

@author: guohan
"""

import os, sys
path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Data_Processing'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

from cleanup_format import read_input, operate_ID_SMILES_Value, operate_Assay_AssayParameter, operate_Operator_Value_Units, write_output
from cleanup_SMILES import cleanup_smiles
from cleanup_duplicates import remove_duplicates
from add_labels import add_label_column
from analysis import analyze_assay, value_distribution
from util import sdf_to_csv, combine_files, split_file, get_subset


def run_cleanup_format(input_file):
    """
    call cleanup_format to clean up format
    :param input_file: str, the name of the input file
    """
    # file 
    df, folder, output_file = read_input(input_file)    
    filename = os.path.join(folder, output_file)
    
    # ID, SMILES, Value
    id_column_name = 'ID'
    smiles_column_name = 'SMILES'
    value_column_name = 'Value'
    df = operate_ID_SMILES_Value(df, id_column_name, smiles_column_name, value_column_name)
    
    # Assay, Assay_Parameter
    assay_column_name = None
    assay = ''
    assayParameter_column_name = None
    assayParameter = ''
    df = operate_Assay_AssayParameter(df, assay_column_name, assay, assayParameter_column_name, assayParameter)
    
    # Operator, Units
    operator_column_name = 'Value'
    operator = None
    value_column_name = 'Value'
    unit_column_name = None
    unit = ''
    df = operate_Operator_Value_Units(df, operator_column_name, operator, value_column_name, unit_column_name, unit)
    
    write_output(df, filename, use_standard_format=True)


def run_cleanup_SMILES(input_file):
    """
    call cleanup_SMILES to clean up SMILES
    :param input_file: str, the name of the input file
    """
    smiles_column_num = 1
    cleanup_chirality = True
    process_disconnection = True
    process_disconnection_method = 'keep_most_atoms'

    cleanup_smiles(input_file, smiles_column_num, cleanup_chirality, process_disconnection, process_disconnection_method)


def run_cleanup_duplicates(input_file):
    """
    call cleanup_duplicates to clean up duplicates
    :param input_file: str, the name of the input file
    """
    by_column = ['Cleaned_SMILES']
    
    remove_duplicates(input_file, by_column, keep = 'first')


def run_add_labels(input_file):
    """
    call add labels
    :param input_file: str, the name of the input file
    """
    task = 'classification'

    value_column_name = 'Value'
    thresholds = [10.0]
    value_label_correlation = 'negative'
    regression_label_function = 'linear'

    add_label_column(input_file, task, value_column_name=value_column_name, thresholds=thresholds,
                     value_label_correlation=value_label_correlation,
                     regression_label_function=regression_label_function)


def run_analysis(input_file, task = ''):
    """
    call analysis
    :param input_file: str, the name of the input file
    """
    # analyze assay
    if task == 'analyze_assay':
        assay_column_name = 'Assay_Parameter'
        value_column_name = 'Value'
        analyze_assay(input_file, assay_column_name, value_column_name)

    # plot value distribution
    elif task == 'value_distribution':
        value_column_name = 'Value'
        range = [-10, 10]
        value_distribution(input_file, value_column_name, useLog=False, range = range)


def run_util(input_file, task = ''):
    """
    call util
    :param input_file: str, the name of the input file
    """
    if task == 'sdf_to_csv':
        ID_name = 'hit_id'
        SMILES_name = 'canonical_smiles'
        library_name = 'test'
        output_file = 'output.csv'
        num = sdf_to_csv(input_file, ID_name, SMILES_name, library_name, output_file, start_id=1)
        print(num)

    # Combination
    elif task == 'combine_files':
        input_file_list = ['tests/example.csv', 'tests/example2.csv', 'tests/example3.csv']
        output_file = 'combination.csv'
        combine_files(input_file_list, output_file = output_file)

    # Splitting
    elif task == 'split_file':
        splitting_idx = 10
        output_file = 'subset.csv'
        split_file(input_file, splitting_idx, output_file = output_file)

    # Get subset
    elif task == 'get_subset':
        num_cpd = 10
        method = 'random'
        output_file = 'subset.csv'
        get_subset(input_file, num_cpd, method = method, output_file = output_file)


if __name__ == '__main__':

    # input_file = 'tests/example.csv'
    # run_cleanup_format(input_file)

    # input_file = 'tests/example_format.csv'
    # run_cleanup_SMILES(input_file)

    # input_file = 'test/example_format_CSP.csv'
    # run_cleanup_duplicates(input_file)

    # input_file = 'tests/example_format_CSP_rmDuplicates.csv'
    # run_add_labels(input_file)

    # input_file = 'tests/example_format.csv'
    # task = 'value_distribution'
    # run_analysis(input_file, task)

    input_file = 'tests/example.csv'
    task = 'combine_files'
    run_util(input_file, task)









