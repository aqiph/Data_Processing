#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 10:03:58 2022

@author: guohan
"""

import os
import pandas as pd
import numpy as np



def read_input(input_file):
    """
    read input file
    :param input_file: str, the filename of the original input file
    """
    folder, basename = os.path.split(os.path.abspath(input_file))
    output_file, fmt = os.path.splitext(basename) # basename without extension
    
    if fmt in {'.csv'}:
        df = pd.read_csv(input_file)
    elif fmt in {'.xlsx'}:
        df = pd.read_excel(input_file)
    else:
        print('Error: Invalid format for the input file')
        return
    
    return df, folder, output_file


def operate_ID_SMILES_Value(df, id_column_name = None, smiles_column_name = 'SMILES', value_column_name = None, dropna_column_names = ['SMILES']):
    """
    rename or add ID, SMILES and Value columns
    :param df: pandas.DataFrame object, input dataframe
    :param id_column_name: str or None, the name of the ID column; if None, generate ints as new IDs.
    :param smiles_column_name: str, the name of the SMILES column
    :param value_column_name: str or None, the name of the Value column; if None, generate an empty column
    :param dropna_column_names: list of str, the column names depend on which to drop nan value
    """    
    columns = df.columns.tolist()
    
    # change ID column name or add ID column
    if id_column_name is None or id_column_name == 'None':
        ids = np.arange(0, df.shape[0], 1, dtype = int)
        df['ID'] = ids
    else:
        assert id_column_name in columns, 'Error: Invalid input ID column name'
        df.rename(columns = {id_column_name:'ID'}, inplace = True)
    df['ID'] = df['ID'].apply(lambda v: to_str(v))
        
    # change SMILES column name 
    assert smiles_column_name in columns, 'Error: Invalid input SMILES column name'
    df.rename(columns = {smiles_column_name:'SMILES'}, inplace = True)
    df['SMILES'] = df['SMILES'].apply(lambda v: to_str(v))
    
    # change Value column name or add Value column
    if value_column_name is None or value_column_name == 'None':
        df['Value'] = np.nan
    else:
        assert value_column_name in columns, 'Error: Invalid input Value column name'
        df.rename(columns = {value_column_name:'Value'}, inplace = True)
    
    # drop nan
    df.dropna(subset = dropna_column_names, how = 'any', inplace = True)
    
    # write to file
    df = df.reset_index(drop = True)    
    print('Clean up ID, SMILES and Value columns, number of rows:', df.shape[0])
    
    return df


def operate_Assay_AssayParameter(df, assay_column_name = None, assay = '', assayParameter_column_name = None, assayParameter = ''):
    """
    rename or add Assay and Assay_Parameter columns
    :param df: pandas.DataFrame object, input dataframe
    :param assay_column_name: str or None, the name of the Assay column; if None, generate Assay column according to 'assay'
    :param assay: str, assay name in the Assay column
    :param assayParameter_column_name: str or None, the name of the Assay_Parameter column; if None, generate Assay_Parameter column according to 'assayParameter'
    :param assayParameter: str, assay parameter name in the Assay_Parameter column
    """
    columns = df.columns.tolist()
    
    # change Assay column name or add Assay column according to 'assay'
    if assay_column_name is None or assay_column_name == 'None':
        df['Assay'] = assay
    else:
        assert assay_column_name in columns, 'Error: Invalid input Assay column name'
        df.rename(columns = {assay_column_name:'Assay'}, inplace = True)
    df['Assay'] = df['Assay'].apply(lambda v: to_str(v))
        
    # change Assay_Parameter column name or add Assay_Parameter column according to 'assayParameter'
    if assayParameter_column_name is None or assayParameter_column_name == 'None':
        df['Assay_Parameter'] = assayParameter
    else:
        assert assayParameter_column_name in columns, 'Error: Invalid input Assay_Parameter column name'
        df.rename(columns = {assayParameter_column_name:'Assay_Parameter'}, inplace = True)
    df['Assay_Parameter'] = df['Assay_Parameter'].apply(lambda v: to_str(v))
    
    # write to file
    df = df.reset_index(drop = True)
    print('Clean up Assay and Assay_Parameter columns, number of rows:', df.shape[0])
    
    return df


def operate_Operator_Value_Units(df, operator_column_name = None, operator = '=', value_column_name = 'Value', unit_column_name = None, unit = ''):
    """
    rename or add Operator and Units columns
    :param df: pandas.DataFrame object, input dataframe
    :param operator_column_name: str or None, the name of the Operator column; if None, generate Operator column according to operator
    :param operator: str, operator name in the Operator column
    :param value_column_name: str, the name of the Value column
    :param unit_column_name: str or None, the name of the Units column; if None, generate Units column according to unit
    :param unit: str, unit name in the Units column
    """
    columns = df.columns.tolist()
    
    # change Operator column name, or add Operator column according to 'operator' or splitting from Value column
    if operator_column_name is None or operator_column_name == 'None':
        df['Operator'] = operator
    else:
        assert operator_column_name in columns, 'Error: Invalid input Operator column name'
        if operator_column_name == value_column_name: # operator and value in the same cell
            df['Value'] = df['Value'].apply(lambda x: split_operator_value(x))
            df[['Operator', 'Value']] = pd.DataFrame(df['Value'].values.tolist())
        else:
            df.rename(columns = {operator_column_name:'Operator'}, inplace = True)
    df['Operator'] = df['Operator'].apply(lambda v: to_str(v))
    
    # change Units column name, or add Units column according to 'unit'
    if unit_column_name is None or unit_column_name == 'None':
        df['Units'] = unit
    else:
        assert unit_column_name in columns, 'Error: Invalid input Units column name'
        df.rename(columns = {unit_column_name:'Units'}, inplace = True)
    df['Units'] = df['Units'].apply(lambda v: to_str(v))
    
    # write to file
    df = df.reset_index(drop = True)
    print('Clean up Operator and Units columns, number of rows:', df.shape[0])
    
    return df


def to_str(value):
    """
    change value to a string if it is not np.nan or empty string, else keep it as np.nan
    """
    s = str(value).strip()
    if s in {'', 'nan'}:
        return np.nan
    return s


def split_operator_value(raw):
    """
    Helper function for cleanup_format
    split operator and value
    :param raw: str, original value
    :return: list, [operator, value]
    """
    raw = str(raw).strip()
    
    if raw[:2] in {'<=', '>='}:
        operator, number = raw[:2], raw[2:]
    elif raw[:1] in {'<', '>', '='}:
        operator, number = raw[:1], raw[1:]
    elif raw[0].isdigit():
        operator, number = '=', raw        
    else:
        return ['=', raw]
    
    number = float(number.replace(',', '').replace(' ', ''))
    return operator, number


def remove_unnamed_columns(df):
    """
    remove unnamed columns
    """
    unnamed_cols = df.columns.str.contains('Unnamed:')
    unnamed_cols_name = df.columns[unnamed_cols]
    df.drop(unnamed_cols_name, axis=1, inplace=True)
    return df


def write_output(df, output_file, use_standard_format=False):
    """
    write to output file
    :param df: pandas.DataFrame object, input dataframe
    :param output_file: str, output file path without extension
    :param use_standard_format: bool, whether to change columns to standard format or not
    """
    df = df.reset_index(drop = True)
    print('Number of rows:', df.shape[0])
    if use_standard_format:
        df = pd.DataFrame(df, columns = ['ID', 'SMILES', 'Assay', 'Assay_Parameter', 'Operator', 'Value', 'Units'])
    else:
        df = remove_unnamed_columns(df)
    df.to_csv('{}_format.csv'.format(output_file))



if __name__ == '__main__':
    input_file = 'tests/example.csv'
    
    df, folder, output_file = read_input(input_file)    
    filename = os.path.join(folder, output_file)
    
    id_column_name = 'Compound ID'
    smiles_column_name = 'SMILES'
    value_column_name = 'EC50 (nM)'
    dropna_column_names = ['SMILES']
    df = operate_ID_SMILES_Value(df, id_column_name, smiles_column_name, value_column_name, dropna_column_names)
    
    assay_column_name = None
    assay = 'Activity'
    assayParameter_column_name = None
    assayParameter = 'EC50'
    df = operate_Assay_AssayParameter(df, assay_column_name, assay, assayParameter_column_name, assayParameter)
    
    operator_column_name = 'Value'
    operator = None
    value_column_name = 'Value'
    unit_column_name = None
    unit = 'nM'
    df = operate_Operator_Value_Units(df, operator_column_name, operator, value_column_name, unit_column_name, unit)
    
    write_output(df, filename, use_standard_format=True)
    
