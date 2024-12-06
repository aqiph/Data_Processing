#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 17:35:12 2022

@author: guohan
"""

import os, sys
import pandas as pd
import numpy as np
import json
from rdkit import Chem
module_path = '/Users/guohan/Documents/Codes/Data_Processing/utils'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')
from tools import remove_unnamed_columns



### Change other format to .csv ###
def sdf_to_csv(input_file, ID_name='id', SMILES_name='smiles', library_name='', output_file=None, start_id = 0, **kwargs):
    """
    Read input .sdf file, convert to .csv file.
    :param input_file: str, path of the input .sdf file.
    :param ID_name: str, name of ID in .sdf file.
    :param SMILES_name: str, name of SMILES in .sdf file.
    :param library_name: str, library name.
    :param output_file: str, name of the output .csv file.
    :param start_id: int, if ID is not specified in .sdf file, ID is generated starting with start_id.
    """
    # output file
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:
        output_file = basename
    output_file = os.path.join(folder, os.path.splitext(output_file)[0])

    # read mol
    ID_list = []
    SMILES_list = []
    property_name = kwargs.get('property_name', None)
    if property_name is not None:
        property_list = []
    sdf_sup = Chem.SDMolSupplier(input_file)

    for i, mol in enumerate(sdf_sup):
        if mol is None:
            continue

        try:
            ID = mol.GetProp(ID_name).strip()
            ID = library_name + ' ' + ID
        except:
            ID = library_name + ' ' + str(i + start_id)
        ID = ID.strip()

        try:
            SMILES = mol.GetProp(SMILES_name).strip()
        except:
            SMILES = Chem.MolToSmiles(mol)

        if property_name is not None:
            try:
                property = mol.GetProp(property_name)
            except:
                property = np.nan
            property_list.append(property)

        ID_list.append(ID)
        SMILES_list.append(SMILES)

    assert len(ID_list) == len(SMILES_list)
    print('The last id is {}'.format(i + start_id))

    # create DataFrame
    df = pd.DataFrame({'ID': ID_list, 'SMILES': SMILES_list})
    if property_name is not None:
        df[property_name] = property_list

    # write output file
    df = df.reset_index(drop=True)

    print('Number of SMILES:', df.shape[0])
    output_file = '{}_{}.csv'.format(output_file, df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)

    return i + start_id


def json_to_csv(input_file, output_file=None, format='json'):
    """
    Read input .json file, convert to .csv file.
    :param input_file: str, path of the input .json file.
    :param output_file: str, name of the output .csv file.
    :param format: str, 'json' or 'jsonl'.
    """
    # output file
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:
        output_file = basename
    output_file = os.path.join(folder, os.path.splitext(output_file)[0])

    # get data
    if format == 'json':
        with open(input_file, 'r', encoding='utf8') as fr:
            data = json.load(fr)
    elif format == 'jsonl':
        data = []
        with open(input_file, 'r', encoding='utf8') as fr:
            for line in fr:
                data.append(json.loads(line))
    else:
        print('Error: Invalid format for the input file')
        return

    df = pd.DataFrame(data)

    # write output file
    print('Number of SMILES:', df.shape[0])
    output_file = '{}_{}.csv'.format(output_file, df.shape[0])
    # df = remove_unnamed_columns(df)
    df.to_csv(output_file)


### Combine files and split a file ###
def combine_files(input_file_list, columns = None, output_file = None):
    """
    Combine multiple data sets.
    :param input_file_list: list of strs, paths of the input files.
    :param columns: columns in the output_file.
    :param output_file: str or None, pre-defined output file name.
    """
    # output name
    folder, basename = os.path.split(os.path.abspath(input_file_list[0]))
    if output_file is None:
        basename_without_ext  = os.path.splitext(basename)[0]
        output_file = os.path.join(folder, basename_without_ext + '_combine.csv')
    else:
        output_file = os.path.join(folder, output_file)
    
    # read files
    df_list = []
    for file in input_file_list:
        df = pd.read_csv(file)
        df_list.append(df)
    
    # concat DataFrames
    df = pd.concat(df_list, ignore_index = True, sort = False)
    
    # select columns
    if columns is not None:
        df = pd.DataFrame(df, columns = columns)
    
    # write to file
    df = df.reset_index(drop = True)
    print('Number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)


def split_file(input_file, splitting_idx, output_file = None):
    """
    Split input_file into two files.
    :param input_file: str, path of the input file.
    :param splitting_idx: int, the number of rows in the first file.
    :param output_file: str or None, pre-defined output file name.
    """
    # output name
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:        
        basename_without_ext = os.path.splitext(basename)[0]
    else:
        basename_without_ext = os.path.splitext(output_file)[0]
    output_file_1 = os.path.join(folder, basename_without_ext + '_1.csv')
    output_file_2 = os.path.join(folder, basename_without_ext + '_2.csv')
    
    # read and split file
    df = pd.read_csv(input_file, index_col = 0)
    df_1 = df.iloc[:splitting_idx, :]
    df_2 = df.iloc[splitting_idx:, :]
    
    # write to file
    df_1 = df_1.reset_index(drop = True)
    df_2 = df_2.reset_index(drop = True)
    
    print('Number of rows in file 1 is {}'.format(df_1.shape[0]))
    print('Number of rows in file 2 is {}'.format(df_2.shape[0]))
    df_1 = remove_unnamed_columns(df_1)
    df_2 = remove_unnamed_columns(df_2)
    df_1.to_csv(output_file_1)
    df_2.to_csv(output_file_2)


### Get subset ###
def get_subset(input_file, num_cpd, method = 'random', output_file = None):
    """
    Get subset.
    :param input_file: str, path of the input file.
    :param num_cpd: int, the number of compounds.
    :param method: str, the way to select subset: random: randomly selected; an int: the first id of continuous rows.
    :param output_file: str or None, pre-defined output file name.
    """
    # output name
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:
        basename_without_ext  = os.path.splitext(basename)[0]
        output_file = os.path.join(folder, basename_without_ext + '_subset.csv')
    else:
        output_file = os.path.join(folder, output_file)
    
    # read files
    df = pd.read_csv(input_file, index_col = 0)
    
    # get subset
    num_rows = df.shape[0]
    
    if method == 'random':
        df_subset = df.sample(n = num_cpd)
        
    elif isinstance(method, int):      
        df_subset = df.iloc[method : (method + num_cpd)]

    # write to file
    df_subset = df_subset.reset_index(drop = True)
    
    print('Number of rows:', df_subset.shape[0])
    df_subset = remove_unnamed_columns(df_subset)
    df_subset.to_csv(output_file)


### Add SMILES to the query file ###
def add_SMILES(input_file_query, id_column_name_query='ID', input_file_SMILES=None):
    """
    Add SMILES for compounds in input_file_query from SMILES in input_file_SMILES.
    :param input_file_query: str, path of the query compounds.
    :param id_column_name_query: str, name of the ID column in input_file_query
    :param input_file_SMILES: str, path of the input file containing ID and SMILES.
    """
    # files
    df_query = pd.read_csv(input_file_query)
    if id_column_name_query != 'ID':
        df_query.rename(columns={id_column_name_query: 'ID'}, inplace=True)
    if input_file_SMILES is None:
        # input_file_SMILES = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
        input_file_SMILES = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGeneralUse/HTS_forGeneralUse_446664.csv'
    df_SMILES = pd.read_csv(input_file_SMILES)
    df_SMILES = pd.DataFrame(df_SMILES, columns=['ID', 'SMILES'])

    # merge
    df = pd.merge(df_query, df_SMILES, how='left', on=['ID'])

    # write output file
    df = df.reset_index(drop=True)
    print('Number of rows in the file:', df.shape[0])
    df = remove_unnamed_columns(df)
    output_file = f'{os.path.splitext(input_file_query)[0]}_SMILES_{df.shape[0]}.csv'
    df.to_csv(output_file)



if __name__ == '__main__':
    
    ### Change other format to .csv ###
    # input_file = 'tests/sdf_to_csv.sdf'
    # ID_name = 'hit_id'
    # SMILES_name = 'canonical_smiles'
    # library_name = 'test'
    # output_file = 'output.csv'
    # num = sdf_to_csv(input_file, ID_name, SMILES_name, library_name, output_file, start_id=1)
    # print(num)

    # input_file = 'tests/example_json_to_csv.json'
    # output_file = 'example_json_to_csv.csv'
    # json_to_csv(input_file, output_file)


    ### Combine files and split a file ###
    # input_file_list = ['tests/example.csv', 'tests/example2.csv', 'tests/example3.csv']
    # output_file = 'combination.csv'
    # combine_files(input_file_list, output_file = output_file)

    # input_file = 'tests/example.csv'
    # splitting_idx = 10
    # output_file = 'subset.csv'
    # split_file(input_file, splitting_idx, output_file = output_file)


    ### Get subset ###
    # input_file = 'tests/example_noIndex.csv'
    # num_cpd = 10
    # method = 'random'
    # output_file = 'subset.csv'
    # get_subset(input_file, num_cpd, method = method, output_file = output_file)


    ### Add SMILES to the query file ###
    input_file_query = 'tests/test_add_SMILES.csv'
    id_column_name_query = 'Compound_ID'
    input_file_SMILES = None
    add_SMILES(input_file_query, id_column_name_query=id_column_name_query, input_file_SMILES=input_file_SMILES)




    
    