#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 17:35:12 2022

@author: guohan
"""

import os
import pandas as pd
from rdkit import Chem



### Helper functions ###

def remove_unnamed_columns(df):
    """
    remove unnamed columns
    """
    unnamed_cols = df.columns.str.contains('Unnamed:')
    unnamed_cols_name = df.columns[unnamed_cols]
    df.drop(unnamed_cols_name, axis=1, inplace=True)
    return df


### Change other format to .csv ###

def sdf_to_csv(input_file, ID_name='id', SMILES_name='smiles', library_name='', output_file=None, start_id = 0):
    """
    read input .sdf file, convert to .csv file
    :param input_file: str, the filename of the input .sdf file
    :param ID_name: str, name of ID in .sdf file
    :param SMILES_name: str, name of SMILES in .sdf file
    :param library_name: str, library name
    :param output_file: str, the filename of the output .csv file
    :param start_id: int, if ID is not specified in .sdf file, ID is generated starting with start_id
    """
    # output name
    folder, basename = os.path.split(os.path.abspath(input_file))
    if output_file is None:
        basename_without_ext = os.path.splitext(basename)[0]
        output_file = os.path.join(folder, basename_without_ext)
    else:
        output_file = os.path.join(folder, os.path.splitext(output_file)[0])

    # read mol
    ID_list = []
    SMILES_list = []
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

        ID_list.append(ID)
        SMILES_list.append(SMILES)

    assert len(ID_list) == len(SMILES_list)
    print('The last id is {}'.format(i + start_id))

    # create DataFrame
    df = pd.DataFrame({'ID': ID_list, 'SMILES': SMILES_list})

    # write output file
    df = df.reset_index(drop=True)

    print('Number of SMILES:', df.shape[0])
    output_file = '{}_{}.csv'.format(output_file, df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)

    return i + start_id


### Combination and splitting ###

def combine_files(input_file_list, columns = None, output_file = None):
    """
    combine multiple data sets
    :param input_file_list: list of strs, the file names of the original input file
    :param columns: columns in the output_file
    :param output_file: str or None, pre-defined output file name
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
    split input_file into two files
    :param input_file: str, file path of the input file
    :param splitting_idx: int, the number of rows in the first file
    :param output_file: str or None, pre-defined output file name
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
    get subset
    :param input_file: str, file path of the input file
    :param num_cpd: int, the number of compounds
    :param method: str, the way to select subset: random: randomly selected; an int: the first id of continuous rows
    :param output_file: str or None, pre-defined output file name
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
    


if __name__ == '__main__':
    
    ### Change other format to .csv ###
    # input_file = 'tests/sdf_to_csv.sdf'
    # ID_name = 'hit_id'
    # SMILES_name = 'canonical_smiles'
    # library_name = 'test'
    # output_file = 'output.csv'
    # num = sdf_to_csv(input_file, ID_name, SMILES_name, library_name, output_file, start_id=1)
    # print(num)


    ### Combination ###
    # input_file_list = ['tests/example.csv', 'tests/example2.csv', 'tests/example3.csv']
    # output_file = 'combination.csv'
    # combine_files(input_file_list, output_file = output_file)
    
    
    ### Splitting ###
    # input_file = 'tests/example.csv'
    # splitting_idx = 10
    # output_file = 'subset.csv'
    # split_file(input_file, splitting_idx, output_file = output_file)
    
    
    ### Get subset ###
    input_file = 'tests/example_noIndex.csv'
    num_cpd = 10
    method = 'random'
    output_file = 'subset.csv'
    get_subset(input_file, num_cpd, method = method, output_file = output_file)




    
    