#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 6 15:21:00 2024

@author: guohan
"""

import os, sys
import pandas as pd
import numpy as np
module_path = '/Users/guohan/Documents/Codes/Data_Processing/utils'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')
from tools import remove_unnamed_columns



### Add SMILES to the query file ###
def add_SMILES(input_file_query, id_column_name_query='ID'):
    """
    Add SMILES for compounds in input_file_query from SMILES in input_file_SMILES.
    :param input_file_query: str, path of the query compounds.
    :param id_column_name_query: str, name of the ID column in input_file_query
    """
    # files
    df_query = pd.read_csv(input_file_query)
    if id_column_name_query != 'ID':
        df_query.rename(columns={id_column_name_query: 'ID'}, inplace=True)
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


### Add library ID to the query file ###
def get_library_match(input_file, input_library, id_column_name='ID'):
    """
    Match compounds in input_file with compounds in input_library.
    :param input_file: str, path of the input file.
    :param input_library: str, path of the input library.
    :param id_column_name: str, name of the ID column in the input_file.
    """
    df = pd.read_csv(input_file)
    df_lib = pd.read_csv(input_library)
    lib_ID = str(os.path.split(input_library)[1]).split('_')[1]
    COLUMNS = df.columns.tolist() + ['Library_ID']

    df['ID_'] = df[id_column_name].apply(lambda id: str(id))
    df_lib['ID_'] = df_lib['ID'].apply(lambda id: str(id))
    df_lib = pd.DataFrame(df_lib, columns=['ID_'])

    df_int = pd.merge(df, df_lib, how='inner', on=['ID_'])
    df_int['Library_ID'] = lib_ID
    df_int = pd.DataFrame(df_int, columns=COLUMNS)
    print(f'Number of cmps in {lib_ID} is {df_int.shape[0]}')

    return df_int


def add_libID(input_file_query, id_column_name_query='ID'):
    """
    Add library ID for compounds in input_file_query.
    """
    # output file path without extension
    output_file = os.path.splitext(os.path.abspath(input_file_query))[0]

    folder = '/Users/guohan/Documents/Projects/Datasets/HTS/'
    Library = ['Disease-Based-Library', 'Fragment-Library', 'High-Throughput-Library',
               'Nature-Product-Library', 'Repurposing-Library']
    Sublibrary = [['Anti-Virus', 'TBAC-Alliance', 'TBAC-Calibr'],
                  ['Fragment-Enamine', 'Fragment-MayBridge'],
                  ['AnalytiCon-NATx', 'ChemBridgeCL20K', 'ChemBridgeCL90K', 'ChemBridgeEXP50K', 'ChemDiv-SMART-TM',
                   'Expresspick-Pfizer', 'LifeChemical-30K', 'MayBridge-50K', 'Selleck-PFZ'],
                  ['AnalytiCon-MEGx-NP', 'MetaSci-HML-NP', 'Pharmacodia-NP', 'TargetMol-NP'],
                  ['FDA', 'Prc&Cli', 'ReFrame']]

    # loop over sublibraries
    dfs = []
    for i, library in enumerate(Library):
        for j, sublibrary in enumerate(Sublibrary[i]):
            # library file
            files = os.listdir(f'{folder}/{library}/{sublibrary}_20230317/final')
            filtered_files = [file for file in files if file.startswith(f'HTS_{sublibrary}_forGeneralUse_')]
            input_library = f'{folder}/{library}/{sublibrary}_20230317/final/{filtered_files[0]}'

            df_lib_match = get_library_match(input_file_query, input_library, id_column_name=id_column_name_query)
            dfs.append(df_lib_match)

    # concat DataFrames
    df_match = pd.concat(dfs, ignore_index=True, sort=False)
    df_match['ID_'] = df_match[id_column_name_query].apply(lambda id: str(id))
    df_match = pd.DataFrame(df_match, columns=['ID_', 'Library_ID'])
    df = pd.read_csv(input_file_query)
    COLUMNS = df.columns.tolist() + ['Library_ID']
    df['ID_'] = df[id_column_name_query].apply(lambda id: str(id))
    df = pd.merge(df, df_match, how='left', on=['ID_'])
    df = pd.DataFrame(df, columns=COLUMNS)

    # write output file
    df = df.reset_index(drop=True)
    print('Number of compounds:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(f'{output_file}_{df.shape[0]}.csv')



if __name__ == '__main__':
    ### Add SMILES to the query file ###
    # input_file_query = 'tests/test_add_SMILES.csv'
    # id_column_name_query = 'Compound_ID'
    # add_SMILES(input_file_query, id_column_name_query=id_column_name_query)


    ### Add library ID to the query file ###
    input_file_query = 'tests/test_add_libID.csv'
    add_libID(input_file_query, id_column_name_query='ID')

