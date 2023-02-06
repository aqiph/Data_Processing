#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 16:48:52 2021

@author: guohan
"""

import os
import pandas as pd



def duplicate_analysis(df, output_file, by_column = ['Cleaned_SMILES'], id_column_name = 'ID', smiles_column_name = 'SMILES'):
    """
    analyze the duplicates in df
    :para df: pandas.DataFrame object, input dataframe
    :para output_file: str, output file path without extension for analysis report, folder/filename_without_extension
    :para by_column: list of str, list of column names according to which to group the input file
    :para id_column_name: str, the name of the ID column
    :para smiles_column_name: str, the name of the SMILES column
    """
    # output file
    output_file = '{}_duplicates.txt'.format(os.path.splitext(output_file)[0])
    output = open(output_file, 'w')
    output.write('*****************************\n')
    
    # group df by specified columns   
    gb = df.groupby(by = by_column, dropna = False)
    
    # ID column name and SMILES column name
    columns = df.columns.tolist()
    assert id_column_name in columns, 'Error: {} not found'.format(id_column_name)
    assert smiles_column_name in columns, 'Error: {} not found'.format(smiles_column_name)
    
    # write down info
    num_groups = 0
    
    for labels, df_group in gb:
        num_cmp = df_group.shape[0]
        
        if num_cmp <= 1:
            continue
        
        output.write('Invariants: {}\n'.format(str(labels)))
        output.write('ID      SMILES\n')
        
        IDs = df_group[id_column_name].tolist()
        SMILESs = df_group[smiles_column_name].tolist()
        if 'Value' in columns:
            values = df_group['Value'].tolist()
        else:
            values = [' ' for _ in range(num_cmp)]
        
        for row in range(num_cmp):
            output.write('{}   {}   {}\n'.format(str(IDs[row]), str(SMILESs[row]), str(values[row])))
        
        output.write('*****************************\n')
        
        num_groups += 1
    
    print('The number of duplicating groups is ', num_groups)
    
    output.close()


def remove_duplicates(input_file, by_column = ['Cleaned_SMILES'], keep = 'first', id_column_name = 'ID', smiles_column_name = 'SMILES'):
    """
    clean up smiles with GChem ChEMBL_Structure_Pipeline, add a new column 'Cleaned_SMILES'
    :para input_file: str, the filename of the input file
    :para by_column: list of str, list of column names according to which to group the input file
    :para keep: str, how to keep data if there are duplicates
    :para id_column_name: str, the name of the ID column
    :para smiles_column_name: str, the name of the SMILES column
    """
    # output file path without extension
    output_file, fmt = os.path.splitext(os.path.abspath(input_file))
    
    if fmt in {'.csv'}:
        df = pd.read_csv(input_file, index_col = 0)
    elif fmt in {'.xlsx'}:
        df = pd.read_excel(input_file)
    else:
        print('Error: Invalid format of the input file')
        return
    
    # analysis
    duplicate_analysis(df, output_file, by_column, id_column_name, smiles_column_name)
    
    # remove duplicates
    print('The number of rows before removing duplicates:', df.shape[0])

    df.drop_duplicates(by_column, keep = keep, inplace = True, ignore_index = True)
    df = df.reset_index(drop = True)
    
    # write to file
    print('The number of rows after removing duplicates:', df.shape[0])
    df.to_csv('{}_rmDuplicates.csv'.format(output_file))



if __name__ == '__main__':
    
    input_file = 'tests/example_format_CSP.csv'
    by_column = ['Cleaned_SMILES']
    
    remove_duplicates(input_file, by_column, keep = 'first', id_column_name = 'ID', smiles_column_name = 'SMILES')
    
    
    
    