#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 10:05:14 2021

@author: guohan
"""

import os, warnings
import pandas as pd
import numpy as np
from rdkit.Chem.SaltRemover import SaltRemover
from chembl_structure_pipeline import *
from chembl_structure_pipeline.checker import *
from rdkit import RDLogger



def cleanup_single_smiles_by_CSP(smiles, cleanup_chirality = False):
    """
    Clean up a single smiles with chembl_structure_pipeline.
    :param smiles: str, smiles.
    :param cleanup_chirality: bool, whether to remove chirality or not.
    :return: cleaned smiles by chembl_structure_pipeline, flag to indicate if this smiles is valid.
    """
    # disable rdkit logging
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)

    # True if the input smiles can be properly converted, False if there is an error
    flag = True
    
    try:
        # get mol object
        mol = Chem.MolFromSmiles(smiles)
        # standardize mol
        mol_std = standardize_mol(mol)
        # get parent
        mol_par = get_parent_mol(mol_std)[0]
        # canonicalize SMILES, remove chirality if required
        if cleanup_chirality:
            smiles = Chem.MolToSmiles(mol_par, isomericSmiles = False)
        else:
            smiles = Chem.MolToSmiles(mol_par)        
        mol = Chem.MolFromSmiles(smiles)
        smiles_canonicalized = Chem.MolToSmiles(mol)
        
    except:
        smiles_canonicalized = smiles
        print(f"Error: Invalid SMILES {smiles}")
        flag = False
        
    return smiles_canonicalized, flag


def cleanup_library_by_CSP(df, smiles_column_name, cleanup_chirality = False):
    """
    Clean up smiles with GChem ChEMBL_Structure_Pipeline, add a new column 'Cleaned_SMILES', remove chirality in SMILES (optional).
    :param df: pandas.DataFrame object, input dataframe.
    :param smiles_column_name: str, the name of the SMILES column.
    :param cleanup_chirality: bool, whether to remove chirality or not.
    """
    # add 'Cleaned_SMILES' to columns
    columns = df.columns.tolist()    
    if 'Cleaned_SMILES' in columns:
        warnings.warn('Cleaned_SMILES already exists!')
    smiles_column_num = columns.index(smiles_column_name)
    columns.insert(smiles_column_num + 1, 'Cleaned_SMILES')
        
    # use ChEMBL_Structure_Pipeline to clean up smiles    
    df['Cleaned_SMILES'] = df[smiles_column_name].apply(lambda smiles: cleanup_single_smiles_by_CSP(smiles, cleanup_chirality))
    df[['Cleaned_SMILES', 'filter']] = pd.DataFrame(df['Cleaned_SMILES'].values.tolist())

    # generate df_cleaned
    df_cleaned = df[df['filter']]
    df_cleaned = pd.DataFrame(df_cleaned, columns = columns)
    df_cleaned = df_cleaned.reset_index(drop = True)
    
    # generate df_error
    df_error = df[~df['filter']]
    df_error = pd.DataFrame(df_error, columns = columns)
    df_error = df_error.reset_index(drop = True)
    
    return df_cleaned, df_error


def cleanup_disconnection_in_single_smiles(smiles, saltRemover, process_disconnection_method):
    """
    Record and process a single disconnected SMILES (containing '.', i.e., polymer, salt, solvent).
    :param smiles: str, SMILES.
    :param saltRemover: rdkit SaltRemover object.
    :param process_disconnection_method: str, method for processing other disconnected SMILES,
    if process_disconnection_method == 'keep_longest', keep the longest part in SMILES;
    if process_disconnection_method == 'keep_most_atoms', keep the part with the most atoms.
    """
    if '.' in smiles:
        mol = Chem.MolFromSmiles(smiles)
        
        # remove simple ions in salts
        mol = saltRemover.StripMol(mol)
        smiles = Chem.MolToSmiles(mol)
        
        # process other disconnected SMILES
        if '.' in smiles and process_disconnection_method == 'keep_longest':
            smiles = max(smiles.split('.'), key = len)

        elif '.' in smiles and process_disconnection_method == 'keep_most_atoms':
            fragment_list = smiles.split('.')
            fragment_length_list = [len(Chem.MolFromSmiles(fragment).GetAtoms()) for fragment in fragment_list]
            smiles = fragment_list[np.argmax(fragment_length_list)]
        
        # canonicalize and check SMILES
        mol = Chem.MolFromSmiles(smiles)
        smiles_canonicalized = Chem.MolToSmiles(mol)
        if len(smiles_canonicalized) == 0:
            return None, True
        return smiles_canonicalized, True
    
    else:
        return smiles, False


def cleanup_disconnection_in_library(df, process_disconnection = False, process_disconnection_method = None):
    """
    Record and process disconnected SMILES (containing '.', i.e., polymer, salt, solvent) (optional).
    :param df: pandas.DataFrame object, input dataframe.
    :param process_disconnection: bool, whether to process disconnected SMILES or not.
    :param process_disconnection_method: str, method for processing other disconnected SMILES,
    if process_disconnection_method == 'keep_longest', keep the longest part in SMILES;
    if process_disconnection_method == 'keep_most_atoms', keep the part with the most atoms.
    """
    columns = df.columns.tolist()
    
    # record and process disconnected SMILES
    if process_disconnection:
        ion_list = "[Li,Na,K,Rb,Cs,Mg,Ca,Sr,Ba,Mn,Fe,Co,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,F,Cl,Br,I]"
        saltRemover = SaltRemover(defnData = ion_list)
    
        df['Cleaned_SMILES'] = df['Cleaned_SMILES'].apply(lambda smiles: cleanup_disconnection_in_single_smiles(smiles, saltRemover, process_disconnection_method))
        df[['Cleaned_SMILES', 'filter']] = pd.DataFrame(df['Cleaned_SMILES'].values.tolist())
    
    # record disconnected SMILES
    else:
        df['filter'] = df['Cleaned_SMILES'].apply(lambda smiles: '.' in smiles)
    
    # generate disconnected SMILES list, df_disconnected
    df_disconnected = df[df['filter']]
    df_disconnected = pd.DataFrame(df_disconnected, columns = columns)
    df_disconnected = df_disconnected.reset_index(drop = True)

    # generate processed SMILES list, df
    df.dropna(subset = ['Cleaned_SMILES'], how = 'all', inplace = True)
    df = pd.DataFrame(df, columns = columns)
    df = df.reset_index(drop = True)
    
    return df, df_disconnected


def remove_unnamed_columns(df):
    """
    Remove unnamed columns.
    """
    unnamed_cols = df.columns.str.contains('Unnamed:')
    unnamed_cols_name = df.columns[unnamed_cols]
    df.drop(unnamed_cols_name, axis=1, inplace=True)
    return df


def cleanup_smiles(input_file, smiles_column_name, cleanup_chirality = False, process_disconnection = False, process_disconnection_method = None):
    """
    Clean up smiles with GChem ChEMBL_Structure_Pipeline, add a new column 'Cleaned_SMILES', remove chirality in SMILES (optional),
    record and process disconnected SMILES (containing '.', i.e., polymer, salt, solvent) (optional).
    :param input_file: str, path of the input file.
    :param smiles_column_name: str, the name of the SMILES column.
    :param cleanup_chirality: bool, whether to remove chirality or not.
    :param process_disconnection: bool, whether to process disconnected SMILES or not.
    :param process_disconnection_method: str, method for processing other disconnected SMILES,
    if process_disconnection_method == 'keep_longest', keep the longest part in SMILES;
    if process_disconnection_method == 'keep_most_atoms', keep the part with the most atoms.
    """
    # output file path without extension
    output_file, fmt = os.path.splitext(os.path.abspath(input_file))
    
    if fmt in {'.csv'}:
        df = pd.read_csv(input_file)
    elif fmt in {'.xlsx'}:
        df = pd.read_excel(input_file)
    else:
        print('Error: Invalid format of the input file')
        return
    
    # clean up smiles with GChem ChEMBL_Structure_Pipeline, add a new column 'Cleaned_SMILES', remove chirality in SMILES (optional)
    df, df_error = cleanup_library_by_CSP(df, smiles_column_name, cleanup_chirality)
    
    # record and process disconnected SMILES (containing '.', i.e., polymer, salt, solvent) (optional)
    df, df_disconnected = cleanup_disconnection_in_library(df, process_disconnection, process_disconnection_method)
    
    # write to file
    df = df.reset_index(drop = True)
    print('Number of rows after removing invalid SMILES:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv('{}_CSP.csv'.format(output_file))

    df_error = df_error.reset_index(drop = True)
    print('Number of error SMILES:', df_error.shape[0])
    df_error = remove_unnamed_columns(df_error)
    df_error.to_csv('{}_SMILES_errors.csv'.format(output_file))

    df_disconnected = df_disconnected.reset_index(drop = True)
    print('Number of SMILES to check:', df_disconnected.shape[0])
    df_disconnected = remove_unnamed_columns(df_disconnected)
    df_disconnected.to_csv('{}_SMILES_to_check.csv'.format(output_file))



if __name__ == '__main__':
    
    input_file = 'tests/example_format.csv'
    smiles_column_name = 'SMILES'
    cleanup_chirality = True
    process_disconnection = True
    process_disconnection_method = 'keep_most_atoms'
    
    cleanup_smiles(input_file, smiles_column_name, cleanup_chirality, process_disconnection, process_disconnection_method)

