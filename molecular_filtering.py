#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 10:24:00 2023

@author: guohan
"""

import os, sys
import argparse
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from tools import remove_unnamed_columns



### Helper functions ###

def check_SMILES(smiles):
    if not smiles:
        print(f"Error: Invalid SMILES {smiles}")
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        print(f"Error: Invalid SMILES {smiles}")
        return False
    if not mol:
        print(f"Error: Invalid SMILES {smiles}")
        return False
    return True


### main ###

def main(args):
    """
    Filter molecules based on given rules
    """
    # output file
    folder, basename = os.path.split(os.path.abspath(args.input_file))
    output_file = args.output_file
    if output_file == '':
        output_file = basename
    output_file = os.path.join(folder, os.path.splitext(output_file)[0])

    # read file
    df = pd.read_csv(args.input_file)
    COLUMNS = df.columns.tolist()

    # convert to mol
    df = pd.DataFrame(df[df[args.smiles_column_name].apply(check_SMILES)])
    df['Mol'] = df[args.smiles_column_name].apply(Chem.MolFromSmiles)

    # molecular filtering
    # molecular weight
    if args.rule_MW_UB is not np.nan:
        df = pd.DataFrame(df[df['Mol'].apply(lambda smiles: ExactMolWt(smiles) <= args.rule_MW_UB)])
        if False:
            df['MW'] = df['Mol'].apply(ExactMolWt)
            COLUMNS.append('MW')

    # write output file
    df = pd.DataFrame(df, columns = COLUMNS)
    df = df.reset_index(drop=True)
    print('Number of SMILES:', df.shape[0])
    output_file = '{}_filtered_{}.csv'.format(output_file, df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(output_file)


def get_parser():
    """
    generate parser
    """
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--input_file', default='tests/molecular_filtering_test.csv', type=str, help='Path of the input file.')
    argparser.add_argument('--output_file', default='', type=str, help='File name of the output file.')
    argparser.add_argument('--smiles_column_name', default='Cleaned_SMILES', type=str, help='Name of the SMILES column.')

    argparser.add_argument('--rule_MW_UB', default=np.nan, type=float, help='Rule: upper bound of molecular weight.')


    args = argparser.parse_args()

    return args


if __name__ == '__main__':
    args = get_parser()
    main(args)