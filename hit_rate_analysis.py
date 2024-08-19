#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 6 14:08:00 2023

@author: guohan
"""

import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_size(12)
path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Data_Processing'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')
from tools import remove_unnamed_columns


def get_hit_count(input_file, input_library, id_column_name='ID', print_hits=False):
    """
    Calculate the number of hits
    :param input_file: str, path of the file containing selected/active compounds
    :param input_library: str, path of the library file
    :param id_column_name: str, name of the ID column in the input_file
    :param print_hits: bool, whether to print the hits in this library or not
    """
    df_active = pd.read_csv(input_file)
    df_library = pd.read_csv(input_library)
    sublibrary = str(os.path.split(input_library)[1]).split('_')[1]
    output = os.path.splitext(os.path.abspath(input_file))[0]
    output = f'{output}_{sublibrary}'

    COLUMNS = df_active.columns.tolist()
    df_active['ID_'] = df_active[id_column_name].apply(lambda id: str(id))
    df_library['ID_'] = df_library['ID'].apply(lambda id: str(id))
    df_library.rename(columns={'ID':'ID_lib', 'SMILES':'SMILES_lib', 'Cleaned_SMILES':'Cleaned_SMILES_lib',
                               'Assay':'Assay_lib', 'Assay_Parameter':'Assay_Parameter_lib', 'Operator':'Operator_lib',
                               'Value':'Value_lib', 'Units':'Units_lib'}, inplace=True)

    df_int = pd.merge(df_active, df_library, how='inner', on=['ID_'])
    df_int = pd.DataFrame(df_int, columns = COLUMNS)
    hit_count = df_int.shape[0]
    print(f'Number of cmps in {input_library} is {hit_count}')
    if print_hits:
        df_int = remove_unnamed_columns(df_int)
        df_int.to_csv(f'{output}_{hit_count}.csv')

    return hit_count


def hit_counts_in_HTS(input_file, id_column_name='ID', print_hits=False):
    """
    Calculate the number of hits and hit rate in each sublibrary of HTS library
    :param input_file: str, path of the file containing selected/active compounds
    :param id_column_name: str, name of the ID column in the input_file
    :param print_hits: bool, whether to print the hits in this library or not
    """
    # output file path without extension
    output_file = os.path.splitext(os.path.abspath(input_file))[0]

    folder = '/Users/guohan/Documents/Projects/Datasets/HTS/'
    Library = ['Disease-Based-Library', 'Fragment-Library', 'High-Throughput-Library',
               'Nature-Product-Library', 'Repurposing-Library']
    Sublibrary = [['Anti-Virus', 'TBAC-Alliance', 'TBAC-Calibr'],
                  ['Fragment-Enamine', 'Fragment-MayBridge'],
                  ['AnalytiCon-NATx', 'ChemBridgeCL20K', 'ChemBridgeCL90K', 'ChemBridgeEXP50K', 'ChemDiv-SMART-TM',
                   'Expresspick-Pfizer', 'LifeChemical-30K', 'MayBridge-50K', 'Selleck-PFZ'],
                  ['AnalytiCon-MEGx-NP', 'MetaSci-HML-NP', 'Pharmacodia-NP', 'TargetMol-NP'],
                  ['FDA', 'Prc&Cli', 'ReFrame']]
    Sublibrary_tot = [[186, 1364, 861],
                      [1800, 2490],
                      [17947, 19999, 90237, 49974, 53200, 4208, 30080, 52531, 99040],
                      [3520, 837, 868, 2130],
                      [2664, 1129, 11599]]

    Sublib_names = []
    Hit_counts = []
    Hit_rates = []
    for i, library in enumerate(Library):
        for j, sublibrary in enumerate(Sublibrary[i]):
            files = os.listdir(f'{folder}/{library}/{sublibrary}_20230317/final')
            filtered_files = [file for file in files if file.startswith(f'HTS_{sublibrary}_forGeneralUse_')]
            input_library = f'{folder}/{library}/{sublibrary}_20230317/final/{filtered_files[0]}'

            hit_count = get_hit_count(input_file, input_library, id_column_name, print_hits=print_hits)
            hit_rate = float(hit_count)/float(Sublibrary_tot[i][j])
            Sublib_names.append(sublibrary)
            Hit_counts.append(hit_count)
            Hit_rates.append(f'{hit_rate:.4f}')

    # write to file
    df = pd.DataFrame({'Library Name': Sublib_names, 'Hit Counts': Hit_counts, 'Hit Rate': Hit_rates})
    print('The number of rows:', df.shape[0])
    df = remove_unnamed_columns(df)
    df.to_csv(f'{output_file}_stat.csv')


def plot_stat(input_file, remove_zero=True):
    """
    plot hit counts and hit rates
    :param input_file: str, path of the file containing library name, hit counts and hit rates
    :param remove_zero: bool, whether to remove rows with zero count or not
    """
    # output file path without extension
    output_file = os.path.splitext(os.path.abspath(input_file))[0]

    df = pd.read_csv(input_file)
    if remove_zero:
        df = pd.DataFrame(df[df['Hit Counts'].apply(lambda hitCounts: int(hitCounts) != 0)])

    # plot
    hit_counts = df['Hit Counts'].tolist()
    library_names = df['Library Name'].tolist()
    hit_rates = df['Hit Rate'].tolist()

    fig, ax1 = plt.subplots(figsize = (8,6))
    ax1.bar(library_names, hit_counts, color = 'blue', edgecolor='black')
    for i, count in enumerate(hit_counts):
        ax1.text(i, count + 0.5, str(count), ha='center')
    ax1.set_xlabel('Library', fontproperties=font)
    ax1.set_ylabel('Counts', fontproperties=font)
    plt.xticks(rotation=80)

    ax2 = ax1.twinx()
    ax2.scatter(library_names, hit_rates, color='red')
    ax2.plot(library_names, hit_rates, color='black', linestyle='--')
    ax2.set_ylabel('Hit Rate', fontproperties=font)

    plt.savefig('{}_plot.pdf'.format(output_file), dpi=300, bbox_inches = 'tight')
    plt.show()
    plt.close()



if __name__ == '__main__':
    input_file = 'tests/test_hit_rate_analysis.csv'
    hit_counts_in_HTS(input_file, id_column_name='Analog_ID')

    # input_file = 'tests/test_hit_rate_analysis_stat.csv'
    # plot_stat(input_file, remove_zero=True)