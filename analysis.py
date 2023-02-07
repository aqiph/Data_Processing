#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:21:09 2021

@author: guohan
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_size(5)



def analyze_assay(input_file, assay_column_name, value_column_name):
    """
    Get unique assays in the 'assay_column_name' column, box plot of each assay
    :param input_file: str, the filename of the input file
    :param assay_column_name: str, name of assay column to be analyzed
    :param value_column_name: str, name of value column to be analyzed
    """
    # output file path without extension
    output_file, fmt = os.path.splitext(os.path.abspath(input_file))

    if fmt in {'.csv'}:
        df = pd.read_csv(input_file, index_col=0)
    elif fmt in {'.xlsx'}:
        df = pd.read_excel(input_file)
    else:
        print('Error: Invalid format of the input file')
        return

    # change assay_column_name to float
    df[value_column_name] = df[value_column_name].apply(to_float)
    
    # get a list of values in the 'assay_column_name' column, _assay_type.csv
    output_file_counts = '{}_{}_values.csv'.format(output_file, assay_column_name)
    df_assay_counts = df[assay_column_name].value_counts().to_frame()
    df_assay_counts.to_csv(output_file_counts)
    
    # box plot of the 'value_column_name' column on the 'assay_column_name' column
    output_file_boxplot = '{}_{}_boxplot.pdf'.format(output_file, assay_column_name)
    bplot = df.boxplot(column = value_column_name, by = assay_column_name)
    plt.xticks(fontproperties=font, rotation='vertical')
    plt.yticks(fontproperties=font)
    plt.ylim([0.0, 2000])
    bplot.figure.savefig(output_file_boxplot, format = 'pdf')

 
def value_distribution(input_file, value_column_name, useLog = False, range = [0, 10]):
    """
    plot distribution for values in the 'value_column_name' column
    :param input_file: str, the filename of the input file
    :param value_column_name: str, name of value column to be analyzed
    :param useLog: bool, whether to use log or raw value
    :param range: list of two ints, range of the plot
    """
    # output file path without extension
    output_file, fmt = os.path.splitext(os.path.abspath(input_file))

    if fmt in {'.csv'}:
        df = pd.read_csv(input_file, index_col=0)
    elif fmt in {'.xlsx'}:
        df = pd.read_excel(input_file)
    else:
        print('Error: Invalid format of the input file')
        return

    # get values
    print('Columns are ', df.columns)
    df[value_column_name] = df[value_column_name].apply(to_float)
    values = df[value_column_name].values.tolist()
    
    if useLog:
        values = np.log10(values)
    
    if True:
        print('maximum:', max(values), 'minimum:', min(values))
    
    # plot distribution
    output_file_distribution = '{}_{}_distribution.pdf'.format(output_file, value_column_name)
    
    plt.figure(1)
    plt.hist(values, 50, range = range)
    plt.xlabel(value_column_name, fontproperties=font)
    plt.ylabel('Counts', fontproperties=font)
    plt.xticks(fontproperties=font)
    plt.yticks(fontproperties=font)

    
    plt.savefig(output_file_distribution, dpi = 300)
    plt.show()
    plt.close()

    return df


def to_float(v):
    try:
        v = float(v)
        return v
    except:
        return np.nan
    
    

if __name__ == '__main__':

    # analyze assay
    # input_file = 'tests/example_format.csv'
    # assay_column_name = 'Assay_Parameter'
    # value_column_name = 'Value'
    # analyze_assay(input_file, assay_column_name, value_column_name)

    # plot value distribution
    input_file = 'tests/example_format.csv'
    value_column_name = 'Value'
    value_distribution(input_file, value_column_name, useLog = False, range = [-10, 10])