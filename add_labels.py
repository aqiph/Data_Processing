#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 11:43:28 2022

@author: guohan
"""

import os
import pandas as pd
import numpy as np
from functools import partial



def add_label_column(input_file, task, **kwargs):
    """
    add labels based on the Value column
    :param input_file: str, the filename of the input file
    :param task: str, task for labeling, including 'classification', 'regression'
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
    
    # add labels
    task = task.lower()
    if task == 'classification':
        value_column_name = kwargs.pop('value_column_name', 'Value')
        thresholds = kwargs.pop('thresholds', [0])
        value_label_correlation = kwargs.pop('value_label_correlation', 'negative')
        label_fn = partial(ClassificationLabeling, value_column_name=value_column_name, thresholds=thresholds, value_label_correlation=value_label_correlation)
    elif task == 'regression':
        value_column_name = kwargs.pop('value_column_name', 'Value')
        regression_label_function = kwargs.pop('regression_label_function', 'linear').lower()
        label_fn = partial(RegressionLabeling, value_column_name=value_column_name, regression_label_function=regression_label_function)
    else:
        'Error: Input function for labeling is not defined'
        return

    df['Label'] = df.apply(label_fn, axis = 1)
    
    # write to file
    df = df.reset_index(drop = True)
    print('Number of rows after adding labels:', df.shape[0])
    if task in {'classification'}:
        print('Number of examples in each class:', df['Label'].value_counts())
    df = remove_unnamed_columns(df)
    df.to_csv('{}_labels.csv'.format(output_file))


# Functions for labeling

def ClassificationLabeling(row, value_column_name, thresholds, value_label_correlation):
    """
    Helper function for add_label_column(). Generate multi-class classification label based on value
    :param row: row of DataFrame
    :param value_column_name: str, name of the value column
    :param thresholds: list of floats, monotonically descending or ascending
    :param value_label_correlation: str, 'positive' indicates positive correlation between values and labels;
    'negative' indicates negative correlation between values and labels.
    """
    value = row[value_column_name]
    value_label_correlation = value_label_correlation.lower()
    label = 0

    # positive correlation: class 0 -- smallest values, class n --- largest values
    if value_label_correlation == 'positive':
        thresholds.sort()
        while label < len(thresholds):
            if value < float(thresholds[label]):
                return label 
            label += 1
        return label
    
    # negative correlation: class 0 --- largest values, class n --- smallest values
    elif value_label_correlation == 'negative':
        thresholds.sort(reverse = True)
        while label < len(thresholds):
            if value >= float(thresholds[label]):
                return label
            label += 1
        return label

    else:
        print('Error: Invalid input correlation between values and labels')
        return np.nan


def RegressionLabeling(row, value_column_name, regression_label_function):
    """
    Helper function for add_label_column(). Generate regression label based on value
    :param row: row of DataFrame
    :param value_column_name: str, name of the value column
    :param regression_label_function: str, name of the function used to convert value to regression label, can be: 'linear', 'quadratic'
    """
    value = row[value_column_name]
    regression_label_fn = regression_label_function.lower()

    if regression_label_fn == 'linear':
        return float(value)

    elif regression_label_fn == 'quadratic':
        return - float(value) ** 2.0 if value <= 0.0 else 0


def remove_unnamed_columns(df):
    """
    remove unnamed columns
    """
    unnamed_cols = df.columns.str.contains('Unnamed:')
    unnamed_cols_name = df.columns[unnamed_cols]
    df.drop(unnamed_cols_name, axis=1, inplace=True)
    return df



if __name__ == '__main__':
    
    input_file = 'tests/example_noIndex_format_CSP_rmDuplicates.csv'
    task = 'classification'

    value_column_name = 'Value'
    thresholds = [10.0]
    value_label_correlation = 'negative'
    regression_label_function = 'linear'

    add_label_column(input_file, task, value_column_name=value_column_name, thresholds=thresholds, value_label_correlation=value_label_correlation,
                     regression_label_function= regression_label_function)
    
