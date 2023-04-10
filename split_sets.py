#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: guohan
"""

import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split



### Splitter ###

class Splitter(object):
    """
    Splitter
    """

    def __init__(self, X, stratify = None, seed = None):
        self.X = X
        self.stratify = stratify
        self.seed = seed

        if self.stratify is not None:
            assert len(self.X) == len(self.stratify), 'Error: X and label should be the same length.'


    def train_test(self, trainset_ratio):
        """
        split X into trainset and testset, if self.stratify is not None, data is split in a stratified fashion
        :param trainset_ratio: float, ratio of trainset
        """
        assert trainset_ratio <= 1.0 and trainset_ratio >= 0.0, 'Error: trainset_ratio should between 0.0 and 1.0.'
        X_train, X_test = train_test_split(self.X, train_size = trainset_ratio, random_state = self.seed, stratify = self.stratify)
        return X_train, X_test


### helper function ###

def indexes_to_bools(indexes, num_tot):
    assert num_tot > np.max(indexes)

    bools = [False] * num_tot
    for i in indexes:
        bools[i] = True

    return bools


def remove_unnamed_columns(df):
    """
    remove unnamed columns
    """
    unnamed_cols = df.columns.str.contains('Unnamed:')
    unnamed_cols_name = df.columns[unnamed_cols]
    df.drop(unnamed_cols_name, axis=1, inplace=True)
    return df


### split sets ###

def main(input_file, stratified_by = None, trainset_ratio = 0.8):
    """
    split the set
    :param input_file: str, file path of the input file
    :param stratified_by: str or None, column based on which to split data in a stratified fashion
    :param trainset_ratio: float, ratio of trainset
    """
    # folder
    folder, _ = os.path.split(os.path.abspath(input_file))

    # read file
    df = pd.read_csv(input_file)
    num_rows = df.shape[0]
    X = np.arange(num_rows)
    if stratified_by is not None:
        stratify = df[stratified_by].values.tolist()
    else:
        stratify = None

    # split sets
    splitter = Splitter(X, stratify = stratify)

    X_train_index, X_test_index = splitter.train_test(trainset_ratio = trainset_ratio)
    print(X_test_index)

    df_train_set = df[indexes_to_bools(X_train_index, num_rows)]
    df_test_set = df[indexes_to_bools(X_test_index, num_rows)]

    # output file
    df_train_set = df_train_set.reset_index(drop=True)
    df_test_set = df_test_set.reset_index(drop=True)
    print('Number of rows trainset is {}'.format(df_train_set.shape[0]))
    print('Number of rows testset is {}'.format(df_test_set.shape[0]))
    df_train_set = remove_unnamed_columns(df_train_set)
    df_test_set = remove_unnamed_columns(df_test_set)
    df_train_set.to_csv(os.path.join(folder, 'trainset.csv'))
    df_test_set.to_csv(os.path.join(folder, 'testset.csv'))



if __name__ == '__main__':
    input_file = 'tests/example_noIndex_split_train_test.csv'
    main(input_file, stratified_by='Label', trainset_ratio = 0.8)