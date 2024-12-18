#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 13:36:00 2023

@author: guohan

"""

import os
import pandas as pd
import numpy as np



def remove_unnamed_columns(df):
    """
    Remove unnamed columns.
    """
    unnamed_cols = df.columns.str.contains('Unnamed:')
    unnamed_cols_name = df.columns[unnamed_cols]
    df.drop(unnamed_cols_name, axis=1, inplace=True)
    return df