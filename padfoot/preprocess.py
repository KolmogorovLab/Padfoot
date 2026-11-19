#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 17:39:44 2025

@author: keskusa2
"""

import pandas as pd

def generate_rm(rm_file, rm_path):
    df = pd.read_csv(rm_file, sep='\t', header=None)
    df_small = df[[4, 5, 6, 10]]
    df_small.to_csv(rm_path, sep='\t', header=False, index=False)
    