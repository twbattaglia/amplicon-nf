#!/usr/bin/env python

import sys
import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# Import counts and create
def parse_pileup(file_loc):
    basename = file_loc.replace('-counts.txt', '')
    df = pd.read_table(file_loc, sep='\t')
    df['SampleId'] = basename
    df = df[['SampleId', '#ID', 'Length', 'Plus_reads']]
    df.columns = ['SampleId', 'Id', 'Length', 'Reads']
    return(df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', dest='input', type=str, nargs='+', help='Set of multiple input files from BBmap output.')
    parser.add_argument('-o', '--output', dest='mode', type=str, default = 'merged-counts.txt', help="Final output of merged tables")
    #parser.add_argument('-l', '--library', dest='mode', type=str, help="Text file of the library")
    args = parser.parse_args()

    # Parse each counts table + concat
    dfs  = [parse_pileup(file) for file in args.input]

    # Pivot wider to have a matrix of counts
    wide_df = pd.concat(dfs).pivot(index='Id', columns='SampleId', values='Reads')

    # Write table
    wide_df.to_csv('merged-counts.txt', sep='\t')

    # Make a plot comparing each pairwise counts
    sns_plot = sns.pairplot(wide_df)
    sns_plot.savefig("merged-counts-plot.png")
