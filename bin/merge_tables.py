#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def parse_pileup(file_loc):
    basename = file_loc.replace('-counts.txt', '')
    df = pd.read_table(file_loc, sep='\t')
    df['SampleId'] = basename
    df = df[['SampleId', '#ID', 'Length', 'Plus_reads']]
    df.columns = ['SampleId', 'Id', 'Length', 'Reads']
    return(df)

def parse_rpkm(file_loc):
    basename = file_loc.replace('-rpkm.txt', '')
    df = pd.read_table(file_loc, sep='\t', skiprows=4)
    df['SampleId'] = basename
    df = df[['SampleId', '#Name', 'Length', 'Reads']]
    df.columns = ['SampleId', 'Id', 'Length', 'Reads']
    return(df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--counts', dest='counts', type=str, nargs='+', help='Set of multiple input files from BBmap alignment.')
    parser.add_argument('-o', '--output', dest='output', type=str, default = 'merged-counts.txt', help="Final output of merged tables")
    parser.add_argument('-r', '--rpkm', dest='rpkm', action='store_true', help="Final output of merged tables")

    args = parser.parse_args()

    # Parse each counts table + concat
    if args.rpkm:
        dfs = [parse_rpkm(file) for file in args.counts]
    else:
        dfs = [parse_pileup(file) for file in args.counts]

    # Pivot wider to have a matrix of counts
    wide_df = pd.concat(dfs).pivot(index='Id', columns='SampleId', values='Reads')

    # Write table
    wide_df.to_csv('merged-counts.txt', sep='\t')

    # Make a plot comparing each pairwise counts
    sns_plot = sns.pairplot(wide_df)
    sns_plot.savefig("merged-counts-plot.png")
