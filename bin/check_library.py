#!/usr/bin/env python

import sys
import os
import argparse
import statistics
from Bio import SeqIO
#from cogent3 import make_unaligned_seqs,get_distance_calculator
#from cogent3.core.alignment import Alignment
#from cogent3.evolve.models import HKY85
import numpy as np
import pandas as pd
#import seaborn as sns
#import matplotlib
#import matplotlib.pyplot as plt
import warnings
#matplotlib.use('Agg')
#matplotlib.style.use('ggplot')
warnings.filterwarnings("ignore")

# Pad sequences that are too short
def pad_read(read, max_length):
    padseq = read.seq + (max_length-len(read.seq)) * '-'
    return str(padseq)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--library', dest='library', type=str, help='Input text file for library')
    parser.add_argument('-m', '--mode', dest='mode', type=str, default = 'fasta', help="Mode for checking. Can be 'txt' or 'fasta'")
    args = parser.parse_args()

    # Parse FASTA file
    if args.mode == "fasta":
        record_dict = SeqIO.to_dict(SeqIO.parse(args.library, "fasta"))
        rec_len = [ len(record.seq) for id, record in record_dict.items() ]

        # Log basic stats about the library
        print('Number of records:', len(rec_len))
        print('Minimum length (nt):', min(rec_len))
        print('Maximum length (nt):', max(rec_len))
        print('Mean length (nt):', round(statistics.mean(rec_len), 2) )
        print('Mode length (nt):', statistics.mode(rec_len) )
        library_mode = statistics.mode(rec_len)

        # Print histogram of length distribution
        #sns_plot = sns.distplot(rec_len, kde=False, rug=True).set_title('Distribution of library read lengths (nt)')
        #fig = sns_plot.get_figure()
        #fig.savefig("library-histogram.png")

        # Pad short sequences
        #pad_record_dict = {}
        #for id, record in record_dict.items():
        #    if len(record.seq) < library_mode:
        #        case = { id: pad_read(record, library_mode) }
        #    elif len(record.seq) > library_mode:
        #        print("Read is longer than the mode lengths")
         #   else:
         #       case = { id: str(record.seq) }
         #   pad_record_dict.update(case)

        # Import & Estimate distances
        #seq = make_unaligned_seqs(pad_record_dict, moltype="dna")
        #aln = Alignment(seq)
        #dists = aln.distance_matrix(calc="hamming", show_progress=True)
        #dist_df = pd.DataFrame(dists)
        #dist_df.head()
        #dist_calc.run( parallel=True )
        #dists = dist_calc.get_pairwise_distances()
        #mycluster = upgma(dist_calc.get_pairwise_distances())
        #mycluster.write('test_upgma.tree')

        #print(dists)
