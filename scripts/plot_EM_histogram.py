#!/usr/bin/env python

from optparse import OptionParser
import sys
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--results_file", dest="results_file")
    parser.add_option("--outdir", dest="outdir")
    (options, args) = parser.parse_args()

    # parse arguments
    name = options.name
    results_file = options.results_file
    outdir = options.outdir
    if outdir == None:
        outdir=""

    # read in file
    results_df = pd.read_csv(results_file, sep=' ')

    # plot
    sns.set_style('whitegrid')
    a = -1
    b = 1

    results_df.plot.bar(x='mu', y='p', legend=False)

    # save plot
    outfile=os.path.join(outdir, name + '.barplot.pdf')
    plt.title("Histogram of binned effect sizes - %s" % name, fontsize=15)

    plt.savefig(outfile)

if __name__== "__main__":
  main()
