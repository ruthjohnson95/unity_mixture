#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)


def simulate_mixture(p_vec, mu_vec, sigma_vec, M, N, V):
    K = len(p_vec)
    k_matrix = np.random.multinomial(1, p_vec, size=M)

    # make large beta_k_matrix
    beta_k_matrix = np.zeros((M,K))

    # sample from K gaussians
    for k in range(0, K):
        beta_k_matrix[:,k] = st.norm.rvs(loc=mu_vec[k],scale=sigma_vec[k],size=M)

    # only pick betas from chosen components
    beta_k = np.multiply(beta_k_matrix,  k_matrix)

    beta = np.sum(beta_k, axis=1)

    # sample GWAS effects
    M_k = np.multiply(M, p_vec)
    sigma_g = np.sum(np.multiply(M_k, sigma_vec))
    sigma_e = (1-sigma_g)/float(N)
    mu = np.matmul(V, beta)
    cov = np.multiply(V, sigma_e)

    beta_hat = st.multivariate_normal.rvs(mu, cov)

    df = {'BETA_STD': beta_hat}
    beta_hat_df = pd.DataFrame(data=df)

    return beta_hat_df


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--p_vec", dest="p_vec", nargs='+', default=".25 .25 .50")
    parser.add_option("--mu_vec", dest="mu_vec", nargs='+', default="-.10 0 .10")
    parser.add_option("--sigma_vec", dest="sigma_vec", nargs="+", default=".001 .001 .001")
    parser.add_option("--bins", dest="bins")
    parser.add_option("--ld_file", dest="ld_file")
    parser.add_option("--M", dest="M", default=10)
    parser.add_option("--N", dest="N", default=100000)
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")

    (options, args) = parser.parse_args()

    # parse input options
    if options.bins is None:
        p_vec= [float(item) for item in options.p_vec.split()]
        mu_vec = [float(item) for item in options.mu_vec.split()]
        sigma_vec = [float(item) for item in options.sigma_vec.split()]
    else:
        bins = int(bins)

    outname = options.name
    ld_file = options.ld_file
    M = int(options.M)
    N = int(options.N)
    seed = int(options.seed)
    outdir = options.outdir

    # check that mixture components are same size
    if len(p_vec) == len(mu_vec) == len(sigma_vec):
        pass
    else:
        logging.info("ERROR: p, mu, sigma lists are NOT same size...exiting")
        exit(1)

    # simulate values
    if ld_file is None:
        V = np.eye(M)
    else:
        V = np.loadtxt(ld_file)

    beta_hats = simulate_mixture(p_vec, mu_vec, sigma_vec, M, N, V)

    # save to outfile
    outfile = os.path.join(outdir, outname + '.' + str(seed) + '.txt')
    beta_hats.to_csv(outfile, index=False)

    logging.info("DONE...simulated files can be found at: %s" % outfile)


if __name__== "__main__":
  main()
