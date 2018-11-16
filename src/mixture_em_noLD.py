#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd

# global variables
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def log_likelihood(beta_tilde, gamma, sigma_e, W):
    M = len(beta_tilde)
    mu = np.multiply(W, gamma)
    cov = np.multiply(np.eye(M), sigma_e)
    print beta_tilde.shape
    print mu.shape
    print cov.shape
    st.multivariate_normal.logpdf(x=beta_tilde, mean=mu, cov=cov)
    return log_like


def initialize_p(K):
    # random draw to initialize values
    p_init = np.random.dirichlet([1]*K,1)
    p_init = p_init.ravel()
    #p_init = [.5, .5]
    return p_init

def EM(p_t, mu_vec, sigma_vec, beta_tilde, sigma_e, its):
    K = len(p_t)
    M = len(beta_tilde)

    # print initial value
    print p_t

    for i in range(0, its):
        # Expectation step
        r_MK = np.empty((M,K))

        for m in range(0, M):
            z = np.empty(K)
            for k in range(0, K):
                like = st.norm.pdf(beta_tilde[m], mu_vec[k], sigma_vec[k] + sigma_e)
       
                if np.nan == like:
                    like = 0 
                    logging.info("Encountered Nan!")
                z[k] = p_t[k]*like

                if np.isnan(z[k]):
                    print "Iteration: %d" % i 
                    print p_t[k]
                    print like
                    print z[k]
                    exit(1)
                #print p_t[k]
                

#            print np.sum(z)
            for k in range(0, K):
                if np.sum(z) > 0:
                    r_mk = z[k]/np.sum(z)
                else:
#                    print "sum z is equal to zero"
                    r_mk = 0 

                r_MK[m,k] = r_mk

        # Maximization step
        for k in range(0, K):
            if np.sum(r_MK) > 0:
                p_t[k] = np.sum(r_MK[:,k])/np.sum(r_MK)
            else:
                p_t[k] = 0 

        print "Iteration %d: " % i 
        print p_t

    return p_t


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--mu_vec", dest="mu_vec",)
    parser.add_option("--sigma_vec", dest="sigma_vec")
    parser.add_option("--ldsc_h2", dest="ldsc_h2")
    parser.add_option("--bins", dest="bins")
    parser.add_option("--ld_half_file", dest="ld_half_file")
    parser.add_option("--N", dest="N", default=100000)
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")
    parser.add_option("--precompute", dest="precompute", default='y')
    parser.add_option("--its", dest="its", default=500)
    (options, args) = parser.parse_args()

    # parse command line args
    seed = int(options.seed)
    its = int(options.its)
    N = int(options.N)
    name = options.name
    precompute = options.precompute
    gwas_file = options.gwas_file
    outdir = options.outdir
    ldsc_h2 = float(options.ldsc_h2)

    if options.bins is None:
        mu_vec = [float(item) for item in options.mu_vec.split(',')]
        sigma_vec = [float(item) for item in options.sigma_vec.split(',')]
    else:
        bins = int(options.bins)

        # get LHS of bins
        a = -1.0
        b = 1.0
        step = (b-a)/float(bins)
        sigma_k = ((step*.5)/float(3))**2
        sigma_vec = np.repeat(sigma_k, bins)
        mu_vec = np.empty(bins)
        for k in range(bins):
            if k == 0:
                mu_vec[k] = a + step*.50
            else:
                mu_vec[k] = mu_vec[k-1] + step
        print mu_vec 
        print sigma_vec

    logging.info("Reading in gwas file: %s" % gwas_file)
    df = pd.read_csv(gwas_file, sep=' ')

    ld_half_file = options.ld_half_file

    if ld_half_file is not None:
        beta_tilde = np.asarray(df['BETA_STD_I'])
    else:
        logging.info("Assuming no LD")
        beta_tilde = np.asarray(df['BETA_STD'])

    # set seed
    np.random.seed(seed)


    if ld_half_file is not None:
        W = np.loadtxt(ld_half_file)
        logging.info("Using ld half file: %s" % ld_half_file)

    # get meta values from data
    K = len(mu_vec)
    M = beta_tilde.shape[0]

    # intialize values of chain
    logging.info("initializing p")
    p_init = initialize_p(K)

    # calculate sigma_e 
    sigma_e = (1-ldsc_h2)/N

    p_est = EM(p_init, mu_vec, sigma_vec, beta_tilde, sigma_e, its)
    print "Estimate: "
    print p_est

if __name__== "__main__":
  main()
