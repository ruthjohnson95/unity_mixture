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


def compute_sigma_e(p_vec, sigma_vec, M, N):
    M_k = np.multiply(M, p_vec)
    sigma_g = np.sum(np.multiply(M_k, sigma_vec))
    sigma_e = (1-sigma_g)/float(N)
    return sigma_e

def compute_km_denom(mu_km_vec, sigma_km_vec, mu_vec, sigma_vec, p_vec):
    sum = 0
    K = len(mu_km_vec)
    for k in range(0,K):
        sigma_k = sigma_vec[k]
        mu_k = mu_vec[k]
        var_term = np.sqrt(sigma_km_vec[k]/sigma_vec[k])
        mean_term = np.exp((1/(2*sigma_km_vec[k]))*(mu_km_vec[k]**2) - (1/(2*sigma_k)*(mu_k**2)))
        sum += p_vec[k]*var_term*mean_term

    return sum


def compute_q_km(k, m, p_vec, mu_vec, sigma_vec, psi_m, A, C_t, gamma_t, sigma_e):

    # holds means and variances across K components (NOT SNPs)
    K = len(mu_vec)
    mu_km_vec = np.empty((K,1))
    sigma_km_vec = np.empty((K,1))
    q_vec=np.zeros(K)

    for k in range(0, K):
        mu_k = mu_vec[k]
        sigma_k = sigma_vec[k]
        C_k = C_t[:,k]
        gamma_k = gamma_t
        mu_km_vec[k], sigma_km_vec[k] = compute_mu_sigma_km(m, mu_k, sigma_k, psi_m, A, gamma_k, sigma_e)

    #print mu_km_vec

    q_km_denom = compute_km_denom(mu_km_vec, sigma_km_vec, mu_vec, sigma_vec, p_vec)

    for k in range(0, K):
        var_term = np.sqrt(sigma_km_vec[k]/sigma_vec[k])
        sigma_k = sigma_vec[k]
        mu_k = mu_vec[k]
        mean_term = np.exp((1/(2*sigma_km_vec[k]))*(mu_km_vec[k]**2) - (1/(2*sigma_k)*(mu_k**2)))

        q_km_num = p_vec[k]*var_term*mean_term

        q_km = q_km_num/q_km_denom
        q_vec[k] = q_km

    return q_vec


def dp_term(m, A, gamma_k):
    sum = 0
    nonzero_inds = np.nonzero(gamma_k)[0]
    for i in nonzero_inds:
        if i != m:
            sum += A[i,m] * gamma_k[i]
    return sum


def compute_mu_sigma_km(m, mu_k, sigma_k, psi_m, A, gamma_k, sigma_e):

    sigma_km = 1/(1/sigma_k + 1/sigma_e * A[m,m])

    term1 = (sigma_km*mu_k)/sigma_k
    term2 = sigma_km*sigma_e
    dp = dp_term(m, A, gamma_k)
    term3 = psi_m - dp
    mu_km = term1 + term2*term3

    #print "mu_km - SNP %d: %.4g" % (m, mu_km )
    return mu_km, sigma_km


def gibbs(p_init, gamma_init, C_init, mu_vec, sigma_vec, W, A, psi, beta_tilde, N, its):

    # get metadata
    K = len(mu_vec)
    M = len(beta_tilde)

    # make lists to hold samples
    p_list = []
    C_list = np.empty((M, K))

    # start chain
    p_t = p_init
    C_t = C_init
    gamma_t = gamma_init

    gamma_km = np.empty(K)

    # start sampler
    logging.info("Starting sampler")
    for i in range(0, its):
        for m in range(0, M):

            for k in range(0, K):
                # compute posterior mean and variance
                sigma_e = compute_sigma_e(p_t, sigma_vec, M, N)
                mu_km, sigma_km = compute_mu_sigma_km(m, mu_vec[k], sigma_vec[k], psi[m], A, gamma_t, sigma_e)

                gamma_km[k] = st.norm.rvs(mu_km, sigma_km)

            # sample mixture assignments
            #print "SNP %d" % m
            q_km = compute_q_km(k, m, p_t, mu_vec, sigma_vec, psi[m], A, C_t, gamma_t, sigma_e)
            #print q_km
            C  = st.multinomial.rvs(n=1, p=q_km, size=1)
            #print C
            gamma_t[m] = np.sum(np.multiply(C, gamma_km))
            #print gamma_t[m]

            # end loop through K clusters
        # end loop through SNPs

        alpha = np.sum(C_t, axis=0)
        p_t = st.dirichlet.rvs(alpha)
        p_t = p_t.ravel()
        p_list.append(p_t)

        logging.info("Iteration %d:" % i)
        print p_t

    # end loop iterations

    BURN = its/4
    p_est = np.mean(p_list[BURN:], axis=0)

    return p_est


def precompute_terms(W, beta_tilde, name, outdir):
    A = np.matmul(np.transpose(W), W)
    M = W.shape[0]
    psi = np.empty((M,1))
    for m in range(0,M):
        psi[m] = np.matmul(np.transpose(beta_tilde), W[:,m])

    # save for future use
    A_file = os.path.join(outdir, name+'.A')
    psi_file = os.path.join(outdir, name+'.psi')
    np.save(A_file, A)
    np.save(psi_file, psi)
    return A, psi


def initialize_p(K):
    # random draw to initialize values
    p_init = np.random.dirichlet([1]*K,1)
    p_init = p_init.ravel()
    #p_init = [.5, .5]
    return p_init


def initialize_C_gamma(p_init, mu_vec, sigma_vec, M):
    # create empty array to hold values
    K = len(mu_vec)
    gamma_init = np.empty((M,K))

    C_init = np.random.multinomial(n=1, pvals=p_init, size=M)

    for k in range(0, K):
        gamma_init[:, k] = st.norm.rvs(mu_vec[k], sigma_vec[k], size=M)

    gamma_init = np.sum(np.multiply(gamma_init, C_init), axis=1)

    return C_init, gamma_init


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--mu_vec", dest="mu_vec",)
    parser.add_option("--sigma_vec", dest="sigma_vec")
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

    if options.bins is None:
        mu_vec = [float(item) for item in options.mu_vec.split(',')]
        sigma_vec = [float(item) for item in options.sigma_vec.split(',')]
    else:
        bins = int(bins)
        # TODO

    logging.info("Reading in gwas file: %s" % gwas_file)
    df = pd.read_csv(gwas_file, sep=' ')
    beta_tilde = np.asarray(df['BETA_STD_I'])

    # set seed
    np.random.seed(seed)

    ld_half_file = options.ld_half_file
    W = np.loadtxt(ld_half_file)
    logging.info("Using ld half file: %s" % ld_half_file)

    # precompute terms or load in
    if precompute == 'y':
        logging.info("Pre-computing terms")
        A, psi = precompute_terms(W, beta_tilde, name, outdir)
    else: # assumes files have already been made
        # check if files exist
        A_file = os.path.join(outdir, name+'.A.npy')
        psi_file = os.path.join(outdir, name+'.psi.npy')
        try:
            logging.info("Reading in pre-computed matrices")
            A = np.load(A_file)
            psi = np.load(psi_file)
        except:
            logging.info("ERROR: could not find pre-computed matrices...re-calculating them which may take awhile")
            A, psi = precompute_terms(W, beta_tilde, name, outdir)


    # get meta values from data
    K = len(mu_vec)
    M = beta_tilde.shape[0]

    # intialize values of chain
    logging.info("initializing p")
    p_init = initialize_p(K)
    print(p_init)

    logging.info("initializing gamma_k")
    logging.info("initializing C_k")
    C_init, gamma_init = initialize_C_gamma(p_init, mu_vec, sigma_vec, M)

    p_est = gibbs(p_init, gamma_init, C_init, mu_vec, sigma_vec, W, A, psi, beta_tilde, N, its)
    print "Estimate: "
    print p_est

if __name__== "__main__":
  main()
