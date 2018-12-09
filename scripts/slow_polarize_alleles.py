#!/usr/bin/env python

import os
import numpy as np 
import re

def print_list(line, f):
    print line 
    for i,c in enumerate(line):
        f.write(c)
        if i != len(line)-1:
            f.write(" ")
        else:
            f.write("\n")

mrca_dir="/u/home/r/ruthjohn/ruthjohn/BayesPred/unity_mixture/misc/chimp_alleles"

gwas_file="/u/home/r/ruthjohn/ruthjohn/UNITY_analyses/data/clean_gwas/BIP_2012.txt" 

# open gwas file 
gwas_f = open(gwas_file, 'r')

# read heading 
header = gwas_f.readline()
header = header.strip().split()
header.append('AA')

chr_ind = header.index('CHR')
A1_ind = header.index('A1')
bp_ind = header.index('BP')
beta_std_ind = header.index('BETA_STD')

# print heading of outfile
out_file = "test_polarize_file.txt"
out_f = open(out_file, "w")
print_list(header, out_f)

for chr in range(1,3): # loop through chr1 and chr2
    
    # get mrca chr file
    mrca_file =  os.path.join(mrca_dir, "chr%i.txt" % chr)
    mrca_f = open(mrca_file, 'r')

    # loop through each file 

    # get elements from gwas file 
    #gwas_line = gwas_f.readline() # first line after header 
    for gwas_line in gwas_f:
        gwas_line = gwas_line.strip().split()
        gwas_chr = gwas_line[chr_ind]
        gwas_bp = gwas_line[bp_ind]
        gwas_A1 = gwas_line[A1_ind]

        # debugging for chr 1 and 2 only!
        if gwas_chr == '2':
            exit(1)

        # get elements from mrca file
        mrca_line = mrca_f.readline() # no header: chri BP AA
        mrca_line = mrca_line.strip().split()

        mrca_chr = re.sub("[^0-9]", "", mrca_line[0])
        mrca_bp = mrca_line[1]
        mrca_AA = mrca_line[2].upper() # make sure comparing upper-case 

        # check if AA not available 
        print mrca_chr, gwas_chr, mrca_bp, gwas_bp
        if mrca_AA != '-' and mrca_chr == gwas_chr and mrca_bp == gwas_bp:
            if mrca_AA == gwas_A1:
                pass
            else:
                gwas_line[beta_std_ind] = gwas_line[beta_std_ind] * -1 
            # print line 
        
            # add AA to line 
            gwas_line.append(mrca_AA)
            print_list(gwas_line, out_f)

        else: # AA not available 
            # move through files
            pass 
        
out_f.close()
