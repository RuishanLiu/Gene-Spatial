import os, sys, pickle, argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as spatial
from scipy import stats, integrate
from tabulate import tabulate
from scipy import stats, integrate
csv.field_size_limit(sys.maxsize)

from functions.load_rawdata import *
from functions.process_spatdata import *
from functions.model_spatdata import *

NUM = 502628 # number of individuals
TOLENGTH = 0.227 # pixel to length (um)
genes_sign = ['SST', 'INS', 'GLUC']  # Genes used as marker for endocrine

def parse_args():
    parser = argparse.ArgumentParser(description='teacher-student for tabular')
    parser.add_argument('--path_file', help='Path to gene expression csv file', type=str)  
    parser.add_argument('--path_nuc', help='Path to nuclei position csv file', type=str)
    parser.add_argument('--path_gene', help='Path to gene conversion csv file', type=str)
    parser.add_argument('--path_save', help='Path to save results', type=str) 
    args = parser.parse_args()
    return args, parser

def main():   
    ##################################### Read Data ######################################
    # class SpatData defined in functions.model_spatdata
    spatdata = SpatData(args.path_file, args.path_gene, args.path_nuc)
    
    ################################ Islets Identification ###############################
    
    threshold = 2 # threshold for islets identification
    spatdata = find_islets(spatdata, genes_sign, threshold=threshold, plot_high=False)
    spatdata = islet_stats(spatdata)
    
    ###################################### Analysis ######################################
    
    # density profile (to islets) grids in pixel - min, max, spacing.
    supports = [0, 7000, 200]
    # density profile (to nuclei) grids in pixel - min, max, spacing.
    supports_nuc = [0, 150, 200]
    spatdata = get_profile(spatdata, supports=supports, supports_nuc=supports_nuc)

    ###################################  Save Results  ###################################
    
    with open(args.path_save + '.pkl', 'wb') as output:
        pickle.dump(spatdata, output, pickle.HIGHEST_PROTOCOL)


if __name__ == __main__:
    args, parser = parse_args()
    main()
        

