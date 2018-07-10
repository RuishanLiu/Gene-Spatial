import numpy as np
import os
import matplotlib.pyplot as plt
import csv
import sys
csv.field_size_limit(sys.maxsize)
import tabulate
# https://github.com/gregbanks/python-tabulate
import scipy.spatial as spatial
from scipy import stats, integrate
from tabulate import tabulate
from scipy import stats, integrate
import pickle

from functions.load_rawdata import *
from functions.process_spatdata import *

NUM = 502628 # number of individuals
TOLENGTH = 0.227 # pixel to length (um)

'''SpatData Class Definition'''

class SpatData:
    def __init__(self, dir_file, dir_filegene, dir_file_nuc):
        ### Import Data
        names = get_items(dir_file)
        print(names)

        data_origin = get_data(dir_file, names)
        data_origin = np.array(data_origin)
        if data_origin[0][0] == '':
            data_origin = data_origin[1:]

        ### Preprocess Data
        # throw away low quality data point
        threshold = 0.55
        prob = 1.

        qual = np.array(data_origin[6][1:], dtype = np.float64)
        ind_corr = []
        for i, q in enumerate(qual):
            if q > threshold:
                if np.random.rand() < prob:
                    ind_corr.append(i+1)
        data = np.array([datum[ind_corr] for datum in data_origin])
        print('Original data number:', data_origin.shape[1])
        print('Data number after processing:', data.shape[1])

        # Get the represented number of correct genes
        names_gene = get_items(dir_filegene)
        data_gene = get_data(dir_filegene, names_gene)
        num_gene = data_gene[2][1:]
        num_gene = np.array(num_gene, dtype = np.int)
        num_gene = np.ndarray.tolist(num_gene)
        name_gene = data_gene[0][1:]
        print(num_gene)
        print(name_gene)

        # index[i]: the data index of all gene i in data

        index = [[] for i in range(len(num_gene))]  # includes the title row

        for i, datum in enumerate(data[5]):
            if i==0:
                continue
            temp_num = int(datum)
            try: # may be not in the list
                ind = num_gene.index(temp_num)
                index[ind].append(i)
            except:
                pass

        # The data number of each gene
        num_ind = []
        for ind in index:
            num_ind.append(len(ind))
        print('data number of each gene:')
        print(num_ind)

        temp = []
        for i, ind in enumerate(num_ind): 
            if ind > 1000:
                temp.append(i)
        print('', temp)
        print(np.array(name_gene)[temp])
        print(num_ind[9])

        ### Nuclei
        names_nuc = get_items(dir_file_nuc)
        print(names_nuc)
        data_nuc_origin = get_data(dir_file_nuc, names_nuc)   
        if data_nuc_origin[0][0] == '':
            data_nuc_origin = data_nuc_origin[1:]
        data_nuc_origin = np.array(data_nuc_origin)
        data_nuc = data_nuc_origin[:, 1:]
       # data_nuc = data_nuc[1:, :]
        print(data_nuc.shape)
        def get_pos_nuc_all(data_nuc):
            return np.array([(float(data_nuc[2][i]), float(data_nuc[3][i])) for i in range(data_nuc.shape[1])])
        point_nuc = get_pos_nuc_all(data_nuc)          
        
        # data_cell[i, j]: the number ithe gene in jth (numbered j+1) parent cells;
        ind_cell = np.array(data[4, :], dtype=np.int64)  # parent cell data
        max_nuc = max(np.array(data_nuc[1], dtype=np.int)) # max(ind_cell)

        data_cell = np.zeros([len(num_gene), max_nuc])
        for i, ind in enumerate(index):
            if ind == []:
                continue
            for j in (ind_cell[np.array(ind)]):
                data_cell[i, j-1] = data_cell[i, j-1] + 1 
                
        ## pos_nuc[0 (1), j] = The x (y) coordinate of jth (numbered as j+1) nuclei
        pos_nuc = np.zeros([2, max_nuc])
        for i in range(data_nuc.shape[1]):
            ind = int(data_nuc[1][i])
            pos_nuc[0][ind-1] = float(data_nuc[2][i])
            pos_nuc[1][ind-1] = float(data_nuc[3][i])         

        self.data = data # all the data after processing
        self.index = index # index[i]: the data index of all gene i in data
        self.name_gene = name_gene # List of all the genes' names ['CBLN1', 'MASP1', 'GJA4', 'MAL2', ...]
        self.data_nuc = data_nuc # 
        self.data_cell = data_cell # data_cell[i, j]: the number ithe gene in jth (numbered j+1) parent cells;
        self.pos_nuc = pos_nuc # pos_nuc[0 (1), j] = The x (y) coordinate of jth (numbered as j+1) nuclei
        self.point_nuc = point_nuc 
        
def get_ind_names(dat, names):
    # get the index for the genes with name names
    return [dat.name_gene.index(name) for name in names]  