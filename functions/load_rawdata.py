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

NUM = 502628 # number of individuals
TOLENGTH = 0.227 # pixel to length (um)

'''Functions for Loading Raw Data'''

def save_data(dir_save, names, data):
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    if type(names) is list:
        for i, name in enumerate(names):
            np.save(os.path.join(dir_save, name), data[i])
    else:
        np.save(os.path.join(dir_save, names), data)
        

def get_items(dir_file):
    # Get the first row in csv file dir_file
    with open(dir_file, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', dialect=csv.excel_tab)
        for i, row in enumerate(reader):
            row_item = row
            break
    return row_item
    

def get_ind(dir_file, names):
    # Find the list index of one item represented by name in csv file dir_file
    # return the index if names is a char; a list of index if names is a list
    row_item = get_items(dir_file)
    if type(names) is list:
        inds = []
        for name in names:
            inds.append(row_item.index(name))
        return inds
    else:
        return row_item.index(names)
    
def get_data(dir_file, names):
    # Get the data of item names for all individuals
    # return the data list if names is char; a list containing data lists if names is a list
    
    if type(names) is list:
        data = [[] for i in names]
    else:
        data = []
        
    inds = get_ind(dir_file, names)
        
    with open(dir_file, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            #if np.mod(i, 10000) == 0:
                #print(i)
            if type(names) is list:
                for j, ind in enumerate(inds):
                    data[j].append(row[ind])
            else:
                data.append(row[inds])
    return data

def generate_data(dir_file, dir_save, names):
    names_new = []
    # only generate data which are not generated before
    for name in names:
        if not os.path.isfile(os.path.join(dir_save, name) + '.npy'):
            names_new.append(name)
            
    # if all items are generated before
    if not names_new:
        return
            
    data_new = get_data(dir_file, names_new)
    save_data(dir_save, names_new, data_new)
    

def get_batch(dir_file, nums):
    # return randomly-chosen num individual batch
    
    if type(nums) is list:
        number_rand = []
        data = []
        k = []
        for num in nums:
            number_rand.append(np.sort(np.random.choice(np.arange(1, NUM+1), num, replace=False)))
            data.append([])
            k.append(0)
    else:
        k = 0
        number_rand = np.sort(np.random.choice(np.arange(1, NUM+1), nums, replace=False))
        
    with open(dir_file, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            if type(nums) is list:
                for j, num in enumerate(nums):
                    if i == number_rand[j][k[j]]:
                        if k[j] < num - 1:
                            k[j] = k[j] + 1
                        data[j].append(row)
            else:
                if i == number_rand[k]:
                    k = k + 1
                    data.append(row)
                
            if np.mod(i, 10000) == 0:
                print(i)   
            if i == 100:
                break
    return data
    
def generate_batch(dir_file, dir_save, nums):
    # generate and save random batch
    nums_new = []
    nums_new_str = []
    for num in nums:
        name = str(num)
        if not os.path.isfile(os.path.join(dir_save, name) + '.npy'):
            nums_new.append(num)
            nums_new_str.append(str(num))           
    data = get_batch(dir_file, nums_new)
    save_data(dir_save, nums_new_str, data)
    