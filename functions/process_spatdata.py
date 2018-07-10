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
from functions.model_spatdata import *

NUM = 502628 # number of individuals
TOLENGTH = 0.227 # pixel to length (um)

'''Functions to process spatdata'''

def find_islets(dat, genes_sign, r_min=50, r_max=500, dr=50, 
                threshold=3, plot_high=False):
    # islets identification       
    def get_nuc_ind(data_cell, number_gene, threshold):
        ind = []
        for i, t in enumerate(data_cell[number_gene]):
            if t >= threshold: # include threshold
                ind.append(i)
        return np.array(ind)

    def plot_nuc(data, title='', window=False, xmin=-5000, xmax=20000, ymin=0, ymax=35000, size=1):
        font_size = 13   
        fig = plt.figure(figsize=(5, 5))
        plt.scatter(data[0], data[1], s = size, edgecolors = 'face', c = 'k', alpha = 1)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)  
        plt.title(title)
        if window:
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

    def get_inds_high_express(inds_nuc_chosen, sign_genes):
        inds_high = []
        sign_window = True
        for sign_gene in sign_genes:
            for ind in inds_nuc_chosen[sign_gene]:
                if ind not in inds_high:
                    if sign_window:
                        x, y = float(dat.data_nuc[2][ind]), float(dat.data_nuc[3][ind])
                        if False: # x<12800 or x>15000 or y<6500 or y>11000:
                            pass
                        else:
                            inds_high.append(ind)
                    else:
                        inds_high.append(ind)
        return np.array(inds_high)           

    ### choose cells have expressions of genes_sign above threshold to identify islets
    inds_nuc_chosen = []
    poses_nuc_chosen = []
    for i in range(len(dat.name_gene)):
        ind_nuc_chosen = get_nuc_ind(dat.data_cell, number_gene=i, threshold=threshold) # include threshold
        if ind_nuc_chosen.shape[0] > 0:
            pos_nuc_chosen = dat.data_nuc[2:4, ind_nuc_chosen]
            if plot_high:
                plot_nuc(pos_nuc_chosen, title=dat.name_gene[i], size=1, window=True, xmin=10000, xmax=18000, ymin=6000, ymax=12000)
        else:
            pos_nuc_chosen = np.array([])
        inds_nuc_chosen.append(ind_nuc_chosen)
        poses_nuc_chosen.append(np.array(pos_nuc_chosen, dtype=np.float))
    inds_nuc_chosen = np.array(inds_nuc_chosen)
    poses_nuc_chosen = np.array(poses_nuc_chosen)      

    ind_sign = get_ind_names(dat, genes_sign)
    inds_high = get_inds_high_express(inds_nuc_chosen, ind_sign)
    print('Number of cells chosen for signed genes:', len(inds_high))


    ### Looking for islets (including smaller ones)
    # Find centers - inds_centers: index in inds_high
    def get_neighbor(data, radius):
        tree = spatial.KDTree(np.array(data))
        neighbors = tree.query_ball_tree(tree, radius)
        return np.array(neighbors)

    def get_pos_nuc_all(data_nuc):
        return np.array([(float(data_nuc[2][i]), float(data_nuc[3][i])) for i in range(data_nuc.shape[1])])        

    def subtract(a, b): # returns: list a - all elements in list b
        temp = []
        for i in a:
            if i not in b:
                temp.append(i)
        return temp

    def get_centers(neighbors, inds_exist, thresh=False):
        num_neighbor = np.array([float(len(neighbor)) for neighbor in neighbors])
        threshold = thresh if thresh else np.mean(num_neighbor)
        inds_centers = []
        temp_exist = inds_exist
        for i in range(len(inds_high)):
            ind_i = np.argmax(num_neighbor[temp_exist])
            ind_i = temp_exist[ind_i]
            if num_neighbor[ind_i] < threshold:
                break
            if subtract(neighbors[ind_i], inds_exist) != []:
                temp_exist = subtract(temp_exist, [ind_i]) 
                continue
            inds_centers.append(ind_i)
            inds_exist = subtract(inds_exist, neighbors[ind_i])   
            temp_exist = subtract(temp_exist, neighbors[ind_i]) 
            if len(temp_exist) == 0:
                break       
        return inds_centers, inds_exist    

    radius = range(r_min, r_max+dr, dr)[::-1]
    single_radius = dr
    neighbors_radius = [get_neighbor(get_pos_nuc_all(dat.data_nuc[:, inds_high]), radius=r) for r in radius]  

    inds_exist = range(len(inds_high))
    inds_circle = []
    radius_circle = []
    for i, r in enumerate(radius):
        neighbors = neighbors_radius[i]
        if i == len(radius)-1:
            inds_centers, inds_exist = get_centers(neighbors, inds_exist, thresh=2)
        else:
            inds_centers, inds_exist = get_centers(neighbors, inds_exist)
        inds_circle = inds_circle + inds_centers
        radius_circle = radius_circle + [r for i in inds_centers] 

    circles = np.array([[float(dat.data_nuc[2][inds_high[center]]), float(dat.data_nuc[3][inds_high[center]]), 
                         radius_circle[i]] for i, center in enumerate(inds_circle)])  

    # Refinement - the position of nucleis are not continuos and the circles may overlap
    radius_set = radius + [max(radius)+r for r in radius]
    neighbors_center = {r: get_neighbor(circles[:, :2], radius=r) for r in radius_set}
    num_neighbors_center = []
    circles_refined = circles
    for i in range(circles.shape[0]):
        r = circles[i, 2]
        for delta_r in radius:
            neighbor = neighbors_center[r+delta_r][i]
            for j in neighbor:
                if j!=i:
                    r_j = circles[j, 2]
                    if r_j > delta_r:
                        circles_refined[i, 2] -= (r_j-delta_r)/2    

    dat.circles = circles_refined
    return dat

def islet_stats(dat, d_densit=200):
    def get_dist(circles, data):
        # return the dist to the circle center
        num_data = dat.data.shape[1]
        num_circles = dat.circles.shape[0]
        dist = np.zeros([num_data, num_circles])
        genes_in = np.zeros(num_data)
        dist_boundary = np.zeros_like(dist)
        for i in range(num_data):
            for j in range(num_circles):
                dist[i, j] = np.linalg.norm(dat.circles[j, :2] - np.array(dat.data[2:4, i], dtype=np.float))
                dist_boundary[i, j] = dist[i, j] - dat.circles[j, 2]
                if dist[i, j] < circles[j, 2]:
                    genes_in[i] = 1
        return dist, genes_in, dist_boundary

    def get_index_in_out(index, genes_in):
        index_out = [[] for i in range(len(index))]
        index_in = [[] for i in range(len(index))]
        for i, inds in enumerate(index):
            for ind in inds:
                if genes_in[ind] == 0:
                    index_out[i].append(ind)
                else:
                    index_in[i].append(ind)
        return index_in, index_out

    def get_distance_distribution(dist_boundary, index_out):
        distri_dist = []
        for i, inds in enumerate(index_out):
            distri_dist.append(np.min(dist_boundary[inds], axis=1))
        return distri_dist 
    
    def get_nuclei_distribution(dat):
        ind_cell = np.array(dat.data[4, :], dtype=np.int64)  # parent cell data
        loc_gene = np.array([dat.data[2:4, i] for i in range(dat.data.shape[1])], dtype=np.float)
        dist_nuc = np.linalg.norm(loc_gene - dat.pos_nuc[:, ind_cell-1].T, axis=1)
        distri_dist = [dist_nuc[ind] for ind in dat.index]
        return distri_dist

    dat.distri_dist_nuc = get_nuclei_distribution(dat)
    dat.dist_center, dat.genes_in, dat.dist_boundary = get_dist(dat.circles, dat.data)
    dat.index_in, dat.index_out = get_index_in_out(dat.index, dat.genes_in)
    dat.distri_dist = get_distance_distribution(dat.dist_boundary, dat.index_out)
    return dat


def get_density(Xs, support):
    # http://pythonhosted.org/PyQt-Fit/KDE_tut.html
    from scipy import stats
    from pyqt_fit import kde, kde_methods  #easy_install distribute; sudo pip install git+https://github.com/Multiplicom/pyqt-fit.git
    densities = []
    for X in Xs:
        if X.shape[0] < 2:
            density = 1. * support
        else:
            est_lin = kde.KDE1D(X, lower=0, method=kde_methods.linear_combination) 
            density = est_lin(support)
           # density[density < 0] = 0
        densities.append(density) 
    return densities

def get_density_kernel(Xs, support):
    from scipy import stats, integrate
    densities = []
    for X in Xs:
        bandwidth = 1.06 * X.std() * X.size ** (-1 / 5.)

        kernels = []
        for x_i in X:
            kernel = stats.norm(x_i, bandwidth).pdf(support)
            kernels.append(kernel)

        density = np.sum(kernels, axis=0)
        density /= integrate.trapz(density, support)
        densities.append(density)
    return densities

def get_corr(dat, densities):
    n = len(dat.name_gene)
    entropy_den = np.zeros([n, n])
    for i in range(n):
        for j in range(n):
            entropy_den[i, j] = stats.entropy(densities[i], qk=densities[j]) 
    return entropy_den

def get_profile(dat, supports=[100, 7000, 200], supports_nuc=[1, 150, 200]):
    support = np.linspace(supports[0], supports[1], supports[2])
    Xs = [dat.distri_dist[ind] for ind in range(len(dat.name_gene))]   
    densities = get_density(Xs, support)
    dat.support = support
    dat.densities = densities    
    dat.corr_den = get_corr(dat, densities)
    
    support_nuc = np.linspace(supports_nuc[0], supports_nuc[1], supports_nuc[2])
    Xs_nuc = [dat.distri_dist_nuc[ind][dat.distri_dist_nuc[ind]<500] for ind in range(len(dat.name_gene))]  # remove outliers
    densities_nuc = get_density(Xs_nuc, support_nuc)
    dat.support_nuc = support_nuc
    dat.densities_nuc = densities_nuc    
    dat.corr_den_nuc = get_corr(dat, densities_nuc)    
    return dat