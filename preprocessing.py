#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:21:34 2019

@author: silvia
"""

import numpy as np
import pandas as pd
from pathlib2 import Path
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import normalize
#organ_data_normed=no.zeros(organ_data.shape[0],organ_data.shape[1])
#for i in range(organ_data.shape[1](:
#    organ_data_normed[:,i]=organ_counts[:,i]/np.repeat(np.sum(organ_counts[:,i]), organ_counts.shape[1])


def load_organ_raw_counts(path_data, organ):
    path_data_=path_data/organ
        
    organ_data=np.load(str(path_data_), allow_pickle=True)
        
    organ_counts=organ_data[1:][:,1:]
    organ_counts=organ_counts.astype('float32')    
    organ_name=organ.split(sep='.')[0]
    organ_name=organ_name.split(sep='-')[0]
    
    n_cells_org =organ_counts.shape[1]
    n_genes=organ_counts.shape[0] 
    
    return organ_counts, organ_name, n_cells_org, n_genes

def save_normed(path_data, organ_name, path_save_normed):
    organ_name=organ_name
    organ=organ_name+'-counts.npy'
    path_save_normed=path_save_normed/'norm_{}.npy'.format(organ_name)
    if not path_save_normed.exists():
        organ_counts, organ_name, n_cells_org, n_genes=load_organ_raw_counts(path_data, organ)
        normed_data= normalize(organ_counts, norm='l1', axis=0)
        np.save(str(path_save_normed), normed_data)

def preprocessed_matrices(path_data, organ_name, path_save_preprocessed, thresh=700):            
    organ=organ_name+'-counts.npy'
    path_save_preprocessed=path_save_preprocessed/'prepro_{}.npy'.format(organ_name)
    organ_counts, organ_name, n_cells_org, n_genes=load_organ_raw_counts(path_data, organ)
    active_genes=np.zeros((n_cells_org))
    for i in range(n_cells_org):
        active_genes[i]=np.argwhere(organ_counts[:,i]>0).shape[0]
    idx=np.argwhere(active_genes>thresh)[0]
    preprocessed_matrix=organ_counts[:,idx]
    preprocessed_matrix= normalize(preprocessed_matrix, norm='l1', axis=0)
    np.save(path_save_preprocessed,preprocessed_matrix)
    
    
    
def QC(path_data, organ_name, path_save):
    organ=organ_name+'-counts.npy'
    organ_counts, organ_name, n_cells_org, n_genes=load_organ_raw_counts(path_data, organ)
    # vediamo quali cellule hanno un alto valore di mean counts 
    #high_counts_cell = np.argwhere(np.mean(organ_counts, axis=1)>)
    cells_mean_value = np.mean(organ_counts, axis=0)
    Max,Min = np.max(cells_mean_value), np.min(cells_mean_value)
    argmax, argmin =np.argwhere(cells_mean_value==Max), np.argwhere(cells_mean_value==Min)
    print(Max,Min)
    print(argmax[0][0], argmin)
    cells_mean_value_normed=np.zeros((cells_mean_value.shape))
    print(cells_mean_value.shape)
    for i in range(organ_counts.shape[1]):
        cells_mean_value_normed[i]=cells_mean_value[i]/np.sum(organ_counts[:,i],axis=0)
    print(cells_mean_value_normed.shape)
    low_expression = np.argwhere(cells_mean_value<10**(-1))
    print(np.argwhere(cells_mean_value<10**(-1)))
    plt.figure(figsize=(10,7))
    for i in low_expression:
        print(i[0])
        n_genes=np.argwhere(organ_counts[:,i[0]]>0).shape[0]
        plt.plot(organ_counts[:,i[0]], 'o', markersize=3.5, label='Cell {}. Expressed genes: {}'.format(i[0],n_genes))  
        plt.yscale('log') 
    plt.xlabel('Genes')
    plt.ylabel('Counts')
    plt.title('{} Low Gene Expression Control'.format(organ_name))
    if low_expression.shape[0]< 10:
        plt.legend()
    plt.savefig(str(path_save/'{}.png'.format(organ_name)))
    plt.close()
    return low_expression[0]

def save_mean_var(path_data, organ_name, path_save_mean, path_save_var):
    organ_name=organ_name
    organ=organ_name+'-counts.npy'
    path_save_mean=path_save_mean/'mean_{}.npy'.format(organ_name)
    path_save_var=path_save_var/'var_{}.npy'.format(organ_name)
    if not path_save_var.exists():
        organ_counts, organ_name, n_cells_org, n_genes=load_organ_raw_counts(path_data, organ)
        cells_mean_value = np.mean(organ_counts, axis=0)
        np.save(str(path_save_mean), cells_mean_value)
        print('done')
        cells_var=np.var(organ_counts, axis=0)  
        print('done')
        np.save(str(path_save_var), cells_var)
        del organ_counts, organ_name, n_cells_org, n_genes, cells_var
        
def plot_mean_variance(path_data, organ_name, path_save_mean, path_save_var, path_save_mean_var):
    organ_name=organ_name
    organ=organ_name+'-counts.npy'
    path_save_mean=path_save_mean/'mean_{}.npy'.format(organ_name)
    path_save_var=path_save_var/'var_{}.npy'.format(organ_name)
    cells_mean_value, cells_var = np.load(path_save_mean, allow_pickle=True), np.load(path_save_var,  allow_pickle=True)   
    print(cells_var[:10])
    plt.figure(figsize=(15,7))
    plt.errorbar(np.arange(0,cells_mean_value.shape[0], 1), cells_mean_value,yerr=(cells_var), marker='o', mfc='red', mec='green') 
    plt.yscale('log')

    plt.xlabel('Cells')
    plt.ylabel('Mean')
    plt.title('{} cells: Mean and Variance'.format(organ_name))
    plt.savefig(str(path_save_mean_var/'{}.png'.format(organ_name)))
    plt.close()


def plot_mean(path_save_mean, path_save_mean_var, path_save_low, organ_name):
    
    mean_data=np.load(str(path_save_mean/'mean_{}.npy'.format(organ_name)), allow_pickle=True)
    
    low_cells=np.argwhere(mean_data<10**(-1))
    path_save_low=path_save_low/'cells_data'
    if low_cells.shape[0]>0:
        np.save(str(path_save_low/'low_cells_idx_{}'.format(organ_name)),low_cells[0])
    
    plt.figure(figsize=(10,7))
    plt.plot(mean_data, 'ro', markersize=2)
    plt.plot(np.repeat(10**(-1), mean_data.shape[0]))
    plt.yscale('log')
    plt.title('Mean gene expression values for {}'.format(organ_name))
    plt.xlabel('cell')
    plt.ylabel('mean counts')
    plt.savefig(str(path_save_mean_var/'meancounts_{}.png'.format(organ_name)))
    plt.close()     
  
def cells_type(path_data, path_annot, path_data_csv, organ_name, path_save_cell_type):  
    organ_csv=organ_name+'-counts.csv'
    organ=organ_name+'-counts.npy'
    path_save=path_save_cell_type/'cell_types_{}.npy'.format(organ_name)
    path_data_csv=path_data_csv/organ_csv   
    organ_data_csv = pd.read_csv(path_data_csv, sep=',',header=None)
        
    cells_id=np.array(organ_data_csv.iloc[0,1:]).astype('str')
    print(cells_id.shape)
    cells_type_organ=np.zeros((cells_id.shape))
    cells_type_organ=cells_type_organ.astype('str')
    annot = pd.read_csv(path_annot, sep=',', header=(0))
    cells_name=np.array(annot['cell']).astype('str')
    cells_types=np.array(annot['cell_ontology_class']).astype('str')
    #print(cells_id)
    
    for i,cell in enumerate(cells_id):
        for j in range(cells_name.shape[0]):
            if cells_id[i]==cells_name[j]:
                cells_type_organ[i]=cells_types[j]
    np.save(path_save,cells_types)
    
    
def main():
    path_data=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/raw_matrices_FACS')
    path_csv=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/00_facs_raw_data/FACS/FACS')
    path_annot=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/00_facs_raw_data/annotations_FACS.csv')
    path_save_low=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/preprocessing/low_gene_expr_control')
    path_save_mean_var=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/preprocessing/mean_variance')
    path_save_mean=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/mean_var/mean')
    path_save_cell_type=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/cells_type')
    path_save_var=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/mean_var/variance')
    path_save_normed=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/normed_matrices')
    path_save_preprocessed=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/preprocessed_matrices')
    organ_name='Aorta'
    organ=organ_name+'-counts.npy'
     # path_save_cell_type, 
    
#    path_save_low_data=path_save_low/'cells_data'
#    low_cells_idx=np.load(str(path_save_low_data/'low_cells_idx_{}.npy'.format(organ_name)), allow_pickle=True)
#    print(low_cells_idx)
#    plt.plot(organ_counts[:,111],'ro', markersize=2.5)
#    plt.yscale('log')
#    plt.ylim(top=1.4)
    #print(np.argwhere(organ_counts[:,111]==))
    listdir_ = os.listdir(path_data)
    #listdir_=listdir_.remove('Large_Intestine-counts.npy')
    for i,organ in enumerate(listdir_):
        if i != 18:
            organ_name=organ.split(sep='.')[0]
            organ_name=organ_name.split(sep='-')[0]
            print(organ_name, i)
            cells_type(path_data, path_annot, path_csv,organ_name, path_save_cell_type)
#            preprocessed_matrices(path_data, organ_name, path_save_preprocessed, thresh=700)
        #plot_mean(path_save_mean, path_save_mean_var,path_save_low, organ_name)
        #plot_mean_variance(path_data, organ_name, path_save_mean, path_save_var,path_save_mean_var)
        
        #low_expression=QC(path_data, organ_name, path_save)
        #np.save(str(path_save/'low_cells_{}.npy'.format(organ_name)), low_expression)
    
if __name__=='__main__':
    main()

