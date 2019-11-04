#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:47:30 2019

@author: silvia
"""

import numpy as np
import pandas as pd
from pathlib2 import Path
import os
import matplotlib.pyplot as plt

def U_plots(path_data, path_save,par):
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    for organ in listdir_:
        
        path_data_=path_data/organ
        
        organ_data=np.load(str(path_data_), allow_pickle=True)
        
        organ_counts=organ_data[1:][:,1:]
        
        n_cells=organ_counts.shape[1]
        n_genes=organ_counts.shape[0]
        
        n_cells_forgenes=np.zeros(n_genes)
        for i in range(n_genes):
            n_cells_forgenes[i]=np.argwhere(organ_counts[i,:]>par).shape[0]
       
        # U plot
        organ_name=organ.split(sep='.')[0]
        plt.figure(figsize=(5,3))
        plt.hist(n_cells_forgenes/np.repeat(n_cells,n_genes), label=organ_name)
    
        plt.xlabel('k')
        plt.ylabel('P(k)')
        plt.legend(loc='upper right')
        plt.title('Probability of gene occurrence in k cells (Threshold: {})'.format(par))
        plt.savefig(str(path_save/'{}.png'.format(organ_name)))
        plt.close()
        
def Rank_plots(path_data, path_save):
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    for organ in listdir_:
        
        path_data_=path_data/organ
        
        organ_data=np.load(str(path_data_), allow_pickle=True)
        
        organ_counts=organ_data[1:][:,1:]
        
        n_cells=organ_counts.shape[1]
        n_genes=organ_counts.shape[0]
        
        plt.plot(np.sort(organ_counts[1,:]), 'o',markersize=1, label='1')
        mean_gene=np.mean(organ_counts, axis=1)
        
        organ_name=organ.split(sep='.')[0]
        
        plt.figure(figsize=(8,5))
        plt.subplot(1,2,1)
        plt.plot(np.sort(mean_gene),'o',markersize=1, label='{}'.format(organ_name))
        plt.legend()
        plt.xlabel('gene')
        plt.ylabel('Mean expression')
        plt.yscale('log')
        
        plt.subplot(1,2,2)
        plt.hist(np.sort(mean_gene), label='{}'.format(organ_name))
        plt.legend()
        
        plt.title('Rank plot and density')
        plt.savefig(str(path_save/'{}.png'.format(organ_name)))
        plt.close()
        
def main():
    path_data=Path('/home/silvia/Desktop/BORSA/SPECIAL_ISSUE/mypyhon/data/raw_matrices_FACS')
    path_save=Path('/home/silvia/Desktop/BORSA/SPECIAL_ISSUE/mypyhon/data/stat_FACS/U')
    path_save_rank=Path('/home/silvia/Desktop/BORSA/SPECIAL_ISSUE/mypyhon/data/stat_FACS/U')
#    U_plots(path_data, path_save,0)
#    U_plots(path_data, path_save,100)
#    U_plots(path_data, path_save,200)
    Rank_plots(path_data, path_save_rank)
    
    
if __name__=='__main__':
    main()