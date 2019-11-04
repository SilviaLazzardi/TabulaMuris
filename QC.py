#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 10:21:53 2019

@author: silvia
"""

import numpy as np
import pandas as pd
from pathlib2 import Path
import os
import matplotlib.pyplot as plt

def plot_qc(path_data, path_save):
    listdir_ = os.listdir(path_data)
    for organ in listdir_:
        
        path_data_=path_data/organ
        print(path_data_)
        
        organ_data=np.load(str(path_data_), allow_pickle=True)
        
        organ_counts=organ_data[1:][:,1:]
        
        # COUNT DEPHT
        organ_counts_cells=np.sum(organ_counts,axis=0)
        # HOW MANY DIFFERENT GENES
        organ_counts_genes=np.zeros((organ_counts.shape[1]))
        
        for i in range(organ_counts.shape[1]):
            organ_counts_genes[i]=np.argwhere(organ_counts[:,i]!=0).shape[0]
        
        count_depth_genes=np.zeros((organ_counts_cells.shape[0],2))
        for i in range(organ_counts_cells.shape[0]):
            count_depth_genes[i,0]=organ_counts_cells[i]
            count_depth_genes[i,1]=organ_counts_genes[i]
        
        idx=np.argsort(count_depth_genes[:,0])[::-1]
        count_depth_genes=count_depth_genes[idx]
        
        organ_name=organ.split(sep='.')[0]
        
        plt.figure(figsize=(10,10))
        title='{} QC'.format(organ_name)
        plt.title(title)
        
        plt.subplot(2,2,1)
        
        plt.hist(organ_counts_cells, bins=np.unique(organ_counts_cells).shape[0])
        plt.xlim(left=-10000, right=np.max(organ_counts_cells)+20)
        plt.ylabel('Frequency')
        plt.xlabel('Count depth')
        
        plt.subplot(2,2,2)
        plt.hist(organ_counts_genes, bins=np.unique(organ_counts_genes).shape[0])
        plt.xlim(left=0, right=np.max(organ_counts_genes)+20)
        plt.ylabel('Frequency')
        plt.xlabel('Number of genes')        

        plt.subplot(2,2,3)
        plt.semilogy(np.sort(organ_counts_cells)[::-1],marker='.')
        plt.ylabel('Count depth') 
        plt.xlabel('Barcode Rank')
        
        plt.subplot(2,2,4)
        plt.plot(count_depth_genes[:,1], count_depth_genes[:,0], 'bo', markersize=2)    
        plt.ylabel('Number of genes')
        plt.xlabel('Count depth')       
        
        plt.suptitle(title, fontsize=20)
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=None)
        plt.savefig(str(path_save/'{}.png'.format(organ_name)))
        plt.close()
def main():
#    import argparse
#    parser = argparse.ArgumentParser()
#    parser.add_argument('--path_data')
#    parser.add_argument('--path_save')
#    args = parser.parse_args()

    # VARIABLES     
    #path_data = args.path_data
    #path_save= args.path_save
    path_data=Path('/home/silvia/Desktop/BORSA/SPECIAL_ISSUE/mypyhon/data/raw_matrices_FACS')
    path_save=Path('/home/silvia/Desktop/BORSA/SPECIAL_ISSUE/mypyhon/QC/FACS')
    plot_qc(path_data, path_save)
    
if __name__=='__main__':
    main()