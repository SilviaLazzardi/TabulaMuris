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
from matplotlib.font_manager import FontProperties
from scipy.stats import norm
fontp=FontProperties()
fontp.set_size('small')

def plot_cell_for_organ(path_data, path_save):
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    n_cells = []
    organs=[]
    
    for organ in listdir_:
        
        path_data_=path_data/organ
        
        organ_data=np.load(str(path_data_), allow_pickle=True)
        
        organ_counts=organ_data[1:][:,1:].astype('float32')
        
        organ_name=organ.split(sep='.')[0]
        organ_name=organ_name.split(sep='-')[0]
        n_cells_org =organ_counts.shape[1]
        n_genes=organ_counts.shape[0]  
        
        n_cells.append(n_cells_org)
        organs.append(organ_name)
   
    plt.figure(figsize=(10,8))
    plt.barh(np.arange(len(organs)),n_cells, align='center', alpha=0.5,color=plt.get_cmap('tab20')(np.arange(20)))
    plt.yticks( np.arange(len(organs)), organs)
    plt.xlabel('Number of cells')
    plt.ylabel('Organs')
    plt.title('Number of cells for organ')
    plt.savefig(str(path_save/'cells_for_organ.png'))
    plt.close()

def U_plots(organ_counts,organ_name, path_save,par):
    n_cells=organ_counts.shape[1]
    n_genes=organ_counts.shape[0]
        
    print(organ_name, n_cells,n_genes)
        
    n_cells_forgenes=np.zeros(n_genes)
    for i in range(n_genes):
        n_cells_forgenes[i]=np.argwhere(organ_counts[i,:]>par).shape[0]
    n_cells_forgenes= n_cells_forgenes.astype(int)
    hist = np.histogram(n_cells_forgenes.astype(int), bins=range(0,len(np.unique(n_cells_forgenes)))) 
    print('Max and min occurrences:', np.max(hist[0]), np.min(hist[0]), 'Par: ', par)
    hist_normed=hist[0]/np.repeat(n_cells,hist[0].shape[0])  
    print(hist_normed.shape)
    plt.plot(hist[1][10:][:-1], hist_normed[10:], label=organ_name)
    plt.xlabel('k')
    plt.ylabel('N(k)/tot cells')
        
def main_Rank_plots(path_data, path_save_rank):
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    mean_genes=[]
    organ_names=[]
    for organ in listdir_[:]:
        path_data_=path_data/organ       
        organ_data=np.load(str(path_data_), allow_pickle=True)      
        organ_counts=organ_data[1:][:,1:].astype('float32')   
        organ_name=organ.split(sep='.')[0]
        organ_name=organ_name.split(sep='-')[0]
        genes_num=organ_counts.shape[1]
        organ_names.append(organ_name)
        mean_gene_organ=np.mean(np.array(organ_counts), axis=1)
        mean_genes.append(mean_gene_organ) 
        plt.plot(np.sort(mean_gene_organ),'o',markersize=1, label='{}'.format(organ_name))
    
    plt.legend(loc='upper left', fancybox=True, prop=fontp)
    plt.xlabel('gene')
    plt.ylabel('Mean expression')
    plt.yscale('log') 
    plt.title('Organs Rank plot')  
    plt.savefig(str(path_save_rank/'Organ_comparison.png'), dpi=300, format='png')
    plt.close()
        
    mean_genes=np.array(mean_genes).reshape(genes_num, len(organ_names))
    plt.figure(figsize=(8,5))
    plt.plot(np.sort(np.mean(mean_genes, axis=1)),'o',markersize=1, label='Global')
    plt.legend()
    plt.xlabel('gene')
    plt.ylabel('Mean expression')
    plt.yscale('log')
    plt.title('Global Rank plot')  
    plt.savefig(str(path_save_rank/'Global.png'))
    plt.close()

def main_save_density_data(path_data, path_save_density, path_save_mean_genes,path_save_pdf_genes, global_=False):
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    plt.figure(figsize=(10,7))
    for organ in listdir_[:]:
        path_data_=path_data/organ       
        organ_data=np.load(str(path_data_), allow_pickle=True)      
        organ_counts=organ_data[1:][:,1:].astype('float32')
        print(organ,'loaded')
        n_cells=organ_counts.shape[1]
        n_genes=organ_counts.shape[0]
        organ_name=organ.split(sep='.')[0]
        organ_name=organ_name.split(sep='-')[0]
        mean_gene_organ=np.mean(np.array(organ_counts), axis=1)/np.repeat(n_cells,n_genes)
        print('mean done')
        mean_gene_organ=np.sort(mean_gene_organ.astype('float32'))[::-1]
        #mean_gene_organ_not0=np.where(mean_gene_organ>0)
        pdf =norm.pdf(mean_gene_organ)
        print('pdf done')
        np.save(path_save_mean_genes/'mean_genes_expr_{}.npy'.format(organ_name), mean_gene_organ)
        np.save(path_save_pdf_genes/'pdf_mean_genes_expr_{}.npy'.format(organ_name),pdf)
#        plt.plot(mean_gene_organ, pdf, 'o', markersize=2, label=organ_name)
#        print(n_cells, n_genes , 'cells, genes', organ_name)
        del organ_data, organ_counts, mean_gene_organ, pdf
#    plt.legend()
#    plt.xlabel('mean gene expression value')
#    plt.ylabel('Density')
#    plt.title('Mean gene expression density plot')
#    plt.savefig(str(path_save_density/'Density_plot_not0.png'))
#    plt.close()

def main_density_plot(path_data,path_save_density, path_save_mean_genes,path_save_pdf_genes):
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    plt.figure(figsize=(10,7))
    for organ in listdir_[:][:1]:
        path_data_=path_data/organ   
        organ_name=organ.split(sep='.')[0]
        organ_name=organ_name.split(sep='-')[0]
        
        organ_data=np.load(str(path_data_), allow_pickle=True)      
        organ_counts=organ_data[1:][:,1:].astype('float32')
        cells_gene_expre=organ_counts[:,1]
        mean_gene_organ_not0=cells_gene_expre[np.where(cells_gene_expre>0)]
        pdf_not0 =norm.pdf(mean_gene_organ_not0)
        pdf_not0=pdf_not0[np.where]
        #path_save_mean_genes_=path_save_mean_genes/'mean_genes_expr_{}.npy'.format(organ_name)
        #path_save_pdf_genes_= path_save_pdf_genes/'pdf_mean_genes_expr_{}.npy'.format(organ_name)
#        mean_gene_organ=np.load(str(path_save_mean_genes_)).astype('float32')
#        #pdf_gene_organ=np.load(str(path_save_pdf_genes_)).astype('float32')
#        mean_gene_organ_not0=mean_gene_organ[np.where(mean_gene_organ>0)]
#        pdf_not0 =norm.pdf(mean_gene_organ_not0)
#        print('pdf done')
#        np.save(path_save_pdf_genes/'pdf_not0_mean_genes_expr_{}.npy'.format(organ_name),pdf_not0)
        
 #       pdf_not0 =norm.pdf(mean_gene_organ_not0)
        plt.plot(mean_gene_organ_not0, pdf_not0, 'o', markersize=2, label=organ_name)
#        print(n_cells, n_genes , 'cells, genes', organ_name)
        #del mean_gene_organ, mean_gene_organ_not0, pdf_not0
    plt.legend()
#    plt.xlabel('mean gene expression value')
#    plt.ylabel('Density')
#    plt.title('Mean gene expression density plot')
#    plt.savefig(str(path_save_density/'Density_plot_not0.png'))
#    plt.close()
    plt.xlabel('gene expression value')
    plt.ylabel('Density')
    plt.title('Gene expression density plot')
    plt.show()



def main_gene_expression_hist_plot(path_data, path_save_density, global_=True):
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    mean_genes=[]
    organ_names=[]
    plt.figure(figsize=(10,7))
    for organ in listdir_[:]:
        path_data_=path_data/organ       
        organ_data=np.load(str(path_data_), allow_pickle=True)      
        organ_counts=organ_data[1:][:,1:].astype('float32')
        genes_num=organ_counts.shape[0]
        organ_name=organ.split(sep='.')[0]
        organ_name=organ_name.split(sep='-')[0]
        genes_num=organ_counts.shape[0]
        organ_names.append(organ_name)
        mean_gene_organ=np.mean(np.array(organ_counts), axis=1)
        mean_genes.append(mean_gene_organ)
        if global_==False:
            plt.figure(figsize=(10,7))
            hist=np.histogram(mean_gene_organ, bins=len(np.unique(mean_gene_organ)))
            prob=hist[0]/np.repeat(np.sum(hist[0]), hist[0].shape[0])
            plt.plot(hist[1][15:][:-1], prob[15:])
            plt.xlabel('mean gene expression')
            plt.ylabel('Probability')
            plt.title('Gene Expression Density {}'.format(organ_name))
            plt.savefig(str(path_save_density/'{}.png'.format(organ_name)))
            plt.close()
    if global_==True:
        mean_genes=np.array(mean_genes).reshape(genes_num,len(listdir_))
        hist=np.histogram(mean_gene_organ, bins=len(np.unique(mean_gene_organ)))
        prob=hist[0]/np.repeat(np.sum(hist[0]), hist[0].shape[0])
        plt.plot(hist[1][15:][:-1], prob[15:])
        plt.xlabel('mean gene expression')
        plt.ylabel('Counts/N cells')
        plt.title('Global Gene Expression Histogram')
        plt.savefig(str(path_save_density/'Global.png'))
        plt.close()

def main_U_plots(path_data, path_save_U, global_fixed_thresh=None):
#    '''
#    This function produces a plot for each organ with the probability of occurrences of a gene in k cells,
#    or rather the number of times that a gene is found in k cells divided for the number of cells, 
#    where a gene is considered 'turned-on' according with the choosen parameter value.
#    '''
   listdir_ = os.listdir(path_data)
   listdir_=listdir_[:-1]
   organ_names=[]
   plt.figure(figsize=(12,9))
   for organ in listdir_[:]:
       plt.figure(figsize=(12,9))
       path_data_=path_data/organ 
       organ_data=np.load(str(path_data_), allow_pickle=True)
       organ_counts=organ_data[1:][:,1:]
       organ_name=organ.split(sep='.')[0]
       organ_name=organ_name.split(sep='-')[0]
       organ_names.append(organ_name)
       if global_fixed_thresh==None:
           U_plots(organ_counts,organ_name, path_save_U,0)
           U_plots(organ_counts,organ_name, path_save_U,20)
           U_plots(organ_counts,organ_name, path_save_U,300)
       #plt.yscale('log')
   if global_fixed_thresh==None:
       plt.legend(title='Parameter value', fancybox=True)
       plt.title('{} Gene occurrence in k cells for different thresholds'.format(organ_name))
       plt.savefig(str(path_save_U/'U_{}.png'.format(organ_name)))
       plt.close()
   else:
       U_plots(organ_counts,organ_name, path_save_U,global_fixed_thresh)
       plt.legend(title='Parameter value', fancybox=True)
       plt.title('Gene occurrence in k cells with threshold {}'.format(global_fixed_thresh))
       plt.savefig(str(path_save_U/'U_Thresh_{}.png'.format(global_fixed_thresh)))
 
def cumulative_gene_expression(path_data, path_save_rank, global_=True):  
    listdir_ = os.listdir(path_data)
    listdir_=listdir_[:-1]
    plt.figure(figsize=(12,9))
    mean_genes=[]
    for organ in listdir_[:]:
        organ_name=organ.split(sep='.')[0]
        organ_name=organ_name.split(sep='-')[0]

        path_data_=path_data/organ 
        organ_data=np.load(str(path_data_), allow_pickle=True)
        n_genes=23433
        organ_counts=organ_data[1:][:,1:]
        organ_counts=organ_counts/organ_counts.shape[1]
        mean_genes_organ=np.mean(organ_counts, axis=1)
        mean_genes.append(mean_genes_organ)
        cumsum= np.cumsum(np.sort(mean_genes_organ)[::-1])
        max_cumsum=np.max(cumsum)
        v_eff=0.72*max_cumsum
        idx= np.argwhere(v_eff-4<cumsum<v_eff+4)[0]
        print(idx)
        plt.plot(np.cumsum(np.sort(mean_genes_organ)[::-1]),label='{} {}'.format(organ_name, v_eff))
        plt.plot(v_eff, idx,'o')
    plt.legend()
    plt.title('Cumulative mean gene expression function')
    plt.xlabel('genes')
    plt.ylabel('mean gene expression/N cells')
    plt.savefig(str(path_save_rank/'Cumulative_gene_expression.png'))
    plt.close()
    
    if global_==True:
        plt.figure(figsize=(12,9))
        mean_genes=np.mean(np.array(mean_genes).reshape(n_genes,len(listdir_)),axis=1)
        plt.plot(np.cumsum(np.sort(mean_genes)[::-1]),label='Global')
        plt.legend()
        plt.title('Global cumulative mean gene expression function')
        plt.xlabel('genes')
        plt.ylabel('mean gene expression/N cells')
        plt.savefig(str(path_save_rank/'Cumulative_gene_expression.png'))
        plt.close()
        
def main():
    path_data=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/raw_matrices_FACS')
    path_save=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/analysis')
    path_mean_data=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/nean_var/CELLS/mean')
    path_save_mean_genes=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/mean_var/GENES/mean')
    path_save_pdf_genes=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/data/mean_var/GENES/pdf')
    path_save_density=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/analysis/Gene_Expression_Density')
    path_save_rank=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/analysis/Rank')
    path_save_U=Path('/home/slazzardi/Desktop/SPECIAL_ISSUE/mypyhon/analysis/U')
    #plot_cell_for_organ(path_data, path_save)
    #main_U_plots(path_data, path_save_U)
    #main_gene_expression_hist_plot(path_data, path_save_density)
    #cumulative_gene_expression(path_data, path_save, global_=True)
    #main_save_density_data(path_data, path_save_density, path_save_mean_genes, path_save_pdf_genes, global_=False)
    main_density_plot(path_data,path_save_density, path_save_mean_genes,path_save_pdf_genes)
    
if __name__=='__main__':
    main()