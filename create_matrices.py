#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 16:14:22 2019

@author: silvia
"""
import numpy as np
import pandas as pd
from pathlib2 import Path
import os

def read_data(path_data):
    data = pd.read_csv(path_data, sep=',', header=(0))
    return data

def save_matrices(path_data, path_save):
    listdir_data=os.listdir(path_data)
    print(len(listdir_data))
    for organ in listdir_data:
        if organ.split(sep='.')[1]=='csv':
            matrix_name=str(organ.split(sep='.')[0])
            path_data_=path_data/organ
            path_matrix = Path(path_save)/'{}.npy'.format(matrix_name)
            if not path_matrix.exists():
                data = pd.read_csv(path_data_, sep=',', header=(0))
                print(data.shape)
                np.save(str(path_save/'{}'.format(matrix_name)), data)
            else:
                print('{} already exists.'.format(matrix_name))
        
def main():
#    import argparse
#    parser = argparse.ArgumentParser()
#    parser.add_argument('--path_data')
#    parser.add_argument('--path_save')
#    args = parser.parse_args()

    # VARIABLES     
    #path_data = args.path_data
    #path_save= args.path_save
    path_data=Path('/home/silvia/Desktop/BORSA/SPECIAL_ISSUE/mypyhon/00_facs_raw_data/FACS/FACS')
    path_save=Path('/home/silvia/Desktop/BORSA/SPECIAL_ISSUE/mypyhon/data/row_matrices')
    save_matrices(path_data, path_save)
    
if __name__=='__main__':
    main()