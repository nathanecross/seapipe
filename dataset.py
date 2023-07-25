#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 12:07:36 2023

@author: nathancross
"""

import copy
from numpy import empty, zeros
from operator import itemgetter
from os import listdir, path, mkdir
from pandas import DataFrame


def check_dataset(in_dir, outfile=False):
    
    print('')
    print(f'Checking dataset in directory: {in_dir}')
    print('')
    
    part = [x for x in listdir(in_dir) if not '.' in x]
    
    # Check number of subdirectories per participant
    nsd = zeros((len(part)))
    for p, pt in enumerate(part):
        
        subs = [x for x in listdir(f'{in_dir}/{pt}') if not '.' in x]
        nsd[p] = len(subs)
        
        if len(subs) < 1:
            print('')
            print(f'WARNING: {pt} has zero visits')
            print('')
    
    if len(set(nsd)) > 1:  
        print('')
        print('WARNING: Not all participants have the same number of visits')
        print('')
    
    # List of participants and number of visits
    subdirs = DataFrame(zeros((len(part),1)),index=part,columns=['#subdirectories'])
    subdirs.iloc[:,0] = nsd
    if outfile:
        subdirs.to_csv(f'{in_dir}/dataset_audit.csv')
    
    return subdirs


class dataset:
        
    def __init__(self, indir):
        self.audit = check_dataset(indir)