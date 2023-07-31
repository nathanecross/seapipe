#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 12:07:36 2023

@author: nathancross
"""
from .spindle import whales
from .utils.audit import check_dataset, make_bids
import logging
import sys
import copy
from numpy import empty, zeros
from operator import itemgetter
from os import listdir, path, mkdir, walk
from pandas import DataFrame



class pipeline:
        
    """Contains specific information and allows the application of methods of 
    analysis, associated with a dataset. 

    Parameters
    ----------
    indir : str 
        name of the root level directory containing the BIDS organised data

    Attributes
    ----------
    rootpath : str
        name of the root level directory
    datapath : str
        name of the directory containing the raw data (recordings and annotations)
    outpath : str
        name of the directory containing the output (analysed) data

    """
        
    def __init__(self, indir, outfile=False):
        
        self.rootpath = indir
        self.datapath = indir + '/DATA/'
        self.outpath = indir + '/OUT/'
        self.audit = check_dataset(indir, outfile)
        
    def list_dataset(self, outfile=False): 
        
        """Prints out all the files inside the directory <in_dir> along with the
        directories 1 and 2 levels above containing the files. You can specify 
        an optional output filename that will contain the printout.
        """
        in_dir = self.datapath
        if outfile:
            logging.basicConfig(level=logging.DEBUG,
                                filemode='w',
                                format="%(message)s", 
                                handlers=[logging.StreamHandler(sys.stdout)])
            logger = logging.getLogger()
            file_log_handler = logging.FileHandler(f'{in_dir}/{outfile}')
            logger.addHandler(file_log_handler)
            stderr_log_handler = logging.StreamHandler()
            logger.addHandler(stderr_log_handler)
            formatter = logging.Formatter('%(message)s')
            file_log_handler.setFormatter(formatter)
            stderr_log_handler.setFormatter(formatter)
        
        else:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(message)s",
                                handlers=[logging.StreamHandler(sys.stdout)])
            logger = logging.getLogger()   
            
        for dirPath, dirNames, fileNames in walk(in_dir):
            try:
                fileNames.remove('.DS_Store')
            except(ValueError):
                pass
            
            if fileNames:
                dir1 = dirPath.split('/')[-2]
                dir2 = dirPath.split('/')[-1]
                logger.info(f"Directory: {dir1}/{dir2}")
                logger.info(f"Files = {fileNames}")
                logger.info('-' * 10)
                
    def make_bids(self, origin = 'SCN'):
        make_bids(self.datapath, origin=origin)
        
    def whale_it(self, xml_dir=False, out_dir=False, part='all', visit='all',
                 method='Lacourse2018', chan=['Cz'], ref_chan=None, rater=None, 
                 cat=(0,0,1,0), stage=['NREM2'], grp_name='eeg', cycle_idx=None, 
                 frequency=(11,16), adap_bands=False, duration=(0.5, 3)):
        
        # Set input/output directories
        in_dir = self.datapath
        if not xml_dir:
            xml_dir = self.datapath   
        if not out_dir:
            out_dir = self.outpath 
            
        whales.whale_it(in_dir, xml_dir, out_dir, method, chan, ref_chan, rater, 
                        cat, stage, grp_name, cycle_idx, frequency, adap_bands, 
                        duration, part, visit)

        