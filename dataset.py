#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 12:07:36 2023

@author: nathancross
"""

from stats import sleepstats
from spindle import whales
from utils.audit import check_dataset, extract_channels, make_bids
from utils.logs import create_logger, create_logger_outfile
from os import mkdir, path, remove, walk
import logging
import sys


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
        if not path.exists(self.outpath):
            mkdir(self.outpath)
        self.outfile = outfile
        self.audit_init = check_dataset(self.datapath, self.outfile)
        
       
    def audit(self, outfile=False):
        
        ''' Audits the dataset for BIDS compatibility.
            Includes option to save the audit to an output file.
        '''

        if not outfile and not self.outfile:
            logger = create_logger()
            self.audit_update = check_dataset(self.datapath, False, logger)
        else:
            if not outfile:
                outfile = self.outfile
            out_dir = f'{self.outpath}/audit'
            if not path.exists(out_dir):
                mkdir(out_dir)
            out = f'{out_dir}/{outfile}'
            if path.exists(out):
                remove(out)
            logger = create_logger()
            self.audit_update = check_dataset(self.datapath, out, logger)
            
        logger.info('')
        logger.info(self.audit_update)
        
        
    def list_dataset(self, outfile=False): 
        
        """Prints out all the files inside the directory <in_dir> along with the
        directories 1 and 2 levels above containing the files. You can specify 
        an optional output filename that will contain the printout.
        """

        if not outfile and not self.outfile:
            logger = create_logger()  
        else:
            if not outfile:
                outfile = self.outfile
            out_dir = f'{self.outpath}/audit'
            if not path.exists(out_dir):
                mkdir(out_dir)
            out = f'{out_dir}/{outfile}'
            if path.exists(out):
                remove(out)
            logger = create_logger_outfile(out)

        logger.propagate = False
        
        for dirPath, dirNames, fileNames in walk(self.datapath):
            try:
                fileNames.remove('.DS_Store')
            except(ValueError):
                pass
            
            if fileNames:
                dir1 = dirPath.split('/')[-3]
                dir2 = dirPath.split('/')[-2]
                dir3 = dirPath.split('/')[-1]
                logger.info(f"Directory: {dir1}/{dir2}/{dir3}")
                logger.info(f"Files •, {fileNames}")
                logger.info('-' * 10)
                
    def make_bids(self, origin = 'SCN'):
        make_bids(self.datapath, origin=origin)
        
    def extract_channels(self, exclude = None):
        extract_channels(self.datapath, exclude=exclude)
        
    def sleepstats():
        return
        
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

        