#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 14:26:53 2025

@author: ncro8394
"""

from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.trans import fetch

file = '/Users/ncro8394/Documents/projects/seapipe/DATA/sub-IN001/ses-V1/eeg/sub-IN001_ses-V1_eeg.edf'
d = Dataset(file)



annot = Annotations('/Users/ncro8394/Documents/projects/seapipe/derivatives/slowwave/sub-IN001/ses-V1/sub-IN001_ses-V1_slowosc.xml')


segs = fetch(d, annot, cat=(0,0,0,0), evt_type = ['Staresina2015'],
             chan_full=)