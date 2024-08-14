#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 10:43:09 2024

@author: ncro8394
"""

from multiprocessing import Pool


def parallel(FUNC, cores, args):

    pool = Pool(cores)
    results = pool.starmap(FUNC, args)
    
    return results