#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:10:35 2024

@author: nathancross
"""

import logging
import sys


def create_logger():
    # Set logging 
    logging.basicConfig(level=logging.DEBUG,
                        format="%(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    logger = logging.getLogger()
    #logger.handlers.clear()
    logger.propagate = False
    
    return logger

def create_logger_outfile(logfile):
    logging.basicConfig(level=logging.DEBUG,
                        filemode='w',
                        format="%(message)s", 
                        handlers=[logging.StreamHandler(sys.stdout)])
    logger = logging.getLogger()
    logger.handlers.pop()
    file_log_handler = logging.FileHandler(f'{logfile}')
    logger.addHandler(file_log_handler)
    stderr_log_handler = logging.StreamHandler()
    logger.addHandler(stderr_log_handler)
    formatter = logging.Formatter('%(message)s')
    file_log_handler.setFormatter(formatter)
    stderr_log_handler.setFormatter(formatter)

    return logger