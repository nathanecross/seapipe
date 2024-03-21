#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:10:35 2024

@author: nathancross
"""

import logging
import sys



class CustomFormatter(logging.Formatter):
    
    cyan = "\x1b[36;20m"
    grey = "\x1b[38;20m"
    yellow = "\x1b[1;33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format1 = "%(message)s"
    format2 = "%(asctime)s - %(name)s - %(message)s "
    format3 = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    format4 = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.INFO: cyan + format1 + reset,
        logging.DEBUG: grey + format2 + reset,
        logging.WARNING: yellow + format3 + reset,
        logging.ERROR: bold_red + format4 + reset,
        logging.CRITICAL: bold_red + format3 + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt,"%Y-%m-%d %H:%M:%S")
        return formatter.format(record)

class CustomLogfileFormatter(logging.Formatter):
    
    format1 = "%(message)s"
    format2 = "%(asctime)s - %(name)s - %(message)s "
    format3 = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    format4 = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.INFO: format1,
        logging.DEBUG: format2,
        logging.WARNING: format3,
        logging.ERROR: format4,
        logging.CRITICAL: format3
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt,"%Y-%m-%d %H:%M:%S")
        return formatter.format(record)

def create_logger(name):
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    
    ch.setFormatter(CustomFormatter())
    
    logger.addHandler(ch)
    logger.propagate = False
    return logger

def create_logger_empty():
    
    # Set logging 
    logging.basicConfig(level=logging.DEBUG,
                        format="",
                        handlers=[logging.StreamHandler(sys.stdout)])
    logger = logging.getLogger()
    logger.handlers.clear()
    logger.propagate = False
    
    return logger

def create_logger_basic():
    
    # Set logging 
    logging.basicConfig(level=logging.DEBUG,
                        format="%(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    logger = logging.getLogger()
    logger.handlers.clear()
    logger.propagate = False
    
    return logger

def create_logger_outfile(logfile, name=None):

    # logging.basicConfig(level=logging.DEBUG,
    #                     format="%(message)s", 
    #                     handlers=[logging.StreamHandler(sys.stdout)])
    if name:
        logger = logging.getLogger(name)
    else:
        logger = logging.getLogger()
    
    logger.handlers.clear()
    
    # Logger to print to outfile
    file_log_handler = logging.FileHandler(f'{logfile}', mode='a')
    file_log_handler.setLevel(logging.DEBUG)
    file_log_handler.setFormatter(CustomLogfileFormatter())
    logger.addHandler(file_log_handler)
    
    # Logger to print to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(CustomFormatter())
    logger.addHandler(ch)
    
    
    #logger.propagate = False

    return logger