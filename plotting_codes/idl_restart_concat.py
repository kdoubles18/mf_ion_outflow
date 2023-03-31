#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 10:56:11 2023

@author: kdoubles
"""

'''
This script is designed to read in IDL output files from multiple restarts 
from the SWMF, remove the headers, and remove any overlapping data from the
restarted files. Make sure files 
'''
from spacepy import pybats
from spacepy.pybats import bats
import argparse as ap
import os, glob
import re

def arg_parse_file_location():
    """
    A parser that allows the user to type in the path of the file that they 
    want to process. User must input path name to the level where the output
    files are.
        
    Args
    ------
    path_name (str): 
        Name of path to the directory of the data.
            
    Return
    -------
    Arguments for the file and the path. 
    
    """
    
    parser = ap.ArgumentParser(description = 'Read in path name')
    parser.add_argument('-path', help = 'The path for the file. Default \
                        directory of this script.', type = str, \
                            default = 'None')
    args = parser.parse_args()

    return args

global file_location
directory = arg_parse_file_location()

mag_grid_list = {'filename': [],'idl_file':[],'mag_grid':[]}
## Mag Grid File
for filename in glob.glob((directory.path + '/' +'mag_grid*.outs')):
    mag_grid_list['filename'].append(filename)
    
for filename in mag_grid_list['filename']:
    mag_grid_list['idl_file'].append(pybats.IdlFile(filename))
    
for filename in mag_grid_list['filename']:
    mag_grid_list['mag_grid'].append(bats.MagGridFile(filename))
    

for n in range(0,len(mag_grid_list['filename'])):

    if mag_grid_list['idl_file'][n].attrs['format'] == 'asc':
        pybats._read_idl_ascii(mag_grid_list['idl_file'][n], 
                        header=mag_grid_list['idl_file'][n]._header,
                        keep_case=mag_grid_list['idl_file'][n]._keep_case)
    elif mag_grid_list['idl_file'][n].attrs['format'] == 'bin':
        pybats._read_idl_bin(mag_grid_list['idl_file'][n], 
                      header=mag_grid_list['idl_file'][n]._header,
                      keep_case=mag_grid_list['idl_file'][n]._keep_case)
    else:
        raise ValueError('Unrecognized file format: {}'.format(
            mag_grid_list['idl_file'][n]._format))