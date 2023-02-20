#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:24:32 2023

@author: kdoubles

Designed to plot data from SWMF .outs outputs for y=0, z=0 and shl
files. 
"""

#imports
import numpy as np
import matplotlib.pyplot as plt
import spacepy.plot as splot
from matplotlib import gridspec
import spacepy.datamodel as dm
import argparse as ap
import os

from spacepy import pybats
from spacepy.pybats import bats

def arg_parse_file_location():
    """
    A parser that allows the user to type in the path of the file that they 
    want to process. User must input path and file name if file is not in same
    working directory as this script.
        
    Args
    ------
    file (str):
        Input the path and file name if file is not in same working
        directory as this script. Otherwise, input the file name.
            
    Return
    -------
    Arguments for the file and the path. 
    
    """
    
    parser = ap.ArgumentParser(description = 'Read in file name')
    parser.add_argument('-file', \
                      help = 'The name of the file that you want to \
                      evaluate', type = str,default='None')
    parser.add_argument('-path', help = 'The path for the file. Default \
                        directory of this script.', type = str, default = \
                        os.getcwd())
    
    args = parser.parse_args()

        

    return args
    return args

def make_BatsFile(file_location):
    """
    Load file from arg parse and translate to Spacepy bats File format.
    Helpful for plotting with the Bats2d options.
        Args
        -------
        file_name (str):
            This comes from the argumentx` parser input.
        Return
        --------
        bats2d file format.
    """
    
    bats_file = bats.Bats2d(file_location)
    
    return bats_file

def grid_plot(i):
    """
    Plot the grid layout of the run. 
    
    Args
    -------
    bats_file ():
        Comes from spacepy object 'spacepy.pybats.bats.Bats2d'
    
    Returns
    -------
    Grid plot produced from Spacepy grid layout plotting function.

    """
    
    bats_file_location = make_BatsFile(file_location)
    grid_plot = bats_file_location.add_grid_plot()
    
    return grid_plot

def save_grid_plot(i, outfile):
    
    fig = grid_plot(i)
    test = print('Writing file:' + outfile)
    plt.savefig(outfile)
    
    return test, fig

    


global file_location
file_path = arg_parse_file_location()
file_location = f'{file_path.path}'+'/'+f'{file_path.file}'

bats_file = make_BatsFile(file_location)

fig = plt.figure(figsize=[8,4])
fig, ax, cont, cbar = bats_file.add_contour('x', 'z', 'opp', target=fig,loc = 121, dolog = True, xlim=[-10,10], ylim=[-10,10],add_cbar=True)






