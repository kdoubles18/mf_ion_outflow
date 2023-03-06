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

global file_location
file_path = arg_parse_file_location()
file_location = f'{file_path.path}'+'/'+f'{file_path.file}'

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
    grid_plot = bats_file_location.add_grid_plot(i)
    
    
    return grid_plot

def save_grid_plot(i, outfile):
    
    fig = grid_plot(i)
    test = print('Writing file:' + file_path.path +  outfile)
    plt.savefig(outfile)
    
    return test, fig

bats_file = make_BatsFile(file_location)

def contour_plot(i,direction_1, direction_2, variable):
    make_BatsFile(file_location).switch_frame(i)
    fig = plt.figure(figsize=[4,4])
    fig, ax, cont, cbar = bats_file.add_contour(direction_1, direction_2, \
                                                variable, target=fig,\
                                                dolog = True, xlim=[-10,10], \
                                                ylim=[-10,10],add_cbar=True)
    return

def save_contour_plots(i,direction_1,direction_2,variable,outfile_2):
    make_BatsFile(file_location).switch_frame(i)
    fig = plt.figure(figsize=[4,4])
    fig, ax, cont, cbar = bats_file.add_contour(direction_1, direction_2, \
                                                variable, target=fig,\
                                                dolog = True, xlim=[-10,10], \
                                                ylim=[-10,10],add_cbar=True)
    test = print('Writing file:' + outfile_2)
    plt.savefig(outfile_2)

    return test

def plot_all_contour_sf(direction_1,direction_2,variable,outfile_2):
    '''
    This function allows the user to plot contour plots of specified variables
    and save these plots to a folder within the path of the file. This plots 
    all iterations of the specified file, which could be many depending on 
    interval of plotting frequency in run. You have been sufficiently warned.
    
        Args
        ---------
        direction_1 (string): selection of x, y, or z direction.
            Please write with ' ' around the direction.
        direction_2 (string): selection of x, y, or z direction.
            Please write with ' ' around the direction.
        variable (string): selection of the variable to plot.
            Please write with ' ' around the variable to contour plot.
    
        Return
        ---------
        Checks to see if there is a folder that already exists with the variable
        name specified. If the folder does not exist, one will be created with 
        the directions selected and the variable name where contour plots of 
        the variable selected will be created. If it does exist, an error 
        will be thrown and the plots will not save. If one wishes to plot a 
        variable again, they will need to remove the folder from the working
        directory.
    '''
    outs_file = make_BatsFile(file_location)
    for i in range(outs_file.attrs['nframe']):
        outs_file.switch_frame(i)
        fig = plt.figure(figsize=[10,10])
        
        fig.suptitle(f"%ith iteration - {variable}" %i)
        
        outs_file.add_contour(direction_1,direction_2,variable,target=fig,\
                              dolog = True, xlim=[-10,10], ylim=[-10,10],add_cbar=True)
                     
        try:
            os.makedirs(f'{direction_1}'+f'{direction_2}'+'_'+f'{variable}')
        except FileExistsError:
            # directory already exists
            pass

        test = print('Writing file:' + outfile_2)
        plt.savefig(file_path.path+'/'+f'{direction_1}'+f'{direction_2}'+'_'+f'{variable}'+'/'+f'{outfile_2}'+'_'+f'{i}')
        
    return

        
    
    



