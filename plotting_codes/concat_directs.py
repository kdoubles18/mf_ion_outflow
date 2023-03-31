#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 20:14:19 2023

@author: kdoubles
"""

'''
This script is used to concatinate files from multiple run directories. The
intended use is to combine directories caused by restarting a run. There are
two instances: 1) when the RESTART does not occur on an interval of plotting 
and 2) when the RESTART does occur on the interval of plotting. The methods 
for concatinating are different so pick the correct one.
'''

def concat_not_interval():
    '''
    This definition will concatinate files where the plotting time scale is
    NOT on the same timing interval as the RESTART file.
    
    Args
    -------
    parent_direct_not (str):
        name of the directory where all the output directories are located.
        Example: /test_directory/, where inside /test_directory/ exists /run_1
        /run2 .... /run_n, with the outputs that need to be concatinated.

    Returns
    -------
    output_direct_not []:
        The directories within parent_direct are concatinated into a single
        directory and saved within the parent_direct.
    '''
    
    
    
    return


def concat_on_interval():
    '''
    This definition will concatinate files where the plotting time scale IS
    on the same timing interval as the RESTART file.
    
    Args
    -------
    parent_direct_on (str):
        name of the directory where all the output directories are located.
        Example: /test_directory/, where inside /test_directory/ exists /run_1
        /run2 .... /run_n, with the outputs that need to be concatinated.

    Returns
    -------
    output_direct_on []:
        The directories within parent_direct are concatinated into a single
        directory and saved within the parent_direct.
    '''
    
    
    return