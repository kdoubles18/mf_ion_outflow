#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 18:00:26 2023

@author: kdoubles
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 19:50:42 2022

@author: kdoubles
"""

__author__ = 'Kailtin Doublestein'
__email__ = 'kdoubles@umich.edu'

import numpy as np
import spacepy
from spacepy import pybats
from spacepy.pybats import bats, rim, IdlFile, calc_wrapper
import spacepy.datamodel as dm
from spacepy.plot import set_target, applySmartTimeTicks
from matplotlib.gridspec import GridSpec
from matplotlib import gridspec
import matplotlib.pyplot as plt
import argparse as ap
import os, glob
from operator import itemgetter, attrgetter
from datetime import datetime

def arg_parse_file_location():
    """
    A parser that allows the user to type in the path of the file that they 
    want to process. User must input path and file name if file is not in same
    working directory as this script.
        
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
    parser.add_argument('-runname', help = 'This is the name of the run. \
                        Default is "None".', type = str, \
                            required = True)
    args = parser.parse_args()

    return args

is_multi = True
if input("Is this a multifluid run? (y/n) ") != "y":
    is_multi = False
    print(is_multi)

global file_location
directory = arg_parse_file_location()
run_name = directory.runname

'''    
#Read in the y=0 .outs
for filename in glob.glob((directory.path + '/' + 'y*.outs')):
    y_2d = bats.Bats2d(filename)
    y_idl = bats.IdlFile(filename)

#Make new folder for plots to be stored in
isExist = os.path.exists(directory.path + '/plots/')
if not isExist:

   # Create a new directory because it does not exist
   os.mkdir(directory.path + '/plots/')
         
#Density for the three types of fluids for this run. All three are in a single
#plot together.
for i in range(0,y_2d.attrs['nframe'],y_2d.attrs['nframe']):
    y_2d.switch_frame(i)
    grid = gridspec.GridSpec(3,1)
    fig = plt.figure(figsize=(4,12))
    fig.subplots_adjust(hspace = 0.3)
    
    t_now = y_2d.attrs['runtimes'][i]
    t_now_hr = t_now/3600
    
    ax1, ax2, ax3 = fig.add_subplot(grid[0,0]),\
                    fig.add_subplot(grid[1,0]), \
                    fig.add_subplot(grid[2,0])
    fig, ax, cont, cbar = y_2d.add_contour('x','z','rho',target = ax1,\
                                xlim=[-10,10], ylim=[-10,10],add_cbar=True)
    fig, ax, cont, cbar = y_2d.add_contour('x','z','hprho',target = ax2,\
                                xlim=[-10,10], ylim=[-10,10],add_cbar=True)
    fig, ax, cont, cbar = y_2d.add_contour('x','z','oprho',target = ax3,\
                                xlim=[-10,10], ylim=[-10,10],add_cbar=True)

    y_2d.add_b_magsphere(target=ax1,colors='Black',DoLast=False)
    y_2d.add_b_magsphere(target=ax2,colors='Black',DoLast=False)
    y_2d.add_b_magsphere(target=ax3,colors='Black',DoLast=False)

    plt.suptitle('Y=0, Density of Fluid Types - \
                 {}'.format(t_now_hr), fontsize=10)
    
    plt.savefig('{}/plots/y0_dens_{}'.format(directory.path,i), dpi=300)
    plt.close()
    print('plotting frame: {}'.format(i))

#Get spherical coordinate values
y_2d['r']  = np.sqrt(y_2d['x']**2+y_2d['z']**2)
y_2d['phi'] = 0
y_2d['theta'] = np.cos(y_2d['z']/np.sqrt(y_2d['x']**2+y_2d['z']**2))**(-1)

#Calculate radial velocity.
ur = y_2d['ux']*np.sin(y_2d['theta'])*np.cos(y_2d['phi']) + \
   y_2d['uy']*np.sin(y_2d['theta'])*np.sin(y_2d['phi']) + \
   y_2d['uz']*np.cos(y_2d['theta'])
y_2d['ur'] = dm.dmarray(ur, {'units':y_2d['ux'].attrs['units']})
'''
'''
#Read in the log file for dst and cpcpn
for filename in glob.glob((directory.path + '/' + 'log*.log')):
    log_file = bats.LogFile(filename)

#Plot dst vs time.
fig = plt.figure(figsize=(4,4))
plt.plot(log_file['time'],log_file['dst'])
plt.ylabel('Dst')
plt.xlabel('Time, t, [s]')
plt.title('Dst & Time')
fig.autofmt_xdate(rotation=45)
plt.savefig('{}/plots/dst_time'.format(directory.path), dpi=300)
plt.close()
#plot cpcpn vs time.
fig = plt.figure(figsize=(4,4))
plt.plot(log_file['time'],log_file['cpcpn'])
plt.ylabel('Cpcp')
plt.xlabel('Time, t, [s]')
plt.title('Cpcp & Time')
fig.autofmt_xdate(rotation=45)
plt.savefig('{}/plots/cpcp_time'.format(directory.path), dpi=300)
plt.close()
'''
#If no folder of plots exists in /IE/ subfolder, then this will make one.
isExist = os.path.exists(directory.path + 'IE/plots/')
if not isExist:

   # Create a new directory because it does not exist
   os.mkdir(directory.path + '/IE/plots/')

#dictionary for iono files and time associated with file.
iono_file_dict = {'iono_file':[],'time':[]}

#find all the files in the /IE/ folder and read them in as a Pybats object,
#then appead them to the dictionary.
for filename in glob.glob((directory.path + 'IE/' + '*.idl')):
    #Print if you want to see that it works 
    print('read file: {}'.format(filename))
    iono = rim.Iono(filename)
    iono_file_dict['iono_file'].append(iono)
    iono_file_dict['time'].append(iono.attrs['time'])

#plotting all the iono files in the dictionary, not sorted out by time.
for i in enumerate(iono_file_dict['iono_file']):
    fig = plt.figure(figsize=(6,7))
    i[1].add_cont('n_jr', target=fig, add_cbar=True)

    # Check that the times match.
    #   Compares the time stored in the dictionary with the attribute.
    #   Can be removed; just a demonstration.
    print('\t Compare: {}  {}'.format(i[1].attrs['time'], iono_file_dict['time'][i[0]]))
    
    plt.title('N Jr - {}'.format(i[1].attrs['time']))
    print('plotting frame {:03d} at time {}'.format(i[0], i[1].attrs['time']))
    plt.savefig('{}/plots/ie_njr_{}.png'.format(directory.path, i[0]), dpi=300)

    #To plot with time in the filename:
    #plt.savefig('{}ie_njr_{}.png'.format(directory, i[1].attrs['time'].strftime("%Y%m%d-%H%M%S")), dpi=300)
    
    plt.close()
