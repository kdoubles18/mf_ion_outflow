#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:27:00 2023

@author: kdoubles
"""
## Might not need all these, idk
import os
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from spacepy import pybats
import spacepy.plot as splot
import spacepy.datamodel as dm
from spacepy.pybats import bats
import warnings
from matplotlib import MatplotlibDeprecationWarning
from matplotlib.gridspec import GridSpec
import argparse as ap
import os, glob
#import shell_plotter.py as Shell

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
    args = parser.parse_args()

    return args

global file_location
directory = arg_parse_file_location()

'''
Convert all files in directory location to usable formats with SpacePy.
'''

## Mag Grid File
for filename in glob.glob((directory.path + '/' +'mag_grid*.outs')):
    mag_grid = bats.MagGridFile(filename)
    
## Geo Index File
for filename in glob.glob((directory.path + '/' +'geo*.log')):
    geo_index = bats.GeoIndexFile(filename)

## Log File
for filename in glob.glob((directory.path + '/' +'log*.log')):
    log_file = bats.BatsLog(filename)
 
## Shell Slice
#for filename in glob.glob((directory.path + '/' +'shl_*.outs')):
#    shell_file = Shell.ShellSlice(filename)

#Total number of frames.
nFrame = mag_grid.attrs['nframe']
print('nFrame: {}'.format(nFrame))

## Dictionary of needed values.
mag_grid_dict = {'times':[],'lat_avg':[],'dbdtn':[],'dbdte':[],'dbdth':[]}

dBdth = np.empty((mag_grid['dBn'].shape[0],mag_grid['dBn'].shape[1],nFrame))

## Iterate through frames in .outs file
##  -> could modify to iterate through .out files too
for iFrame in range(nFrame):

    mag_grid.switch_frame(iFrame)
    
    ## get simulation time
    t_sim = mag_grid.attrs['runtimes'][iFrame]

    if iFrame % 100 == 0: # avoid screen barf
        print('iFrame: {}'.format(iFrame))
        print('tSimulation: {}'.format(t_sim))
        
    ## Calc dBh
    mag_grid.calc_h()

    ## Average over latitude
    mag_grid['lat_avg'] = np.mean(mag_grid['dBh'], axis=0, dtype=np.float64)
    
    ## Append to lists for plotting
    mag_grid_dict['times'].append(t_sim)
    mag_grid_dict['lat_avg'].append(mag_grid['lat_avg'])
  
delta_t_0 = mag_grid.attrs['runtimes'][1]+mag_grid.attrs['runtimes'][0]
dBdth[:,:,0] = (mag_grid['dBh'][2]+4*mag_grid['dBh'][1]-
                3*mag_grid['dBh'][0])/(delta_t_0)    
mag_grid_dict['dbdth'].append(dBdth[:,:,0])

#Central difference approximation
for iFrame in range(1,nFrame-1):
    t_1 = mag_grid.attrs['runtimes'][iFrame-1]
    t_2 = mag_grid.attrs['runtimes'][iFrame+1]
    delta_t = t_2 - t_1
    
    if iFrame % 100 == 0: # avoid screen barf
        print('iFrame: {}'.format(iFrame))
        print('dbdth size: {}'.format(len(mag_grid_dict['dbdth'])))

    dBdth[:,:,iFrame] = (mag_grid['dBh'][:,:])/(2*delta_t)
    mag_grid_dict['dbdth'].append(dBdth[:,:,iFrame])
        
delta_t_0 = (mag_grid.attrs['runtimes'][nFrame-2]+
             mag_grid.attrs['runtimes'][nFrame-1])
dBdth[:,:,-1] = (3*mag_grid['dBh'][-1]-4*mag_grid['dBh'][-2]+
                 mag_grid['dBh'][-3])/(delta_t_0)    
mag_grid_dict['dbdth'].append(dBdth[:,:,-1])

for iFrame in range(nFrame):
    mag_grid['lat_avg_dbdth'] = np.mean(dBdth, axis=0, dtype=np.float64)
    mag_grid['lon_avg_dbdth'] = np.mean(dBdth, axis=1, dtype=np.float64)

'''
Begin plotting all the plots I need for analysis.
'''
## Plot Northern Hemisphere
fig, ax = plt.subplots(figsize=[14,8])
for iLat in np.arange(85,171,10): ## start from Equator, do every 10 deg
	ax.plot(mag_grid_dict['times'],
        	np.asarray(mag_grid_dict['lat_avg'])[:,iLat],
                label = 'Lat = {}'.format(iLat-85))
ax.legend()
ax.set_xlabel('Simulation Time [s]')
ax.set_ylabel('dBh [nT]')
ax.set_title('Northern Hemisphere - 10 eV IBC')
plt.savefig('{}/lat_avg_northern_10eV.png'.format(directory.path),dpi=300)
#plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(14,8))
for iLat in np.arange(0,86,10): ## end at Equator, do every 10 deg
	ax.plot(mag_grid_dict['times'],
        	np.asarray(mag_grid_dict['lat_avg'])[:,iLat],
                label = 'Lat = {}'.format(iLat-85))
ax.legend()
ax.set_xlabel('Simulation Time [s]')
ax.set_ylabel('dBh [nT]')
ax.set_title('Southern Hemisphere - 10 eV IBC')
plt.savefig('{}/lat_avg_southern_10eV.png'.format(directory.path),dpi=300)
#plt.show()
plt.close()

#Some subplot stuff
fig = plt.figure(figsize=(14,20), constrained_layout=True)
gs0 = fig.add_gridspec(9,1, figure=fig)
plt.suptitle('Run 02 - dBh, 10 eV IBC',fontsize=20)
for n in np.arange(0,9):  
    ax0 = fig.add_subplot(gs0[n]) 
    iLat = 10*n + 85 ## or something; I can't test my algebra
    ax0.plot(mag_grid_dict['times'],
                 np.asarray(mag_grid_dict['lat_avg'])[:,iLat],
                 label = 'Lat = {}'.format(iLat-85))
    ax0.legend()

ax0.set_xlabel('Simulation time [s]',fontsize=20)
ax0.set_ylabel('dBh [nT]', fontsize=20)
plt.savefig('{}/lat_avg_N_separated_10eV.png'.format(directory.path),dpi=300)    
#plt.show()
plt.close()

fig = plt.figure(figsize=(14,20), constrained_layout=True)
gs1 = fig.add_gridspec(9,1, figure=fig)
plt.suptitle('Run 02 - dBh, 10 eV IBC',fontsize=20)
for n in np.arange(0,9):  
    ax1 = fig.add_subplot(gs1[n]) 
    iLat = 10*n ## or something; I can't test my algebra
    ax1.plot(mag_grid_dict['times'],
                 np.asarray(mag_grid_dict['lat_avg'])[:,iLat],
                 label = 'Lat = {}'.format(iLat-85))
    ax1.legend()

ax1.set_xlabel('Simulation time [s]',fontsize=20)
ax1.set_ylabel('dBh [nT]', fontsize=20)
plt.savefig('{}/lat_avg_S_separated_10eV.png'.format(directory.path),dpi=300)    
#plt.show()
plt.close()
            
## Calc dBh
geo_index.fetch_obs_ae()
geo_index.fetch_obs_kp()

fig = plt.figure(figsize=[8,14],constrained_layout=True)
gs = GridSpec(3,1, figure=fig)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2,0])

ax1.plot(geo_index['time'],geo_index['Kp'])
ax1.set_ylabel('Kp')
ax1.tick_params(labelbottom=False)

ax2.plot(geo_index['time'],geo_index['AU'], label='AU',marker='o')
ax2.plot(geo_index['time'],geo_index['AL'], label='AL',marker='*')
ax2.set_ylabel('AE')
ax2.tick_params(labelbottom=False)
ax2.legend()

ax3.plot(log_file['time'],log_file['dst'])
ax3.set_xlabel('time')
ax3.set_ylabel('Dst, [nT]')
ax3.tick_params(rotation=30)

plt.suptitle('28 amu/cc Density - 10 eV Temperature IBC')

plt.savefig('{}/gen_'.format(directory.path),dpi=300) 
#plt.show()
plt.close()

## plot dBdth
plt.figure()
for iLat in np.arange(85,171,10): ## start from Equator, do every 10 deg
	plt.plot(mag_grid_dict['times'],
        	np.asarray(mag_grid['lat_avg_dbdth'][iLat,:]),
                label = 'Lat = {}'.format(iLat-85))
plt.xlim(1,35999)
plt.show()

## dB_h / dt 
#Adapted from Pulkkinen et al. 2013, where B_x and B_y are B_n and B_e.
'''
    ## No times were saved in the file you gave me.
    ## If you want actual times, do this:
    #epoch = dt.datetime(YYYY, MM, DD,
    #                    HH, MM, SS, MSC) ## starttime of the run
    #t_sim = epoch + dt.timedelta(seconds=t_sim))
    
    #dt = np.array([x.total_seconds() for x in np.diff(geo_index['time'])])

'''
