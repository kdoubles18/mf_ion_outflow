#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created Thurs 16 Mar 2023
@author: tkeebler
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
warnings.filterwarnings("ignore",category=MatplotlibDeprecationWarning)

## To make Dan happy or something, idk
splot.style('default')

## Change to match your own directory structure
directory = '//Users/kdoubles/Data/run_results/run_01_full'
filename = 'mag_grid_n00009569_00225137.outs'

## Read file, check number of frames.
mag_grid = bats.MagGridFile('{}/{}'.format(directory,filename))
nFrame = mag_grid.attrs['nframe']
print('nFrame: {}'.format(nFrame))

## Dictionary to hold stuff since I'm lazy
## and rationalizing the familiar is easy
mag_grid_dict = {'times':[],'lat_avg':[]}

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

    ## No times were saved in the file you gave me.
    ## If you want actual times, do this:
    #epoch = dt.datetime(YYYY, MM, DD,
    #                    HH, MM, SS, MSC) ## starttime of the run
    #t_sim = epoch + dt.timedelta(seconds=t_sim))

    ## Append to lists for plotting
    mag_grid_dict['times'].append(t_sim)
    mag_grid_dict['lat_avg'].append(mag_grid['lat_avg'])


## For testing:
#print(np.shape(mag_grid['dBh']))
#print(mag_grid.attrs)
#print(np.shape(mag_grid['lat_avg']))
#print(np.shape(mag_grid_dict['times']))
#print(np.shape(np.asarray(mag_grid_dict['lat_avg'])[:,0]))

## Plot Northern Hemisphere
fig, ax = plt.subplots(figsize=[14,8])
for iLat in np.arange(85,171,10): ## start from Equator, do every 10 deg
	ax.plot(mag_grid_dict['times'],
        	np.asarray(mag_grid_dict['lat_avg'])[:,iLat],
                label = 'Lat = {}'.format(iLat-85))
ax.legend()
ax.set_xlabel('Simulation Time [s]')
ax.set_ylabel('dBh [nT]')
ax.set_title('Northern Hemisphere - Nominal')
plt.savefig('{}/lat_avg_northern.png'.format(directory),dpi=300)
plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(14,8))
for iLat in np.arange(0,86,10): ## end at Equator, do every 10 deg
	ax.plot(mag_grid_dict['times'],
        	np.asarray(mag_grid_dict['lat_avg'])[:,iLat],
                label = 'Lat = {}'.format(iLat-85))
ax.legend()
ax.set_xlabel('Simulation Time [s]')
ax.set_ylabel('dBh [nT]')
ax.set_title('Southern Hemisphere - Nominal')
plt.savefig('{}/lat_avg_southern.png'.format(directory),dpi=300)
#plt.show()
plt.close()


#Some subplot cool stuff
fig = plt.figure(figsize=(14,20), constrained_layout=True)
gs0 = fig.add_gridspec(9,1, figure=fig)
plt.suptitle('Run 0q - dBh, Nominal',fontsize=20)
for n in np.arange(0,9):  
    ax0 = fig.add_subplot(gs0[n])
    iLat = 10*n + 85 ## or something; I can't test my algebra
    ax0.plot(mag_grid_dict['times'],
                 np.asarray(mag_grid_dict['lat_avg'])[:,iLat],
                 label = 'Lat = {}'.format(iLat-85))
    ax0.legend()

ax0.set_xlabel('Simulation time [s]',fontsize=20)
ax0.set_ylabel('dBh [nT]', fontsize=20)
plt.savefig('{}/lat_avg_N_separated.png'.format(directory),dpi=300)    
#plt.show()
plt.close()






