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
    parser.add_argument('-runname', help = 'This is the name of the run. \
                        Default is "None".', type = str, \
                            required = True)
    args = parser.parse_args()

    return args

if input("Is this multifluid (y/n)") != "y":
    is_multi = False

global file_location
directory = arg_parse_file_location()
run_name = directory.runname

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

## Dictionary to hold stuff since I'm lazy
## and rationalizing the familiar is easy
mag_grid_dict = {'times':[],
                 'dBh':[],'dBh_lat_avg':[]}
## Iterate through frames in .outs file
##  -> could modify to iterate through .out files too
for iFrame in range(nFrame):

    mag_grid.switch_frame(iFrame)
    
    ## get simulation time
    t_sim = mag_grid.attrs['runtimes'][iFrame]
<<<<<<< HEAD

=======
    
>>>>>>> 3ecf3d8 (updating code to allow input naming process)
    if iFrame % 100 == 0: # avoid screen barf
        print('iFrame: {}'.format(iFrame))
        print('tSimulation: {}'.format(t_sim))
        
    ## Calc dBh
    mag_grid.calc_h()
<<<<<<< HEAD

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

    mag_grid_dict['dBh'].append(mag_grid['dBh'])
    
    ## Compute latitude-averaged dBh
    mag_grid['dBh_lat_avg'] = np.mean(mag_grid['dBh'], axis=0, dtype=np.float64)
    
    ## Append to lists for plotting
    mag_grid_dict['times'].append(t_sim)
    mag_grid_dict['dBh_lat_avg'].append(mag_grid['dBh_lat_avg'])
#print(np.shape(mag_grid_dict['dBh']))
## (601, 360, 171)

## Compute time derivative of dBh.
mag_grid_dict['dBhdt'] = []
for i, timestamp in enumerate(mag_grid_dict['times']):
    if i == 0:
        ## Compute dBh/dt using forward difference
        dBh = -mag_grid_dict['dBh'][i+2] + \
              4*mag_grid_dict['dBh'][i+1] - \
              3*mag_grid_dict['dBh'][i]
        dT = mag_grid_dict['times'][i+2] - mag_grid_dict['times'][i]
        mag_grid_dict['dBhdt'].append(dBh/dT)
    elif i == np.shape(mag_grid_dict['times'])[0]-1:
        ## Compute dBh/dt using backwards difference
        dBh = 3*mag_grid_dict['dBh'][i] - \
              4*mag_grid_dict['dBh'][i-1] + \
              mag_grid_dict['dBh'][i-2]
        dT = mag_grid_dict['times'][i] - mag_grid_dict['times'][i-2]
        mag_grid_dict['dBhdt'].append(dBh/dT)
    else:
        ## Compute dBh/dt using central difference
        dBh = mag_grid_dict['dBh'][i+1] - mag_grid_dict['dBh'][i-1]
        dT = mag_grid_dict['times'][i+1] - mag_grid_dict['times'][i-1]
        mag_grid_dict['dBhdt'].append(dBh/dT)
        
## Compute latitude averages of dBh/dt
mag_grid_dict['dBhdt_lat_avg'] = []
for i,timestamp in enumerate(mag_grid_dict['times']):
    mag_grid_dict['dBhdt_lat_avg'].append(np.mean(mag_grid_dict['dBhdt'][i],
                                                  axis=0, dtype=np.float64))

## Compute dBh/dt at a point.
## Arbitrarily select Fairbanks, AK [65,-148]
lat = 45
lon = -75
## this math will break if the maggrid resolution changes
i_lat = lat + 85
if lon < 0: i_lon = lon + 360
mag_grid_dict['dBh_point'] = []
for i,timestamp in enumerate(mag_grid_dict['times']):
    mag_grid_dict['dBh_point'].append(mag_grid_dict['dBh'][i][i_lon][i_lat])
mag_grid_dict['dBhdt_point'] = []
for i, timestamp in enumerate(mag_grid_dict['times']):
    if i == 0:
        ## Compute dBh/dt using forward difference
        dBh = -mag_grid_dict['dBh_point'][i+2] + \
              4*mag_grid_dict['dBh_point'][i+1] - \
              3*mag_grid_dict['dBh_point'][i]
        dT = mag_grid_dict['times'][i+2] - mag_grid_dict['times'][i]
        mag_grid_dict['dBhdt_point'].append(dBh/dT)
    elif i == np.shape(mag_grid_dict['times'])[0]-1:
        ## Compute dBh/dt using backwards difference
        dBh = 3*mag_grid_dict['dBh_point'][i] - \
              4*mag_grid_dict['dBh_point'][i-1] + \
              mag_grid_dict['dBh_point'][i-2]
        dT = mag_grid_dict['times'][i] - mag_grid_dict['times'][i-2]
        mag_grid_dict['dBhdt_point'].append(dBh/dT)
    else:
        ## Compute dBh/dt using central difference
        dBh = mag_grid_dict['dBh_point'][i+1] - mag_grid_dict['dBh_point'][i-1]
        dT = mag_grid_dict['times'][i+1] - mag_grid_dict['times'][i-1]
        mag_grid_dict['dBhdt_point'].append(dBh/dT)

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

#Make new folder for plots to be stored in
isExist = os.path.exists(directory.path + '/plots/')
if not isExist:

   # Create a new directory because it does not exist
   os.mkdirs(directory.path + '/plots/')

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

plt.suptitle('Geo Indexes - 5 amu IBC Single Fluid')

plt.savefig('{}/gen_indexes'.format(directory.path),dpi=300) 
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
## dB_h / dt 
#Adapted from Pulkkinen et al. 2013, where B_x and B_y are B_n and B_e.
'''
    ## No times were saved in the file you gave me.
    ## If you want actual times, do this:
    #epoch = dt.datetime(YYYY, MM, DD,
    #                    HH, MM, SS, MSC) ## starttime of the run
    #t_sim = epoch + dt.timedelta(seconds=t_sim))
    
    #dt = np.array([x.total_seconds() for x in np.diff(geo_index['time'])])

plt.plot(log_file['time'],log_file['rho'],label = 'total')
#plt.plot(log_file['time'],log_file['rhosw'],label='sw')
#plt.plot(log_file['time'],log_file['rhoion'],label='ion')
plt.xticks(rotation=45)
plt.legend()


marker = ["8","o","s","H","X","D","*","d",">"]
## Plot Northern Hemisphere dBh
fig, ax = plt.subplots(figsize=(14,8))
m=0
for iLat in np.arange(85,171,10): ## start from Equator, do every 10 deg
        ax.plot(mag_grid_dict['times'],
                np.asarray(mag_grid_dict['dBh_lat_avg'])[:,iLat],
                label = 'Lat = {}'.format(iLat-85),marker=marker[m],
                markersize=12,markevery=50)
        m = m+1
ax.legend()
ax.set_xlabel('Simulation Time [s]',fontsize=16)
ax.set_ylabel('dBh [nT]',fontsize=16)
ax.set_title('dBh, Northern Hemisphere - {}'.format(run_name),fontsize=20)
plt.savefig('{}/plots/lat_avg_northern.png'.format(directory.path),dpi=300)
#plt.show()
plt.close()

## Plot Southern Hemisphere dBh
m=0
fig, ax = plt.subplots(figsize=(14,8))
for iLat in np.arange(0,86,10): ## end at Equator, do every 10 deg
        ax.plot(mag_grid_dict['times'],
                np.asarray(mag_grid_dict['dBh_lat_avg'])[:,iLat],
                label = 'Lat = {}'.format(iLat-85),marker=marker[m],
                markersize=12,markevery=50)
        m = m+1
ax.legend()
ax.set_xlabel('Simulation Time [s]',fontsize=16)
ax.set_ylabel('dBh [nT]',fontsize=16)
ax.set_title('dBh, Southern Hemisphere - {}'.format(run_name),fontsize=20)
plt.savefig('{}/plots/lat_avg_southern.png'.format(directory.path),dpi=300)
#plt.show()
plt.close()

## Plot Northern Hemisphere dBh/dt
m=0
fig, ax = plt.subplots(figsize=(14,8))
for iLat in np.arange(85,171,10): ## start from Equator, do every 10 deg
        ax.plot(mag_grid_dict['times'],
                np.asarray(mag_grid_dict['dBhdt_lat_avg'])[:,iLat],
                label = 'Lat = {}'.format(iLat-85),marker=marker[m],
                markersize=12,markevery=50)
ax.legend()
ax.set_xlabel('Simulation Time [s]',fontsize=16)
ax.set_ylabel('dBh/dt [nT/s]',fontsize=16)
ax.set_title('dBh/dt, Northern Hemisphere - {}'.format(run_name),fontsize=20)
plt.savefig('{}/plots/dbhdt_lat_avg_northern.png'.format(directory.path),dpi=300)
#plt.show()
plt.close()

#Plot dBh separated by latitude - Northern Hemisphere
m=0
fig = plt.figure(figsize=(14,20), constrained_layout=True)
gs0 = fig.add_gridspec(9,1, figure=fig)
plt.suptitle('dBh, Northern Hemisphere - {}'.format(run_name),fontsize=20)
for n in np.arange(0,9):  
    ax0 = fig.add_subplot(gs0[n]) 
    iLat = 10*n + 85 ## or something; I can't test my algebra
    ax0.plot(mag_grid_dict['times'],
                 np.asarray(mag_grid_dict['dBh_lat_avg'])[:,iLat],
                 label = 'Lat = {}'.format(iLat-85))
    ax0.legend()

ax0.set_xlabel('Simulation time [s]',fontsize=20)
ax0.set_ylabel('dBh [nT]', fontsize=20)
plt.savefig('{}/plots/lat_avg_N_separated_10eV.png'.format(directory.path),dpi=300)    
#plt.show()
plt.close()

#Plot dBh separated by latitude - Southern Hemisphere
m=0
fig = plt.figure(figsize=(14,20), constrained_layout=True)
gs1 = fig.add_gridspec(9,1, figure=fig)
plt.suptitle('dBh, Southern Hemisphere - {}'.format(run_name),fontsize=20)
for n in np.arange(0,9):  
    ax1 = fig.add_subplot(gs1[n]) 
    iLat = 10*n ## or something; I can't test my algebra
    ax1.plot(mag_grid_dict['times'],
                 np.asarray(mag_grid_dict['dBh_lat_avg'])[:,iLat], 
                 label = 'Lat = {}'.format(iLat-80))
    ax1.legend()

ax1.set_xlabel('Simulation time [s]',fontsize=20)
ax1.set_ylabel('dBh [nT]', fontsize=20)
plt.savefig('{}/plots/lat_avg_S_separated_10eV.png'.format(directory.path),dpi=300)    
#plt.show()
plt.close()
    
## Trying to do a 20-minute mask
## Using Dan's whiteboard scribble as reference
## Trying to preserve variable names too.
#tstart = mag_grid_dict['times'][0]
mask_minutes = 20
dt_seconds = mag_grid_dict['times'][2] - mag_grid_dict['times'][1]
mask_di = mask_minutes*60/dt_seconds # indices in one mask
#tstop = tstart + dt.timedelta(minutes = mask_minutes)
max20 = []
for i in range(int(np.shape(mag_grid_dict['times'])[0]/mask_di)):
    max20.append(max(mag_grid_dict['dBhdt_point']
                     [int(i*mask_di):int((i+1)*mask_di)]))     
'''
## Plot dBh at a point
fig, ax = plt.subplots(figsize=(14,8))
ax.plot(mag_grid_dict['times'],
        np.asarray(mag_grid_dict['dBh_point']))
ax.set_xlabel('Simulation Time [s]')
ax.set_ylabel('dBh [nT]')
ax.set_title('Location: [{},{}] - Ottawa (OTT)'.format(lat,lon),
             fontname = "Helvetica", fontsize = 20,fontweight='bold')
plt.savefig('{}/plots/dbh_point.png'.format(directory.path),dpi=300)
#plt.show()
plt.close()

## Plot dBh/dt at a point
fig, ax = plt.subplots(figsize=(14,8))
ax.plot(mag_grid_dict['times'],
        np.asarray(mag_grid_dict['dBhdt_point']))
ax.set_xlabel('Simulation Time [s]')
ax.set_ylabel('dBh/dt [nT/s]')
ax.set_title('Location: [{},{}] - Ottawa (OTT)'.format(lat,lon),
             fontname = "Helvetica", fontsize = 20,fontweight='bold')
plt.savefig('{}/dbhdt_point.png'.format(directory),dpi=300)
#plt.show()
plt.close()

## Plot mask of dBh/dt at a point
fig, ax = plt.subplots(figsize=(14,8))
ax.plot(max20)
ax.set_xlabel('Mask Point')
ax.set_ylabel('max(dBh/dt) [nT/s]')
ax.set_title('Location: [{},{}]'.format(lat,lon))
plt.savefig('{}/dbhdt_mask.png'.format(directory),dpi=300)
plt.show()
plt.close()
'''
