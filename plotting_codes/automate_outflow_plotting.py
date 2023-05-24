#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 15:04:22 2023

@author: kdoubles
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from matplotlib import gridspec
import argparse as ap
import os, glob

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

if input("Is this a multifluid run? (y/n) ") != "y":
    is_multi = False

global file_location
directory = arg_parse_file_location()
run_name = directory.runname

## Mag Grid File
for filename in glob.glob((directory.path + '/' +'mag_grid*.outs')):
    mag_grid = bats.MagGridFile(filename)
    
## Geo Index File
for filename in glob.glob((directory.path + '/' +'geo*.log')):
    geo_index = bats.GeoIndexFile(filename)

## Log File
for filename in glob.glob((directory.path + '/' +'log*.log')):
    log_file = bats.BatsLog(filename)

for filename in glob.glob((directory.path + '/' + 'y*.outs')):
    y_2d = bats.Bats2d(filename)
    y_idl = bats.IdlFile(filename)
    
for filename in glob.glob((directory.path + '/' + 'z*.outs')):
    z_2d = bats.Bats2d(filename)
    
nFrame = mag_grid.attrs['nframe']
print('nFrame: {}'.format(nFrame))

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
    
    epoch = dt.datetime(2000, 1, 1, 0, 0, 0, 0)
    times = epoch + dt.timedelta(seconds=t_sim)
    
    if iFrame % 100 == 0: # avoid screen barf
        print('iFrame: {}'.format(iFrame))
        print('tSimulation: {}'.format(t_sim))
        
    ## Calc dBh
    mag_grid.calc_h()
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

#Make new folder for plots to be stored in
isExist = os.path.exists(directory.path + '/plots/')
if not isExist:

   # Create a new directory because it does not exist
   os.mkdir(directory.path + '/plots/')

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

plt.savefig('{}/plots/gen_indexes'.format(directory.path),dpi=300) 
#plt.show()
plt.close()

## plot dBn with rho
plt.figure()
plt.plot(log_file['time'],log_file['rho'],label = 'total')
#plt.plot(log_file['time'],log_file['rhosw'],label='sw')
#plt.plot(log_file['time'],log_file['rhoion'],label='ion')
plt.xticks(rotation=45)
plt.legend()
plt.close()

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

y_2d['r']  = np.sqrt(y_2d['x']**2+y_2d['z']**2)
y_2d['phi'] = 0
y_2d['theta'] = np.cos(y_2d['z']/np.sqrt(y_2d['x']**2+y_2d['z']**2))**(-1)
'''
Calculate radial velocity.


ur = y_idl['Ux']*np.sin(y_idl['theta'])*np.cos(y_idl['phi']) + \
   y_idl['Uy']*np.sin(y_idl['theta'])*np.sin(y_idl['phi']) + \
   y_idl['Uz']*np.cos(y_idl['theta'])

y_idl['ur'] = dm.dmarray(ur, {'units':y_idl['Ux'].attrs['units']})

print(y_idl['r'])

#Stream scatter for velocity field

fig = plt.figure(figsize=[10,4])
y_2d.switch_frame(600)
fig, ax, cont, cbar = y_2d.add_stream_scatter('ux', 'uz',target = fig, loc = 121,\
                                xlim=[-10,10], ylim=[-10,10],colors='Gray')

#mhd_cont.add_b_magsphere(target=ax,colors='Black',DoLast=False)

plt.xlabel(r'$R_e$')
plt.ylabel(r'$R_e$')
plt.title('Y=0, Streamline Velocity with Radial Velocity', fontsize=10)
'''
'''
## These additional plots are produced when prompted question is answered 'y'.
if is_multi == True:
    #Plot all the .out files in the y=0 or z=0 outputs.
    grid=gridspec.GridSpec(2,2)

    for i in range(y_2d.attrs['nframe']):
        y_2d.switch_frame(i)
        fig, ax = plt.subplots(2,2,sharey=True, sharex=True, figsize=[10,10])
        
        #fig = mhd_cont.add_contour('x', 'z', 'oprho', target = fig, loc = 121, dolog = True,\
        #                                add_cbar=True, xlim=[-50,20], ylim=[-30,30],cmap='viridis')
        t_now = y_2d.attrs['time']
        fig.suptitle('{}th iteration - {} '.format(i,t_now))
        ax[0,0] = y_2d.add_contour('x', 'z', 'opp', target=ax[0,0],
                                   xlim=[-20,20],ylim=[-20,20],add_cbar=True,
                                   title='Op Pressure')
        ax[0,1] = y_2d.add_contour('x', 'z', 'p', target=ax[0,1],
                                   add_cbar=True,clabel=None, xlim=[-20,20],
                                   ylim=[-20,20],title='Pressure')
        ax[1,0] = y_2d.add_contour('x', 'z', 'oprho', target=ax[1,0],
                                   add_cbar=True,clabel=None,
                                   xlim=[-20,20],ylim=[-20,20],
                                   zlim=[0,30],title='Op Density')
        ax[1,1] = y_2d.add_contour('x', 'z', 'rho', target=ax[1,1],
                                   add_cbar=True,clabel=None,xlim=[-20,20],
                                   ylim=[-20,20],zlim=[0,30],title='Density')

if is_multi == False:
    for i in range(y_2d.attrs['nframe']):
        y_2d.switch_frame(i)
        fig, ax,  cont, cbar = y_2d.add_contour('x', 'z', 'rho',xlim=[-20,20],
                                                ylim=[-20,20],add_cbar=True,
                                   title='Density, {}'.format(times))
        plt.savefig('{}/plots/dens_contour_{}'.format(directory.path,directory.runname),dpi=300)    
        plt.close()
'''
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
