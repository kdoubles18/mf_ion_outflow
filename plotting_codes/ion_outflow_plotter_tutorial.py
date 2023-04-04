#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 11:45:44 2022

@author: kdoubles
"""

#imports
import numpy as np
import pickle
import matplotlib.pyplot as plt
import spacepy.plot as splot
from matplotlib import gridspec
import spacepy.datamodel as dm
from matplotlib.gridspec import GridSpec

from spacepy import pybats
from spacepy.pybats import bats

path_outs = '/Users/kdoubles/Downloads/drive-download-20230403T171817Z-001/y=0_mhd_1_n00000100_00321530.outs'
#%%
#Reads in the file so Python knows how to plot it
mhd = pybats.IdlFile(path_outs)

#Attrs shows what is present in the file.
print(mhd.keys())
print(mhd.attrs['time'])
#Bats2d allows you to plot different types of plots. 
mhd_cont = bats.Bats2d(path_outs)
mhd_cont.calc_b()

#Plot grid resolution of run
print(mhd_cont['rho'].mean())
mhd_cont.switch_frame(180)
fig, ax = mhd_cont.add_grid_plot()

#%%
'''
Set up to plot y=0 cut. Change third value to plot different attribute in the
#add_contour function. Create a new figure showing the contour of the cuts.
'''
fig = plt.figure(figsize=[8,4])
mhd_cont.switch_frame(180)
fig, ax, cont, cbar = mhd_cont.add_contour('x', 'z', 'rho', target=fig,loc = 121 ,dolog = True, xlim=[-10,10], ylim=[-10,10],add_cbar=True)
plt.title('Y=0, Density, Single Fluid')


mhd_cont.add_b_magsphere(target=ax,colors='Black',DoLast=False)

fig.tight_layout()

#%%

mhd['r']  = np.sqrt(mhd['x']**2+mhd['z']**2)
mhd['phi'] = 0
mhd['theta'] = np.cos(mhd['z']/np.sqrt(mhd['x']**2+mhd['z']**2))**(-1)
'''
Calculate radial velocity.
'''

ur = mhd['Ux']*np.sin(mhd['theta'])*np.cos(mhd['phi']) + \
   mhd['Uy']*np.sin(mhd['theta'])*np.sin(mhd['phi']) + \
   mhd['Uz']*np.cos(mhd['theta'])

mhd_cont['ur'] = dm.dmarray(ur, {'units':mhd['Ux'].attrs['units']})

print(mhd['r'])
#%%
#Stream scatter for velocity field

fig = plt.figure(figsize=[10,4])
mhd_cont.switch_frame(170)
fig, ax, cont, cbar = mhd_cont.add_stream_scatter('ux', 'uz',target = fig, loc = 121,\
                                xlim=[-10,10], ylim=[-10,10],colors='Gray')

#mhd_cont.add_b_magsphere(target=ax,colors='Black',DoLast=False)

plt.xlabel(r'$R_e$')
plt.ylabel(r'$R_e$')
plt.title('Y=0, Streamline Velocity with Radial Velocity', fontsize=10)
#%%
mhd_cont.add_pcolor('x','z','ur',target=ax,add_cbar=True,cmap='coolwarm')
#mhd_cont.add_b_magsphere(target=ax,colors='Black',DoLast=False)

#%%
mhd_cont.calc_jxb()
mhd_cont.calc_gradP()

#%%
for i in np.arange(130,mhd_cont.attrs['nframe'],5):
    mhd_cont.switch_frame(i)
    fig, ax = plt.subplots(2,2,figsize=[10,10])
    t = mhd_cont.attrs['runtimes'][i]
    
    fig.suptitle(f"%i simulation time" %t)
    ax[0,0] = mhd_cont.add_contour('x', 'z', 'rho', target=ax[0,0],xlim=[-10,10],ylim=[-10,10],add_cbar=True,title='Density')
    ax[0,1] = mhd_cont.add_contour('x', 'z', 'rhosw', target=ax[0,1],add_cbar=True,clabel=None,xlim=[-10,10],ylim=[-10,10],title='Rho Sw')
    ax[1,0] = mhd_cont.add_contour('x', 'z', 'rhoion', target=ax[1,0],add_cbar=True,clabel=None,xlim=[-10,10],ylim=[-10,10],title='Rho Iono')
    ax[1,1] = mhd_cont.add_contour('x', 'z', 'p', target=ax[1,1],add_cbar=True,clabel=None,xlim=[-10,10],ylim=[-10,10],title='Pressure')
    fig.tight_layout()


#%%
#Plot all the .out files in the y=0 or z=0 outputs.
grid=gridspec.GridSpec(2,2)

for i in range(mhd.attrs['nframe']):
    mhd_cont.switch_frame(i)
    fig, ax = plt.subplots(2,2,sharey=True, sharex=True, figsize=[10,10])
    
    #fig = mhd_cont.add_contour('x', 'z', 'oprho', target = fig, loc = 121, dolog = True,\
    #                                add_cbar=True, xlim=[-50,20], ylim=[-30,30],cmap='viridis')
    t_now = mhd.attrs['time']
    fig.suptitle(f"%ith iteration - 3 min interval" %i)
    ax[0,0] = mhd_cont.add_contour('x', 'z', 'opp', target=ax[0,0],xlim=[-20,20],ylim=[-20,20],add_cbar=True,title='Op P')
    ax[0,1] = mhd_cont.add_contour('x', 'z', 'p', target=ax[0,1],add_cbar=True,clabel=None,xlim=[-20,20],ylim=[-20,20],title='P')
    ax[1,0] = mhd_cont.add_contour('x', 'z', 'oprho', target=ax[1,0],add_cbar=True,clabel=None,xlim=[-20,20],ylim=[-20,20],zlim=[0,30],title='Op Dens')
    ax[1,1] = mhd_cont.add_contour('x', 'z', 'rho', target=ax[1,1],add_cbar=True,clabel=None,xlim=[-20,20],ylim=[-20,20],zlim=[0,30],title='Dens')

#%%
'''
Plotting log file values & mag grid values
'''

mag_outs_log = bats.BatsLog('/Users/kdoubles/Data/test/log_n000010.log')

mag_outs_log['b'] = np.sqrt(mag_outs_log['bx']**2.0 + mag_outs_log['by']**2.0 + mag_outs_log['bz']**2.0)


print(mag_outs_log.keys())
print(mag_outs_log.values())


fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2,figsize=[10,10])
fig.suptitle('Run 01 - Log File Results')
ax0.plot(mag_outs_log['time'],mag_outs_log['rho'])
ax1.plot(mag_outs_log['time'],mag_outs_log['b'])
ax2.plot(mag_outs_log['time'],mag_outs_log['dst'])
ax3.plot(mag_outs_log['time'],mag_outs_log['p'])
ax0.set_title('Density')
ax1.set_title('B')
ax2.set_title('Dst')
ax3.set_title('Pressure')
ax0.set_ylabel('Number Density')
ax1.set_ylabel('nT')
ax2.set_ylabel('Dst')
ax3.set_ylabel('pressure')

plt.setp(ax0.get_xticklabels(),rotation=30,ha='right')
plt.setp(ax1.get_xticklabels(),rotation=30,ha='right')
plt.setp(ax2.get_xticklabels(),rotation=30,ha='right')
plt.setp(ax3.get_xticklabels(),rotation=30,ha='right')



               
fig.tight_layout()

#%%
from numpy.lib.recfunctions import append_fields

mag_grid = bats.MagFile('/Users/kdoubles/Data/outflow_runs/run_dates/Run_02_full/mag_grid_n00008000_00951393.outs')
#mag_grid_temp_50 = bats.MagGridFile('/Users/kdoubles/Data/run_results/SF_20230306/mag_grid_n00008000_00480948.outs')


lons = mag_grid['Lon']
lats = mag_grid['Lat']
dLon = lons[1]-lons[0]
dLat = lats[1]-lats[0]

xloc = (lons >= 1) & (lons <= 360)
yloc = (lats >= 10) & (lats <= 85)

mag_grid.calc_h()
mag_grid['dBh_arr'] = np.zeros((360,171,mag_grid.attrs['nframe']))

dBh = np.recarray((mag_grid['dBh_arr'],mag_grid.attrs['runtimes']))
#dBh = np.array(dBh.reshape((mag_grid['dBh'].shape[0],mag_grid['dBh'].shape[1],mag_grid.attrs['nframe'])))

#for i in range(mag_grid.attrs['nframe']):
    

#%%
log_file_name = '/Users/kdoubles/Downloads/drive-download-20230403T180214Z-001/log_n000010.log'
log_file = bats.BatsLog(log_file_name)

plt.figure()
plt.plot(log_file['runtime'],log_file['rho'],label='rho')
#plt.plot(log_file['runtime'],log_file['swrho'], label = 'swrho')
#plt.plot(log_file['runtime'],log_file['hprho'], label = 'hprho')
#plt.plot(log_file['runtime'],log_file['oprho'], label = 'oprho')
plt.legend()

plt.figure(figsize=(14,8))
plt.plot(log_file['runtime'],log_file['rhoflx_R=3.0'],label='R = 3.0')
plt.plot(log_file['runtime'],log_file['rhoflx_R=3.5'],label='R = 3.5')
plt.plot(log_file['runtime'],log_file['rhoflx_R=5.0'],label='R = 5.0')
plt.axvline(x = (log_file['runtime'][-1]/2), color = 'k', label = 'switch in IMF Bz')
plt.title('150 amu - 90% H+, 10% O+')
plt.legend()

#%%
rho_sw = mhd_cont['rhosw']
rho_sw_plt = np.arange(0,mhd_cont['rhosw'].shape[0],343)

rho_ion = mhd_cont['rhoion']
rho_ion_plt = np.arange(0,mhd_cont['rhoion'].shape[0],343)

rho_sw_verb = mhd_cont['rhosw'][0]*rho_sw_plt

plt.figure()
plt.plot(mhd_cont['rho'],mhd_cont['rhoion'])
plt.xlabel('Run Times')
plt.ylabel('')

#%%

pathname = '/Users/kdoubles/Downloads/drive-download-20230403T180214Z-001/'
filename = 'geoindex_n00000000.log'
geo_index = bats.GeoIndexFile(pathname+filename)

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

plt.savefig('{}/gen_'.format(pathname),dpi=300) 

#%%
log_file_name = '/Users/kdoubles/Downloads/drive-download-20230403T175725Z-001/log_n000010.log'
log_file = bats.BatsLog(log_file_name)

plt.figure()
plt.plot(log_file['runtime'],log_file['rho'],label='rho')
#plt.plot(log_file['runtime'],log_file['swrho'], label = 'swrho')
#plt.plot(log_file['runtime'],log_file['hprho'], label = 'hprho')
#plt.plot(log_file['runtime'],log_file['oprho'], label = 'oprho')
plt.legend()

plt.figure(figsize=(14,8))
plt.plot(log_file['runtime'],log_file['rhoflx_R=3.0'],label='R = 3.0')
plt.plot(log_file['runtime'],log_file['rhoflx_R=3.5'],label='R = 3.5')
plt.plot(log_file['runtime'],log_file['rhoflx_R=5.0'],label='R = 5.0')
plt.axvline(x = (log_file['runtime'][-1]/2), color = 'k', label = 'switch in IMF Bz')
plt.title('2 amu - 90% H+, 10% O+')
plt.legend()


