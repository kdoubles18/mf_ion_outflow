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

path_outs = '/Volumes/hard_drive_mac/run_data/run_nominal_IBC_sf/GM/y=0_mhd_1_n00008000_00951401.outs'
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
mhd_cont.switch_frame(0)
fig, ax = mhd_cont.add_grid_plot()

#%%
'''
Set up to plot y=0 cut. Change third value to plot different attribute in the
#add_contour function. Create a new figure showing the contour of the cuts.
'''
fig = plt.figure(figsize=[8,4])
mhd_cont.switch_frame(120)
fig, ax, cont, cbar = mhd_cont.add_contour('x', 'z', 'rho', target=fig,
                                           loc = 121 ,dolog = False, 
                                           xlim=[-10,10], ylim=[-10,10],
                                           add_cbar=True)
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
#%%
#Stream scatter for velocity field

fig = plt.figure(figsize=[10,4])
mhd_cont.switch_frame(120)

#fig, ax, cont, cbar = mhd_cont.add_contour('x', 'z', 'ur', target=fig,
#                                           loc = 121,dolog = True, 
#                                           xlim=[-10,10], ylim=[-10,10],
#                                           add_cbar=True,cmap='coolwarm')

fig, ax, cont, cbar = mhd_cont.add_stream_scatter('ux', 'uz',target = fig, loc = 121,
                            xlim=[-5,5], ylim=[-5,5],colors='Gray')

#mhd_cont.add_b_magsphere(target=ax,colors='Black',DoLast=False)
mhd_cont.add_pcolor('x','z','ur',target=ax,add_cbar=True,cmap='bwr',zlim=[-85,85])
plt.Circle((0,0),3.0,fill=False)


plt.xlabel(r'$R_e$')
plt.ylabel(r'$R_e$')
plt.title('Y=0, Streamline Velocity with Radial Velocity', fontsize=10)

#%%
mhd_cont.calc_jxb()
mhd_cont.calc_gradP()

#%%
for i in np.arange(100,mhd_cont.attrs['nframe'],5):
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

path_outs = '/Volumes/hard_drive_mac/run_data/run_nominal_IBC_sf/GM/z=0_mhd_2_n00008000_00951401.outs'

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
mhd_cont.switch_frame(0)
fig, ax = mhd_cont.add_grid_plot()

#%%
'''
Set up to plot y=0 cut. Change third value to plot different attribute in the
#add_contour function. Create a new figure showing the contour of the cuts.
'''
fig = plt.figure(figsize=[8,4])
mhd_cont.switch_frame(120)
fig, ax, cont, cbar = mhd_cont.add_contour('x', 'y', 'rho', target=fig,
                                           loc = 121 ,dolog = False, 
                                           xlim=[-10,10], ylim=[-10,10],
                                           add_cbar=True)
plt.title('Y=0, Density, Single Fluid')


mhd_cont.add_b_magsphere(target=ax,colors='Black',DoLast=False)

fig.tight_layout()

#%%

mhd['r']  = np.sqrt(mhd['x']**2+mhd['y']**2)
mhd['phi'] = 0
mhd['theta'] = np.cos(mhd['y']/np.sqrt(mhd['x']**2+mhd['y']**2))**(-1)
'''
Calculate radial velocity.
'''

ur = mhd['Ux']*np.sin(mhd['theta'])*np.cos(mhd['phi']) + \
   mhd['Uy']*np.sin(mhd['theta'])*np.sin(mhd['phi']) + \
   mhd['Uz']*np.cos(mhd['theta'])

mhd_cont['ur'] = dm.dmarray(ur, {'units':mhd['Ux'].attrs['units']})
#%%
#Stream scatter for velocity field

fig = plt.figure(figsize=[10,4])
mhd_cont.switch_frame(120)

#fig, ax, cont, cbar = mhd_cont.add_contour('x', 'z', 'ur', target=fig,
#                                           loc = 121,dolog = True, 
#                                           xlim=[-10,10], ylim=[-10,10],
#                                           add_cbar=True,cmap='coolwarm')

fig, ax, cont, cbar = mhd_cont.add_stream_scatter('ux', 'uy',target = fig, loc = 121,
                            xlim=[-10,10], ylim=[-10,10],colors='Gray')

#mhd_cont.add_b_magsphere(target=ax,colors='Black',DoLast=False)
mhd_cont.add_pcolor('x','y','ur',target=ax,add_cbar=True,cmap='bwr',zlim=[-85,85])
plt.Circle((0,0),3.0,fill=False)


plt.xlabel(r'$R_e$')
plt.ylabel(r'$R_e$')
plt.title('Y=0, Streamline Velocity with Radial Velocity', fontsize=10)

#%%
mhd_cont.calc_jxb()
mhd_cont.calc_gradP()


