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


class ShellSlice(IdlFile):
    '''
    Shell slices are special MHD outputs where the domain is interpolated
    onto a spherical slice in 1, 2, or 3 dimensions.  Some examples
    include radial or azimuthal lines, spherical shells, or 3D wedges.

    The *Shell* class reads and handles these output types.
    '''

    def __init__(self, filename, radius, format='binary', *args, **kwargs):

        from spacepy.pybats import parse_filename_time

        IdlFile.__init__(self, filename, format=format, keep_case=False)

        # Extract time from file name:
        i_iter, runtime, time = parse_filename_time(self.attrs['file'])
        if 'time' not in self.attrs: self.attrs['time'] = time
        if 'iter' not in self.attrs: self.attrs['iter'] = i_iter

        ### Create some helper variables for plotting and calculations
        d2r = np.pi/180. # Convert degrees to radians.

        # Get grid spacing.  If npoints ==1, set to 1 to avoid math errors.
        self['r'] = dm.dmfilled((1,), fillval=radius, attrs={'units': 'R'})
        self.drad = 0
        self.dlon = (self['lon'][-1] - self['lon'][0])/max(self['grid'][0]-1,1)
        self.dlat = (self['lat'][-1] - self['lat'][0])/max(self['grid'][1]-1,1)

        self.dphi   = d2r*self.dlon
        self.dtheta = d2r*self.dlat

        # Get spherical, uniform grid in units of r_body/radians:
        self.lon, self.r, self.lat = np.meshgrid(
            np.array(self['lon']), np.array(self['r']), np.array(self['lat']))
        self.phi   = d2r*self.lon
        self.theta = d2r*(90-self.lat)

        #self['x']= self.r*np.sin(self.lat)*np.cos(self.lon)
        #self['y']= self.r*np.sin(self.lat)*np.sin(self.lon)
        #self['z'] = self.r*np.cos(self.lat)

    @calc_wrapper
    def calc_urad(self):
        '''
        Calculate radial velocity.
        '''

        ur = self['ux']*np.sin(self.theta)*np.cos(self.phi) + \
             self['uy']*np.sin(self.theta)*np.sin(self.phi) + \
             self['uz']*np.cos(self.theta)

        self['ur'] = dm.dmarray(ur, {'units':self['ux'].attrs['units']})

    @calc_wrapper
    def calc_radflux(self, var, conv=1000.*(100.0)**3):
        '''
        For variable *var*, calculate the radial flux of *var* through each
        grid point in the slice as self[var]*self['ur'].
        Resulting value stored as "var_rflx".
        '''

        #if var+'_rflx' in self:
        #    return

        # Make sure we have radial velocity.
        # if 'ur' not in self:
        self.calc_urad()

        # Calc flux:
        self[var+'_rflx'] = self[var] * self['ur'] * conv

    @calc_wrapper
    def calc_radflu(self, var):
        '''
        For variable *var*, calculate the radial fluence, or the
        spatially integrated radial flux through 2D surfaces of
        constant radius.

        Resulting variable stored as "var_rflu".  Result will be an array
        with one value for each radial distance within the object.
        '''

        # Need at least 2D in angle space:
        if self.dphi==0 or self.dtheta==0:
            raise ValueError('Fluence can only be calculated for 2D+ surfaces.')

        # Trim flux off of val name:
        if '_rflx' in var: var = var[:-5]

        # Convenience:
        flux = var + '_rflx'
        flu  = var + '_rflu'

        # Make sure flux exists:
        if flux not in self: self.calc_radflux(var)

        # Create output container, one point per radial distance:
        self[flu] = np.zeros( self['grid'][0] )

        # Integrate over all radii.
        # Units: convert R to km and cm-3 to km.
        for i, R in enumerate(self['r']):
            self[flu] = np.sum(
                (R* 6371.0)**2 * self[flux][i,:,:] * \
                np.sin(self.theta[i,:,:])  * \
                self.dtheta * self.dphi * 1000.**2 )

    def add_cont_shell(self, value, radius, irad=0, target=None, loc=111,
                       zlim=None, dolabel=True, add_cbar=False,
                       dofill=True, nlev=51, colat_max=90,
                       rotate=np.pi/2, dolog=False, direction=False, 
                       clabel=None,latticks=15., yticksize=14, extend='both',
                       **kwargs):
        '''
        For slices that cover a full hemisphere or more, create a polar
        plot of variable *value*.

        Extra keywords are sent to the matplotlib contour command.
        '''

        from numpy import pi
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.patches import Circle
        from matplotlib.colors import (LogNorm, Normalize, SymLogNorm)
        from matplotlib.ticker import (LogLocator, LogFormatter,
                                       LogFormatterMathtext, MultipleLocator)

        fig, ax = set_target(target, figsize=(10,10), loc=loc, polar=True)
        ax.set_title(f'{value} at {radius} Re')

        # Isolate the values we want to plot. Convert to numpy array
        # to avoid conflict between dmarray and contour plots.
        # Reduce dimensionality as needed.
        if self[value].ndim == 3:
            z = np.array(self[value][irad, :, :])
        else:
            z = np.array(self[value])

        # Get max/min if none given.
        if zlim is None:
            zlim = [0, 0]
            zlim[0] = z.min()
            zlim[1] = z.max()

            if abs(zlim[0]) > zlim[1]:
                zlim[1] = abs(zlim[0])

            if abs(zlim[1]) > abs(zlim[0]):
                zlim[0] = -zlim[1]
            #For log space, no negative zlimits.
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        # Create levels and set norm based on dolog.
        if dolog:  # Log space!
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            z=np.where(z>zlim[0], z, 1.01*zlim[0])
            norm=LogNorm()
            ticks=LogLocator()
            fmt=LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            norm=None
            ticks=None
            fmt=None

        # Select proper contour function based on fill/don't fill.
        if dofill:
            func = ax.contourf
        else:
            func = ax.contour

        # Plot result.  Rotate "rotate" radians to get sun in right spot.
        # Plot against colatitude to arrange results correctly.
        cnt = func(self.phi[irad,:,:]+rotate, 90-self.lat[irad,:,:],
                   np.array(z), levs, norm=norm,
                   extend=extend, **kwargs)

        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(cnt, ax=ax, ticks=ticks, format=fmt, shrink=.85)
            if clabel==None:
                clabel="{} ({})".format(value, self[value].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.

        # Adjust latitude
        ax.set_ylim( [0, colat_max] )

        # Adjust atitude and add better labels:
        ax.yaxis.set_major_locator(MultipleLocator(latticks))
        ax.set_ylim([0,colat_max])
        ax.set_yticklabels('')
        opts = {'size':yticksize, 'rotation':-45, 'ha':'center', 'va':'center'}
        for theta in np.arange(90-latticks, 90-colat_max, -latticks):
            txt = '{:02.0f}'.format(theta)+r'$^{\circ}$'
            ax.text(pi/4., 90.-theta, txt, color='w', weight='extra bold',**opts)
            ax.text(pi/4., 90.-theta, txt, color='k', weight='light', **opts)

        # Use MLT-type labels.
        lt_labels = ['Noon', '18', '00',   '06']
        xticks    = [     0, pi/2,   pi, 3*pi/2]
        xticks = np.array(xticks) + rotate

        # Apply x-labels:
        ax.set_xticks(xticks)
        ax.set_xticklabels(lt_labels)
        ax
        return fig, ax, cnt, cbar


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
radius = '30'
for filename in glob.glob((directory.path + '/' +'shl_mhd_4*.out*')):
    shell_file = ShellSlice(filename, 3.0)

for filename in glob.glob((directory.path + '/' + 'log*.log')):
    log_file = bats.LogFile(filename)
    
for filename in glob.glob((directory.path + '/' +'geo*.log')):
    geo_index = bats.GeoIndexFile(filename)

nFrame = shell_file.attrs['nframe']
shell_file_dict = {'nFrame': [], 'fluence': [], 'ur': [], 'temp_hp': []}

#Make new folder for plots to be stored in
isExist = os.path.exists(directory.path + '/plots/')
if not isExist:

   # Create a new directory because it does not exist
   os.mkdir(directory.path + '/plots/')
def shell_plot_dens():
    for iFrame in range(0,nFrame):
        shell_file.switch_frame(iFrame)
        shell_file.calc_urad()
        shell_file.calc_radflux('rho')
        plt.figure()
        fig = shell_file.add_cont_shell('rho_rflx', 2.7, add_cbar=True,
                                        clabel='Density flux (Mp/cc)',cmap='bwr')
        shell_file.calc_radflu('rho')
        shell_file_dict['fluence'].append(shell_file['rho_rflu'])
        shell_file_dict['ur'].append(shell_file['ur'])
        t_now = shell_file.attrs['runtime']
        t_now_hr = t_now/3600
        plt.title('Density Flux, {} hr simulation time, {}'.format(
            t_now_hr,run_name))
        plt.savefig('{}/plots/shl_slice_{}'.format(directory.path,iFrame), dpi=300)
        plt.close()
    
    return fig

def plot_fluence_log():

    #Plot the fluence across the runtime in log scale
    shell_file.calc_radflu('rho')
    shell_file_dict['fluence'].append(shell_file['rho_rflu'])
    plt.figure(figsize=[8,4])
    plt.plot(shell_file.attrs['runtimes'],shell_file_dict['fluence'])
    plt.yscale('log')
    plt.xlabel('Simulation time, [s]')
    plt.ylabel('Fluence, [particles*m/s]')
    plt.title('Fluence - {}'.format(run_name))
    plt.savefig('{}/plots/fluence_time_series_IBC_log_{}'.format(directory.path,
                                                                 radius),dpi=300)
    plt.close()

    return 

def plot_fluence():
    #Plot the fluence across the runtime with liear scale
    for iFrame in range(0,nFrame):
        shell_file.switch_frame(iFrame)    
        shell_file.calc_radflu('rho')
        shell_file_dict['fluence'].append(shell_file['rho_rflu'])
    plt.figure(figsize=[8,4])
    plt.plot(shell_file.attrs['runtimes'],shell_file_dict['fluence'])
    plt.xlabel('Simulation time, [s]')
    plt.ylabel('Fluence, [particles*m/s]')
    plt.title('Fluence - {}'.format(run_name))
    plt.savefig('{}/plots/fluence_time_series_IBC_{}'.format(directory.path,
                                                             radius),dpi=300)
    plt.close()
    
    return 

def plot_density_flux_cont():
    #Plotting the density, parallel velocity, and temperature for each run
    for iFrame in range(0,nFrame):
        shell_file.switch_frame(iFrame)
        shell_file.calc_radflux('rho')
    
        plt.figure()
        shell_file.add_cont_shell('rho_rflx', 3.0, add_cbar=True,
                                        clabel='Density [Mp/cc]',cmap='bwr')
        t_now = shell_file.attrs['runtime']
        t_now_hr = t_now/3600
        plt.title('Density Flux, {} hr simulation time, {}'.format(
            t_now_hr,run_name))
        plt.savefig('{}/plots/dens_flux_cont_{}_{}'.format(directory.path,iFrame,
                                                         radius), dpi=300)
        plt.close()

def plot_density_cont():
    #Plotting the density, parallel velocity, and temperature for each run
    for iFrame in range(0,nFrame):
        shell_file.switch_frame(iFrame)
    
        plt.figure()
        fig = shell_file.add_cont_shell('rho', 3.0, add_cbar=True,
                                        clabel='Density [Mp/cc]',cmap='YlGn',
                                        zlim=[shell_file['rho'][iFrame].min(),
                                              shell_file['rho'][iFrame].max()])
        t_now = shell_file.attrs['runtime']
        t_now_hr = t_now/3600
        plt.title('Density, {} hr simulation time, {}'.format(
            t_now_hr,run_name))
        plt.savefig('{}/plots/dens_contour_{}_{}'.format(directory.path,iFrame,
                                                         radius), dpi=300)
        plt.close()

def plot_ur_cont():
    for iFrame in range(nFrame):
        shell_file.switch_frame(iFrame)
        plt.figure()
        fig = shell_file.add_cont_shell('ur', 3.0, add_cbar=True,
                                        clabel='Radial Velocity, [km/s]',
                                        cmap='bwr',
                                        zlim=[-1*shell_file_dict['ur'][iFrame].max(),
                                              shell_file_dict['ur'][iFrame].max()])
        t_now = shell_file.attrs['runtime']
        t_now_hr = t_now/3600
        plt.title('Radial Velocity [km/s], {} hr simulation time, {}'.format(
            t_now_hr,run_name))
        plt.savefig('{}/plots/v_para_contour_{}_{}'.format(directory.path,iFrame,
                                                           radius), dpi=300)
        plt.close()

def plot_temp_cont():
    for iFrame in range(0,nFrame):
        shell_file.switch_frame(iFrame)
        shell_file['temp_hp'] = (shell_file_dict['ur'][iFrame]/14.1)**2
        shell_file_dict['temp_hp'].append(shell_file['temp_hp'])
        print("Temp Hp appended at frame: {}".format(iFrame))
    
        plt.figure()
        fig = shell_file.add_cont_shell('temp_hp', 3.0, add_cbar=True,
                                        clabel='Temperature, [eV]',cmap='cividis',
                                        zlim=[shell_file_dict['temp_hp'][iFrame].min(),
                                              shell_file_dict['temp_hp'][iFrame].max()])
        t_now = shell_file.attrs['runtime']
        t_now_hr = t_now/3600
        plt.title('Temperature, {} hr simulation time, {}'.format(t_now_hr,
                                                                run_name))
        plt.savefig('{}/plots/temp_contour_{}_{}'.format(directory.path,iFrame,
                                                         radius), dpi=300)
        plt.close()
        
def plot_basic_run_parameters():
    grid = gridspec.GridSpec(3,1)
    fig = plt.figure(figsize=(6,8),constrained_layout=True)
    ax3, ax4, ax5 = fig.add_subplot(grid[0,0]),\
            fig.add_subplot(grid[1,0]),\
                fig.add_subplot(grid[2,0])
    
    geo_index.fetch_obs_ae()
    geo_index.fetch_obs_kp()
    
    ax3.plot(geo_index['time'],geo_index['Kp'])
    ax3.set_ylabel('Kp')
    ax3.tick_params(labelbottom=False)

    ax4.plot(geo_index['time'],geo_index['AU'], label='AU',marker='o', 
             markersize=6, markevery=20)
    ax4.plot(geo_index['time'],geo_index['AL'], label='AL',marker='*', 
             markersize=6, markevery=20)
    ax4.set_ylabel('AE')
    ax4.tick_params(labelbottom=False)
    ax4.legend()

    ax5.plot(log_file['time'],log_file['dst'])
    ax5.set_xlabel('time')
    ax5.set_ylabel('Dst, [nT]')
    ax5.tick_params(rotation=30)

    plt.suptitle('Indexes: {}'.format(directory.runname))

    plt.savefig('{}/plots/gen_indexes'.format(directory.path),dpi=300) 
    #plt.show()
    plt.close()
    
    
    fig.autofmt_xdate(rotation=45)
    plt.savefig('{}/plots/basic_run_params_{}'.format(directory.path,radius), dpi=300)
    
def plot_outflow_sf_deliver():
    plot_basic_run_parameters()
    plot_density_cont()
    plot_density_flux_cont()
    plot_fluence()
    
    return

def plot_outflow_fluence_grouped():
    fluence_dict = {'fluence_v1': [], 'fluence_v2':[], 'fluence_v3':[], 
                    'fluence_v4':[], 'fluence_nominal':[]}
    
    path_v1 = '/Volumes/hard_drive_mac/run_data/run_10eV_IBC_sf/GM/shl_mhd_4_n00008000_00951259.outs'
    path_v2 = '/Volumes/hard_drive_mac/run_data/run_50eV_IBC_sf/GM/shl_mhd_4_n00008000_00951309.outs'
    path_v3 = '/Volumes/hard_drive_mac/run_data/run_75eV_IBC_sf/GM/shl_mhd_4_n00008000_00951343.outs'
    path_v4 = '/Volumes/hard_drive_mac/run_data/run_150eV_IBC_sf/GM/shl_mhd_4_n00008000_00951193.outs'
    path_nominal = '/Volumes/hard_drive_mac/run_data/run_nominal_IBC_sf/GM/shl_mhd_4_n00008000_00951401.outs'

    IBC_v1 = ShellSlice(path_v1, 3.0)
    IBC_v2 = ShellSlice(path_v2, 3.0)
    IBC_v3 = ShellSlice(path_v3, 3.0)
    IBC_v4 = ShellSlice(path_v4, 3.0)
    IBC_nominal = ShellSlice(path_nominal, 3.0)

    for iFrame in range(0,21):
        IBC_v1.switch_frame(iFrame)    
        IBC_v2.switch_frame(iFrame)
        IBC_v3.switch_frame(iFrame)
        IBC_v4.switch_frame(iFrame)
        IBC_nominal.switch_frame(iFrame)
        
        IBC_v1.calc_radflu('rho')
        IBC_v2.calc_radflu('rho')
        IBC_v3.calc_radflu('rho')
        IBC_v4.calc_radflu('rho')
        IBC_nominal.calc_radflu('rho')
        
        fluence_dict['fluence_v1'].append(IBC_v1['rho_rflu'])
        fluence_dict['fluence_v2'].append(IBC_v2['rho_rflu'])
        fluence_dict['fluence_v3'].append(IBC_v3['rho_rflu'])
        fluence_dict['fluence_v4'].append(IBC_v4['rho_rflu'])
        fluence_dict['fluence_nominal'].append(IBC_nominal['rho_rflu'])
    
    t_sim_hr = IBC_v3.attrs['runtimes']/3600
    
    fig = plt.figure(figsize=[12,10])
    fig, ax = plt.subplots()
    plt.plot(t_sim_hr,fluence_dict['fluence_v1'],'b', 
             marker = 'o', markersize = 6, label='10 eV')
    plt.plot(t_sim_hr,fluence_dict['fluence_v2'],'g', 
             marker = '*', markersize = 6, label='50 eV')
    plt.plot(t_sim_hr,fluence_dict['fluence_v3'],'r', 
             marker = '+', markersize = 6, label='75 eV')
    plt.plot(t_sim_hr,fluence_dict['fluence_v4'],'c', 
             marker = 'x', markersize = 6, label='150 eV')
    plt.plot(t_sim_hr,fluence_dict['fluence_nominal'],
             'k-', markersize = 8, label='Nominal Conditions')
    
    line1 = 0
    line2 = 5
    line3 = 10
    
    plt.xlabel('Simulation time, t, [hours]')
    plt.ylabel('Fluence, [particles*m/s]')
    plt.text(2.5, -.5e26,'+3 nT')
    plt.text(6.5, -.5e26,'-10 nT')
    ax.set_xlim(0,10)
    ax.axvspan(line1, line2, alpha=0.1, color='green')
    ax.axvspan(line2, line3, alpha=0.1, color='red')
    ax.grid()
    #plt.title('Fluence - Parallel Velocity IBCs 3.0 Re')
    plt.legend()
    
    return
    
    

'''
print(is_multi)

if is_multi == True:
    for iFrame in range(0,shell_file.attrs['nframe']):
        shell_file.switch_frame(iFrame)    
        grid = gridspec.GridSpec(1,2)
        fig = plt.figure(figsize=(12,8))
        fig.subplots_adjust(wspace = 1,hspace = 1)
        ax1, ax2 = fig.add_subplot(grid[0,0], polar=True),\
            fig.add_subplot(grid[0,1], polar=True)
        
        ax1 = shell_file.add_cont_shell('hprho', 3.0, add_cbar=True,
                            target=ax1,zlim=[-shell_file['hprho'][iFrame].max(),
                            shell_file['hprho'][iFrame].max()])
        ax2 = shell_file.add_cont_shell('oprho', 3.0, add_cbar=True,
                            target = ax2, zlim=[-shell_file['oprho'][iFrame].max(),
                                 shell_file['oprho'][iFrame].max()])
        #plt.suptitle('{} ,Multifluid Outflow Plots, 3 Re'.format(
        #    directory.runname, radius))
        fig.tight_layout()
        plt.savefig('{}/plots/density_mf_{}_{}'.format(directory.path,iFrame,
                                                         radius), dpi=300)
        plt.close()
         
    for iFrame in range(shell_file.attrs['nframe']):
        shell_file.switch_frame(iFrame)
        shell_file.calc_radflux('hprho')

        shell_file['hprho_rflx'] = shell_file['hprho'] * shell_file['ur'] *\
            1000.*(100.0)**3

        shell_file['oprho_rflx'] = (shell_file['oprho'] * shell_file['ur'] *\
            1000.*(100.0)**3)/16
        
        t_now = shell_file.attrs['runtime']
        t_now_hr = t_now/3600
            
        grid = gridspec.GridSpec(1,3)
        fig = plt.figure(figsize=(24,8))
        fig.subplots_adjust(wspace = 1,hspace = 2)
        ax1, ax2, ax3 = fig.add_subplot(grid[0,0], polar=True),\
            fig.add_subplot(grid[0,1], polar=True), \
                fig.add_subplot(grid[0,2], polar=True)
        
        ax1 = shell_file.add_cont_shell('hprho_rflx', 3.0, add_cbar=True,
                            target=ax1,zlim=[-1e13,1e13],cmap='bwr')
        ax2 = shell_file.add_cont_shell('oprho_rflx', 3.0, add_cbar=True,
                            target = ax2,zlim=[-1e13,1e13],cmap='bwr')
        ax3 = shell_file.add_cont_shell('rho_rflx', 3.0, add_cbar=True,
                            target=ax3,zlim=[-1e11,1e11],cmap='bwr')
        plt.suptitle('{} ,Multifluid Outflow Plots, 3 Re - {} hours'.format(
            directory.runname, t_now_hr))
        plt.savefig('{}/plots/shell_slice_mf_{}_{}'.format(directory.path,iFrame,
                                                         radius), dpi=300)
        plt.close()
'''