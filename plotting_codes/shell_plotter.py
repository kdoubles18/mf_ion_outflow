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
        if flu in self: return

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
                       rotate=np.pi/2, dolog=False, direction=False, clabel=None,
                       latticks=15., yticksize=14, extend='both', **kwargs):
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
                   extend=extend, **kwargs,cmap='bwr')

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

filename = '/Users/kdoubles/Data/Outflow_Runs/Run_Dates/MF_20221119_ta/shl_mhd_8_n00005000_00667147.outs'
SMG_run = ShellSlice(filename, 3)
SMG_run.calc_urad()
SMG_run.calc_radflux('rho')
SMG_run.switch_frame(44)
SMG_run.add_cont_shell('rho_rflx', 3, add_cbar=True,clabel='Density flux (Mp/cc)')
plt.title('Rho Flux at 3 Re')