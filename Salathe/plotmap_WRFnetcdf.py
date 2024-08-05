#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: esalathe

Plot a map of WRF data
"""

# Import numpy and pyplot routines

from numpy import (
        linspace,array, log,exp,sin,cos,sqrt, pi,e, 
        zeros, ones, amin,amax, argmax, arange
        )
from matplotlib.pyplot import plot, figure, show

import numpy as np
import matplotlib.pyplot as plt

# Define the map plotting routine. uses the cartopy package.

def WRFplot(plotvar,lats,lons,  vmin,vmax, maptitle,varname, ColMap, smflg=1, domain='auto'):
    from matplotlib.pyplot import colorbar, axes, pcolormesh, colorbar, title

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.cm import get_cmap
    import cartopy
    import cartopy.crs as crs
    from cartopy.feature import NaturalEarthFeature

    # Set the cartopy mapping object for the WRF domain
    #  (Taken from wrf getcartopy and cartopy_xlim)
    cart_proj = cartopy.crs.LambertConformal(
            central_longitude=-121.0, 
            central_latitude=45.665584564208984, 
            false_easting=0.0, 
            false_northing=0.0, 
            standard_parallels=[30.,60.], 
            globe=None, 
            cutoff=-30)
    ax = axes(projection=cart_proj)

    # Set map limits based on domain (ideally use lat-lon and transform...)
    if domain=='pnw02':
        ax.set_xlim([-875806.9669240027, 1056192.549175313])
        ax.set_ylim([-733768.6404772081, 730230.3670079684])
    elif domain=='pnw01':
        ax.set_xlim([-3.7e6, 1.6e6])
        ax.set_ylim([-2.15e6, 2.3e6])
    elif domain=='west02':
        ax.set_xlim([-8.8e5, 1.2e6])
        ax.set_ylim([-1.58e6, 7.3e5])

    # Use field min/max if -999
    if vmin==-999: vmin=plotvar.min()
    if vmax==-999: vmax=plotvar.max()        
    
    # Color in the data on the map with smoothing

    if smflg == 1:
        smooth = 'gouraud'
    else:
        smooth = 'nearest'

    pcolormesh(lons,lats,
                   plotvar,vmin=vmin,vmax=vmax,
                   transform=crs.PlateCarree(),
                   shading=smooth,
                   cmap=get_cmap(ColMap)
                   )
    
    # Add a color bar
    cbar=colorbar(ax=ax, shrink=.6)#, orientation='horizontal')
    cbar.set_label(varname)


    # Add contour lines
    # plt.contour(lons,lats,plotvar, vmin=vmin,vmax=vmax, colors='gray', linewidths=0.5, transform=crs.PlateCarree())


    # Download and add the states and coastlines
    states = NaturalEarthFeature(category='cultural', scale='50m',
                             facecolor='none',
                             name='admin_1_states_provinces_lines')
    ax.add_feature(states, linewidth=.5, edgecolor='black')
    borders = NaturalEarthFeature(category='cultural', scale='50m',
                             facecolor='none',
                             name='admin_0_boundary_lines_land')
    ax.add_feature(borders, linewidth=.75, edgecolor='black')
    ax.coastlines('50m', linewidth=0.8)
    

    
    # Add gridlines
    #ax.gridlines(color='black', linestyle='dotted')
    
    # Add a title
    title(maptitle)


"""
Once everything is defined, we can do some tasks with the WRF data. Here we just select one year and plot it
"""

# 1) open netcdf file. Neet to install the netCDF4 module
from netCDF4 import Dataset
ncFile = Dataset("ccsm4-wrf_1970-2099_T2MAX_extr.nc", "r", format="NETCDF4")

# 2) read in the data to plot. Note that the [:] converts the netCDF object to a numpy array.
T90=ncFile.variables["T2MAX90"][:]
lats = ncFile.variables["XLAT"][:]
lons = ncFile.variables["XLONG"][:]
yr0 = 1970 # starting year in the dataset (ideally this would be stored in the netCDF file)

# 3) Select a year to plot and pull just this one time from the full data array

year = 1990
iyr = 1990 - yr0
plotvar=T90[iyr,:,:]-273.15

# 4) Plot the map
Tmin=10 # min for contour lines
Tmax=40 # max for contour lines

# Use function defined above to draw the map

figure(figsize=(15,10))
WRFplot(plotvar,lats,lons, Tmin,Tmax, "CCSM4-WRF 90th Percentile Tamx "+str(year), "TMAX90 in °C", "RdYlBu_r")
show()
