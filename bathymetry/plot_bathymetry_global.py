#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:41:33 2019

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
import matplotlib
import cartopy.crs as ccrs
import seaborn as sns
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import math
from numba import jit
import cartopy.mpl.ticker as cticker
from matplotlib.lines import Line2D
import cmocean.cm as cm
import cartopy.feature as cfeature
from scipy.ndimage import gaussian_filter
from copy import copy

def mycolormap():
    import matplotlib.colors as mcolors
    # sample the colormaps that you want to use. Use 128 from each so we get 256
    # colors in total
    colors2 = cm.deep(np.linspace(0., 1, 128))
    colors1 = plt.cm.Greys(np.linspace(0, 1, 128))[100:102]

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    return mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

#%% load sediment sample site locations
ff = np.load('paleolocs_38Ma.npz')
paleolons = ff['longitudes'] 
paleolonstext = copy(paleolons)
paleolats = ff['latitudes']
sitenames = ff['names']
#%%

sns.set(context='paper',style="whitegrid",font="Arial")
fs = 19
font = {'size'   : fs}
matplotlib.rc('font', **font)

#variables
sp = 6
dd = 10
projection = ccrs.PlateCarree(0)
exte = [1, 360, -75, -10]
cmap = mycolormap()
vsbath = [-0.08,6]
rr = 1 # determines the resolution reduction for the plotting
contours = [-5e-5,-3e-5,-2e-5,-1e-5,0,1e-5, 2E-5, 3E-5, 5E-5]#[0,1.25E-5, 2E-5, 2.5E-5, 3E-5,3.75E-5, 5E-5]
#%% Load files
minlat = -84
maxlat = 84
minlat38 = -84
dirRead = '/Users/nooteboom/Documents/PhD/Eocene_POP/OUPUT/gridf/'
# Load 38MA
nc38 = Dataset(dirRead + 'kmt_tx0.1_POP_EO38_poleland.nc')
nc38grid = Dataset(dirRead + 'grid_coordinates_pop_tx0.1_38ma.nc')
lats38 = nc38grid['U_LAT_2D'][:]
lons38 = nc38grid['U_LON_2D'][:] + 25
minlon = np.min(lons38)
maxlon = np.max(lons38)
exte38 = [minlon, maxlon, minlat38, maxlat]
idx = np.where(np.logical_and(lats38[:,0]>=minlat38, lats38[:,0]<=maxlat))

bath38 = (-1*nc38['elevation'][:] / 1000)[idx[0]]
fH38 = (-2 * 7.2921 * 10**(-5) * np.sin(lats38[idx[0]]*np.pi/180) / bath38)
lats38 = lats38[idx[0]]
lons38 = lons38[idx[0]]
#%% start figure
latticks = [-100,-75,-50,-25, 0, 25, 50, 75, 100]#[-75,-65,-55,-45,-35,-25, 0, 25, 50, 75, 100]
fig = plt.figure(figsize=(16,9))

#%% subplot (a)
ax2 = plt.subplot(111, projection=projection)
plt.title('38Ma bathymetry', fontsize=fs)

g = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlocator = mticker.FixedLocator([-180,-90, -0, 90, 180])
g.xlabels_top = False
g.ylabels_right = False
g.xlabels_bottom = True
g.xlabels_left = True
g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.ylocator = mticker.FixedLocator(latticks)
ax2.set_extent(exte38, ccrs.PlateCarree())

im2 = plt.pcolormesh(lons38[::rr,::rr], lats38[::rr,::rr], bath38[::rr,::rr], cmap=cmap, 
                     vmin=vsbath[0], vmax=vsbath[1],
                     transform=ccrs.PlateCarree()
                     )
ax2.scatter(paleolons,paleolats, marker='+', s=1000, color='red', 
            zorder=3000)

tc = 'k'
bb =dict(facecolor='white', alpha=0.8)
for tx in range(len(paleolonstext)):
    if sitenames[tx] in ['U1404','925']:
        ax2.text(paleolonstext[tx]+4,paleolats[tx]-8, sitenames[tx], fontsize=15,
             horizontalalignment='left', color=tc, bbox=bb)    
    elif sitenames[tx] in ['SSQ']:   
        ax2.text(paleolonstext[tx]-15,paleolats[tx]-8, sitenames[tx], fontsize=15,
             horizontalalignment='left', color=tc, bbox=bb)    
    else:
        ax2.text(paleolonstext[tx]+4,paleolats[tx]+4, sitenames[tx], fontsize=15,
             horizontalalignment='left', color=tc, bbox=bb)

#%% final
fig.subplots_adjust(right=0.95)
cbar_ax = fig.add_axes([0.95, 0.1, 0.1, 0.8])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'vertical', fraction = 0.8,
                    aspect=18)
cbar.ax.xaxis.set_label_position('bottom')

#fig.subplots_adjust(bottom=0.17)
#cbar_ax = fig.add_axes([0.135, 0., 0.8, 0.07])
#cbar_ax.set_visible(False)
#cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'horizontal', fraction = 1.2,
#                    aspect=18)
#cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.set_xlabel('km', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)


plt.savefig('Bathymetry.png', dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()