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
import cmocean.cm as cm
import matplotlib.colors as mcolors
#%%
def z_masked_overlap(axe, X, Y, Z, source_projection=None):
    """
    for data in projection axe.projection
    find and mask the overlaps (more 1/2 the axe.projection range)

    X, Y either the coordinates in axe.projection or longitudes latitudes
    Z the data
    operation one of 'pcorlor', 'pcolormesh', 'countour', 'countourf'

    if source_projection is a geodetic CRS data is in geodetic coordinates
    and should first be projected in axe.projection

    X, Y are 2D same dimension as Z for contour and contourf
    same dimension as Z or with an extra row and column for pcolor
    and pcolormesh

    return ptx, pty, Z
    """
    if not hasattr(axe, 'projection'):
        return Z
    if not isinstance(axe.projection, ccrs.Projection):
        return Z

    if len(X.shape) != 2 or len(Y.shape) != 2:
        return Z

    if (source_projection is not None and
            isinstance(source_projection, ccrs.Geodetic)):
        transformed_pts = axe.projection.transform_points(
            source_projection, X, Y)
        ptx, pty = transformed_pts[..., 0], transformed_pts[..., 1]
    else:
        ptx, pty = X, Y


    with np.errstate(invalid='ignore'):
        # diagonals have one less row and one less columns
        diagonal0_lengths = np.hypot(
            ptx[1:, 1:] - ptx[:-1, :-1],
            pty[1:, 1:] - pty[:-1, :-1]
        )
        diagonal1_lengths = np.hypot(
            ptx[1:, :-1] - ptx[:-1, 1:],
            pty[1:, :-1] - pty[:-1, 1:]
        )
        to_mask = (
            (diagonal0_lengths > (
                abs(axe.projection.x_limits[1]
                    - axe.projection.x_limits[0])) / 2) |
            np.isnan(diagonal0_lengths) |
            (diagonal1_lengths > (
                abs(axe.projection.x_limits[1]
                    - axe.projection.x_limits[0])) / 2) |
            np.isnan(diagonal1_lengths)
        )

        # TODO check if we need to do something about surrounding vertices

        # add one extra colum and row for contour and contourf
        if (to_mask.shape[0] == Z.shape[0] - 1 and
                to_mask.shape[1] == Z.shape[1] - 1):
            to_mask_extended = np.zeros(Z.shape, dtype=bool)
            to_mask_extended[:-1, :-1] = to_mask
            to_mask_extended[-1, :] = to_mask_extended[-2, :]
            to_mask_extended[:, -1] = to_mask_extended[:, -2]
            to_mask = to_mask_extended
        if np.any(to_mask):

            Z_mask = getattr(Z, 'mask', None)
            to_mask = to_mask if Z_mask is None else to_mask | Z_mask

            Z = np.ma.masked_where(to_mask, Z)

        return ptx, pty, Z


#%%

sns.set(context='paper',style="whitegrid",font="Arial")
fs = 17
font = {'size'   : fs}
matplotlib.rc('font', **font)

#variables
sp = 6
dd = 10
projection = ccrs.PlateCarree(180)
exte = [200, 250, -75, -50]
cmap = 'viridis'#cm.deep
paleoloc = [360-34.88303588,33.02928823]

rr = 1 # determines the resolution reduction for the plotting
#%% Load files
contours = [-6e-5,-5e-5,-4e-5,-3e-5,-2e-5,-1e-5,0,1e-5, 2E-5, 3E-5, 5E-5]

#NorthAtlantic : 
minlat = 0
maxlat = 75
minlat38 = 0
minlon=280
maxlon=359
vsbath = [-0.04,6]
colors1 = cm.deep(np.linspace(0, 1, 100))
colors2 = matplotlib.cm.get_cmap('Greys')(np.linspace(0, 1, 83))[60:61]
colors = np.vstack((colors2, colors1))
cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
#cmap = cm.deep

# Load 38MA
dirRead = '/Users/nooteboom/Documents/PhD/Eocene_POP/OUPUT/gridf/'
nc38 = Dataset(dirRead + 'kmt_tx0.1_POP_EO38_poleland.nc')
nc38grid = Dataset(dirRead + 'grid_coordinates_pop_tx0.1_38ma.nc')
lats38 = nc38grid['U_LAT_2D'][:]
lons38 = nc38grid['U_LON_2D'][:] + 25
lons38[lons38>180] -= 360
exte38 = [minlon, maxlon, minlat38, maxlat]
idx = np.where(np.logical_and(lats38[:,0]>=minlat38, lats38[:,0]<=maxlat))

bath38 = (-1*nc38['elevation'][:] / 1000)[idx[0]]
lats38 = lats38[idx[0]]
lons38 = lons38[idx[0]]
#%% start figure
latticks = [0,25, 50, 75]
fig = plt.figure(figsize=(14,14))
#%% subplot (a)
ax2 = plt.subplot(111, projection=projection)
plt.title('38Ma bathymetry', fontsize=fs)

g = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
g.xlabels_top = False
g.ylabels_right = False
g.xlabels_bottom = True
g.xlabels_left = True
g.xlabel_style = {'fontsize': fs}
g.ylabel_style = {'fontsize': fs}
g.xformatter = LONGITUDE_FORMATTER
g.yformatter = LATITUDE_FORMATTER
g.ylocator = mticker.FixedLocator(latticks)
g.xlocator = mticker.FixedLocator([240, 280,320, 360])
ax2.set_extent(exte38, ccrs.PlateCarree())

im2 = plt.pcolormesh(lons38[::rr,::rr], lats38[::rr,::rr], bath38[::rr,::rr], cmap=cmap, 
                     vmin=vsbath[0], vmax=vsbath[1],
                     transform=ccrs.PlateCarree()
                     )

plt.scatter([paleoloc[0]],[paleoloc[1]], marker='P', s=1000, color='red', 
            zorder=3000,transform=ccrs.PlateCarree())
#%% final
fig.subplots_adjust(right=0.95)
cbar_ax = fig.add_axes([0.95, 0.1, 0.1, 0.8])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'vertical', fraction = 0.8,
                    aspect=18)
cbar.ax.xaxis.set_label_position('bottom')

cbar.ax.set_xlabel('km', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)


plt.savefig('FH_bathymetry_NA.png', dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()