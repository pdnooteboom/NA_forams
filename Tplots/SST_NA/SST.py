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

def create_locidx(LATS,LONS,minlat,maxlat,minlon,maxlon):
    bolat = np.logical_and(LATS[:,0]>=minlat-5,LATS[:,0]<=maxlat+1)
    if(minlon<180 and maxlon>180):
        bolon = np.logical_and(LONS[50,:]<=minlon,
                               LONS[50,:]<=maxlon+5) 
    else:
        bolon = np.logical_and(LONS[50,:]>=minlon,LONS[50,:]<=maxlon+5)    
    return np.where(bolat), np.where(bolon)

def around_zerolon(lons, lats, field):
    idn = np.where(lons[50]<0)
    idn2 = np.where(lons[50]>=0)
    lats = np.concatenate((lats[:,idn2[0]], lats[:,idn[0]]), axis=1)
    field = np.concatenate((field[:,idn2[0]], field[:,idn[0]]), axis=1)
    return lons, lats, field       
#%%
var = 'SST'# 'T300' #

sns.set(context='paper',style="whitegrid",font="Arial")
fs = 30
font = {'size'   : fs}
matplotlib.rc('font', **font)

lonticks = [-100,-75,-50,-25,0,25,50,75]

#variables
gw = True

exte = [1, 360, -75, -10]
cmap = cm.thermal
#cmap='nipy_spectral'
if(var=='SST'):
    vsbath = [5,32]
    contours = [12.5,17.5,22.5,27.5,30]
elif(var=='T300'):
#    vsbath = [10,22.5]
#    contours = [12.5,15,17.5,20]
    vsbath = [5,32]
    contours = [12.5,17.5,22.5,27.5,30]
    
minlat = 0
maxlat = 75
minlat38 = 0
minlon=280-360
maxlon=359-360
paleoloc = [-34.88303588,33.02928823]

if((minlon<180 and maxlon>180)):
    projection = ccrs.PlateCarree(180)
else:
    projection = ccrs.PlateCarree(180)
    projection2 = ccrs.PlateCarree(180)

assert ~(minlon<-180 and maxlon>-180), 'pick minlon<180 and maxlon>180'   
#%% Load files
# Load 38MA
dirReadHR = '/Volumes/HD/Eocene/output/time_mean/'
nc38 = Dataset(dirReadHR + 'pop_38Ma_avg6years.nc')
dirRead = '/Users/nooteboom/Documents/PhD/Eocene_POP/OUPUT/gridf/'
nc38grid = Dataset(dirRead + 'grid_coordinates_pop_tx0.1_38ma.nc')
dz = nc38grid['w_dep'][1:] - nc38grid['w_dep'][:-1]
lats38 = nc38['ULAT'][:]
lons38 = nc38['ULONG'][:]+25
exte38 = [minlon, maxlon, minlat38, maxlat]
idxlat, idxlon = create_locidx(lats38,lons38,minlat,maxlat,minlon,maxlon)


if(var=='SST'):
    bath38 = (nc38['TEMP'][0])[np.min(idxlat[0]):np.max(idxlat[0]), np.min(idxlon[0]):np.max(idxlon[0])]
elif(var=='T300'): # temperature at 317.5 m depth
    bath38 = (nc38['TEMP'][16])[np.min(idxlat[0]):np.max(idxlat[0]), np.min(idxlon[0]):np.max(idxlon[0])]

lats38 = lats38[np.min(idxlat[0]):np.max(idxlat[0]), np.min(idxlon[0]):np.max(idxlon[0])]
lons38 = lons38[np.min(idxlat[0]):np.max(idxlat[0]), np.min(idxlon[0]):np.max(idxlon[0])]
#%% start figure
class nf(float):
    def __repr__(self):
        s = f'{self:.1f}'
        return f'{self:.0f}' if s[-1] == '0' else s
fig = plt.figure(figsize=(15,13))

#% subplot (a)
ax2 = plt.subplot(111, projection=projection)
if(var=='SST'):
    plt.title('SST 0.1$^{\circ}$ resolution', fontsize=fs)
elif(var=='T300'): # temperature at 317.5 m depth
    plt.title('317m T, 0.1$^{\circ}$ resolution', fontsize=fs)


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
g.xlocator = mticker.FixedLocator(lonticks)
ax2.set_extent(exte38, ccrs.PlateCarree())

im2 = plt.pcolormesh(lons38, lats38, bath38, cmap=cmap, 
                     vmin=vsbath[0], vmax=vsbath[1],
                     transform=ccrs.PlateCarree()
                     )

X, Y, masked_MDT = z_masked_overlap(ax2, lons38, lats38, 
                                    bath38,
                                    source_projection=ccrs.Geodetic())

CS = ax2.contour(X-180, Y, masked_MDT, contours, 
            colors='k', linestyles='-',
                     transform=ccrs.PlateCarree(),zorder=300
                     )
#CS.levels = [nf(val) for val in CS.levels]
#fmt = '%r '
#ax2.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)

land = np.full(bath38.shape, np.nan); land[bath38==0] = 1;
plt.pcolormesh(lons38, lats38, land,  
               vmin=0, vmax=1.6,cmap='binary',
                     transform=ccrs.PlateCarree(),zorder=2
                     )

plt.scatter([paleoloc[0]],[paleoloc[1]], marker='P', s=1000, color='k', 
            zorder=3000,transform=ccrs.PlateCarree())

#% final
fig.subplots_adjust(right=0.95)
cbar_ax = fig.add_axes([0.95, 0.1, 0.1, 0.8])
cbar_ax.set_visible(False)
cbar = fig.colorbar(im2, ax=cbar_ax, orientation = 'vertical', fraction = 0.8,
                    aspect=18)
cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.set_xlabel('$^{\circ}$C', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)
cbar.add_lines(CS)

if(var=='SST'):
    plt.savefig('meanSST_NA.png', dpi=300,bbox_inches='tight',pad_inches=0)
elif(var=='T300'): # temperature at 317.5 m depth
    plt.savefig('meanT300_NA.png', dpi=300,bbox_inches='tight',pad_inches=0)

plt.show()