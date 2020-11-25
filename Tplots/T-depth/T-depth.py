#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 15:00:43 2020

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
from netCDF4 import Dataset

paleoloc = [-34.88303588, 33.02928823] # the paleolocation

if(False): #if True, then perform the interpolation. Takes a while.
    # Load 38MA
    dirReadHR = '/Volumes/HD/Eocene/output/time_mean/'
    nc38 = Dataset(dirReadHR + 'pop_38Ma_avg6years.nc')
    dirRead = '/Users/nooteboom/Documents/PhD/Eocene_POP/OUPUT/gridf/'
    nc38grid = Dataset(dirRead + 'grid_coordinates_pop_tx0.1_38ma.nc')
    depth_t = (nc38grid['w_dep'][1:] + nc38grid['w_dep'][:-1]) / 2
    lats38 = nc38['ULAT'][1300:2400,200:1000]
    lons38 = nc38['ULONG'][1300:2400,200:1000]+25
    
    # the interpolation method used:
    method = 'nearest' # 'linear'
    
    # interpolate the temperature to the paleolocation for every depth level 
    points = np.concatenate((lons38.flatten()[:,np.newaxis], lats38.flatten()[:,np.newaxis]), axis=1)
    Tdep = np.zeros(len(depth_t))
    for z in range(len(depth_t)): # loop all depth levels
        print(z/len(depth_t))
        Tdep[z]= griddata(points, nc38['TEMP'][z,1300:2400,200:1000].flatten(), 
                            ([paleoloc[0]], [paleoloc[1]]), method=method)[0]
    
    np.savez('T-dep', Tdep=Tdep, depth_t=depth_t)

#%% plot figure
fs = 25

import seaborn as sns
sns.set(style='darkgrid', context='paper')

ff = np.load('T-dep.npz')
Tdep = ff['Tdep']; depth_t = ff['depth_t'] / 1000;
Tdep[Tdep==0] = np.nan # do not plot land numbers

fig = plt.figure(figsize=(8,10))
plt.gca().invert_yaxis()
plt.plot(Tdep, depth_t, '-o', linewidth=3, markersize=10)
plt.xlabel('temperature ($^{\circ}$C)', fontsize=fs)
plt.ylabel('depth (km)', fontsize=fs)
plt.xticks(fontsize=fs-2)
plt.yticks(fontsize=fs-2)
plt.savefig('T-dep.png', dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()