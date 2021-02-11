#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:11:18 2021

@author: nooteboom
"""

import numpy as np
import matplotlib.pylab as plt
from netCDF4 import Dataset
from math import sin, cos, sqrt, atan2, radians
from numba import jit
import pandas as pd

@jit(nopython=True)
def distance(lon1, lat1, lon2, lat2):
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    # approximate radius of earth in km
    R = 6373.0
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    
    return R * c  # in km
    
@jit(nopython=True)
def calculate_distance(loc, lons, lats, bath):
    assert len(loc)==2, 'loc is a list which should contain the longitude and latitude of one location'
    assert len(lons.shape)==len(lats.shape)
    assert len(lats.shape)==1

    sd = 10**10
    sd2 = 10**10
    bath0 = np.nan
    for i in range(len(lons)):
        if(bath[i]<0):
            dis = distance(loc[0], loc[1], lons[i], lats[i])
            if(dis<sd):
                sd = dis
        else:
            dis = distance(loc[0], loc[1], lons[i], lats[i])
            if(dis<sd2):
                sd2 = dis
                bath0 = bath[i]   
    return sd, bath0 # in km

#%% load sediment sample site locations
ff = np.load('paleolocs_38Ma.npz')
paleolons = ff['longitudes'] 
paleolons[paleolons<180] += 180
paleolats = ff['latitudes']
sitenames = ff['names']
#%%
minlat = -80
maxlat = 80
dirRead = '/Users/nooteboom/Documents/PhD/Eocene_POP/OUPUT/gridf/'

# Load 38MA bathyemtry
nc38 = Dataset(dirRead + 'kmt_tx0.1_POP_EO38_poleland.nc')
nc38grid = Dataset(dirRead + 'grid_coordinates_pop_tx0.1_38ma.nc')
lats38 = nc38grid['U_LAT_2D'][:]
lons38 = nc38grid['U_LON_2D'][:] + 25
minlon = np.min(lons38)
maxlon = np.max(lons38)
exte38 = [minlon, maxlon, minlat, maxlat]
idx = np.where(np.logical_and(lats38[:,0]>=minlat, lats38[:,0]<=maxlat))

bath38 = (-1*nc38['elevation'][:] / 1000)[idx[0]]
fH38 = (-2 * 7.2921 * 10**(-5) * np.sin(lats38[idx[0]]*np.pi/180) / bath38)
lats38 = lats38[idx[0]]
lons38 = lons38[idx[0]]
#%%
shortestdistances= []
baths= []
#idx = np.where(bath38<0)
for j in range(len(paleolons)):
    lo = paleolons[j]
    la = paleolats[j]
    dis, ba = calculate_distance([lo, la], 
                                lons38.flatten(),
                                lats38.flatten(), bath38.data.flatten())
    shortestdistances.append(dis)
    baths.append(ba)
#%% write as csv
df = pd.DataFrame(np.swapaxes(np.array([paleolons.tolist(), paleolats.tolist(), 
                            sitenames.tolist(), shortestdistances, baths]),0,1),
                   columns=['paleolongitudes', 'paleolatitudes', 'names',
                            'shortest distance to land (km)', 'paleobathymetry (km)'])
df.to_csv('shortest_distances.csv')
    