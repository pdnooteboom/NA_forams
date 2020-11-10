#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:38:15 2020

@author: nooteboom
"""
import numpy as np
import shapefile


if(__name__=='__main__'):
    Ma = 38 # Million years ago
    DirR = '' + 'sitelocs/reconstructed_%dMa.shp'%(Ma)# # Where the file is written + name of file
    sf = shapefile.Reader(DirR)
    
    features = sf.shapeRecords()
    lons = np.zeros(len(features))
    lats = np.zeros(len(features))
    names = np.zeros(len(features)).astype(str)
    for i in range(len(features)):
        feature = features[i]
        first = feature.shape.__geo_interface__
        lons[i] = first['coordinates'][0]    
        lats[i] = first['coordinates'][1]
        names[i] = feature.record[5]
        
    print(lats)
    print(lons)
    print(names)
    # Write the sample locations as an .npz file
    np.savez('paleoloc_%dMa.npz'%(Ma),
             latitudes=lats, longitudes=lons, names=names)