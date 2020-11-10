#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:13:33 2020

@author: nooteboom
"""
from copy import copy
import numpy as np
import shapefile

def check_len(dirr):
    for di, d in enumerate(dirr.keys()):
        if(di==0):
            d0 = copy(d)
        assert len(dirr[d]) == len(dirr[d0]), d+' is not of the same size as plateid'

if(__name__=='__main__'):
    # from Bijl et al. 2011
    # 36-28 Ma
    sites = { 
        # Put here the plateIDs, which will be use to rotate the points eventually
        'plateid': np.array([101]),
            # The names of the sites:
            'names': np.array(['U1411']),
            # the longitudes (deg E) in the present-day:
            'plon_PD': np.array([-49.]),
            # the latitudes (deg N) in the present-day
            'plat_PD': np.array([41.618333]),
            }
    check_len(sites)
    
    #%%
    DirW = '' + 'sitelocs' # Where the file is written + name of file
    trees_shp = shapefile.Writer(DirW)
    trees_shp.autoBalance = 1
    #The fields in the shapefile:
    trees_shp.field("PLATEID1", "C")
    trees_shp.field("NAME", "C")
    trees_shp.field("GPGIM_TYPE", "C")
    trees_shp.field("FROMAGE", "C")
    trees_shp.field("TOAGE", "C")    
    

    counter = 1 # counter keeps track of the amount of sediment sites
    for i in range(len(sites['plon_PD'])): 
        #Write the properties of a single sediment site:
        plateid1 = sites['plateid'][i]
        name = sites['names'][i]
        gpgim_type = 'gpmi:UnclassifiedFeature'
        toage = '-999'
        fromage = '600'        
        
        trees_shp.point(float(sites['plon_PD'][i]),float(sites['plat_PD'][i]))
        
        trees_shp.record(plateid1, name, gpgim_type, fromage, toage)

        print("Feature " + str(counter) + " added to Shapefile.")
        counter += 1
        
    trees_shp.close()
