#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 09:35:05 2019

@author: nooteboom
"""

from parcels.scripts import convert_npydir_to_netcdf as cp
from glob import glob

dirr = '/home/nooteb/petern/parcels/Ilja/output/'

#temp_dirs= [
#]

temp_dirs = glob(dirr + 'out-*')
print(temp_dirs)
for i in range(len(temp_dirs)):
    if(i%10==0):
        print(i/len(temp_dirs))
    cp.convert_npydir_to_netcdf(temp_dirs[i])
