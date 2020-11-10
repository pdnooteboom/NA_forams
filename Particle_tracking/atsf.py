# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:31:22 2017

@author: nooteboom
"""

from parcels import (FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable, Field)
from datetime import timedelta as delta
from datetime import  datetime
import numpy as np
import math
from glob import glob
import sys
import kernels as ke
import setfieldset as sf
import numpy as np

import warnings
import xarray as xr
#warnings.simplefilter("ignore", category=xr.SerializationWarning)

dirread_grid = '/projects/0/palaeo-parcels/Eocene/grids/'
dirread_U = '/projects/0/acc/pop/parcels/tx0.1/p21a.EO38Ma.tx0.1.2pic_control/tavg/'
dirread_T = '/projects/0/acc/pop/parcels/tx0.1/p21a.EO38Ma.tx0.1.2pic_control/movie/'

config = '2pic'
sp = 200. #The sinkspeed m/day
dd = float(sys.argv[1]) #The dwelling depth (m). if dd==-1, then track particle along mixed layer depth.
ls = 30 # days life span

dirwrite = 'output/'

# The 38Ma paleolocation of U1411
lonCOR = 25 # correction for longitude
latlon = [33.02928823, -34.88303588]
gs = 7 # Grid size
lonsz, latsz = np.mgrid[-gs:gs, -gs:gs]
lonsz = lonsz.flatten() + latlon[1]; latsz = latsz.flatten() + latlon[0];
lonsz -= lonCOR

print('amount of releas locations: ', len(latsz))
assert ~(np.isnan(latsz)).any(), 'locations should not contain any NaN values'
dep = dd * np.ones(latsz.shape)

times = np.array([datetime(2009, 12, 30) - delta(days=x) for x in range(0,2*int(365),5)])
time = np.empty(shape=(0));lons = np.empty(shape=(0));lats = np.empty(shape=(0));
for i in range(len(times)):
    lons = np.append(lons,lonsz)
    lats = np.append(lats, latsz)
    time = np.append(time, np.full(len(lonsz),times[i])) 
#%%
def run_particles(dirwrite,outfile,lonss,latss,dep):
    ufiles = sorted(glob(dirread_U + 't.p21a.EO38Ma.tx0.1.2pic_control.0036????.nc'))#[:150]
    bfile = dirread_grid+'kmt_tx0.1_POP_EO38.nc'
    print('set the fields')
    if(dd==-1):
        fieldset = sf.set_the_ML_fieldset(ufiles, bfile)
    else:
        fieldset = sf.set_the_fieldset(ufiles, bfile)
    fieldset.add_constant('dwellingdepth', np.float(dd))
    fieldset.add_constant('sinkspeed', sp/86400.)
    fieldset.add_constant('maxage', ls*86400)
    fieldset.add_constant('surface', 2.5)
    print('field set')
    if(dd==-1):
        class DinoParticle(JITParticle):
            temp = Variable('temp', dtype=np.float32, initial=np.nan)
            age = Variable('age', dtype=np.float32, initial=np.nan)
            salin = Variable('salin', dtype=np.float32, initial=np.nan)
            MLr = Variable('MLr', dtype=np.float32, initial=0)
    else:
        class DinoParticle(JITParticle):
            temp = Variable('temp', dtype=np.float32, initial=np.nan)
            age = Variable('age', dtype=np.float32, initial=np.nan)
            salin = Variable('salin', dtype=np.float32, initial=np.nan)
        
    pset = ParticleSet.from_list(fieldset=fieldset, pclass=DinoParticle, lon=lonss.tolist(), lat=latss.tolist(), 
                       time = time)

    pfile = ParticleFile(dirwrite + outfile, pset, 
                         outputdt=delta(days=1))

    if(dd==-1):
        kernels = ke.initialML + pset.Kernel(AdvectionRK4_3D)+ ke.SinkML + ke.SampleTemp + ke.AgeML 
    else:
        kernels = ke.initial + pset.Kernel(AdvectionRK4_3D)+ ke.Sink + ke.SampleTemp + ke.Age 
    pset.execute(kernels, runtime=delta(days=365*5), dt=delta(minutes=-5), output_file=pfile, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: ke.DeleteParticle,ErrorCode.Delete: ke.DeleteParticle})

    print('Execution finished')

outfile = 'Forams_dd'+str(int(dd)) +'_sp'+str(int(sp)) + '_%s'%(config)
run_particles(dirwrite,outfile,lons,lats,dep)

