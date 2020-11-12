import numpy as np
from netCDF4 import Dataset

dd = 50
sp = 200
config = '2pic'

nop = 146

lonCOR = 25
ncr = Dataset('Forams_dd%d_sp%d_%s.nc'%(dd, sp, config))
obs = len(ncr['lat'][0])
lats = ncr['lat'][:,0]
lons = ncr['lon'][:,0] + lonCOR
ulon = np.unique(lons)
ulat = np.unique(lats)
ulons, ulats = np.meshgrid(ulon, ulat)

if(dd==-1):
    resk = ['lon', 'lat', 'z', 'temp', 'salin', 'age', 'time', 'MLr']
else:
    resk = ['lon', 'lat', 'z', 'temp', 'salin', 'age', 'time']
res = np.full((len(resk), ulons.shape[0], ulons.shape[1], nop, obs), np.nan)

for ri in range(len(resk)):
    if(True):
        for lo in range(ulons.shape[0]):
            for la in range(ulats.shape[1]): 
                 idx = np.where(np.logical_and(lons==ulon[lo], lats==ulat[la]))[0]
                 if(resk[ri]=='lon'):
                     ch = ncr[resk[ri]][idx] + lonCOR
                 else:
                     ch = ncr[resk[ri]][idx]
                 mas = ch.mask
                 ch = ch.data
                 ch[mas] = np.nan
                 for j in range(ch.shape[0]):
                     res[ri,lo, la, j] = ch[j]

# Write the netcdf file
ncw = Dataset('GridForam_dd%d_sp%d_%s.nc'%(dd,sp,config), 'w')
lat = ncw.createDimension('latd', len(ulat))
lon = ncw.createDimension('lond', len(ulon))
traj = ncw.createDimension('traj', nop)
obsd = ncw.createDimension('obs', obs)

vlon = ncw.createVariable('lon', np.float32, ('lond','latd','traj','obs',))
vlon[:] = res[0]
vlat = ncw.createVariable('lat', np.float32, ('lond','latd','traj','obs',))
vlat[:] = res[1]
vz = ncw.createVariable('z', np.float32, ('lond','latd','traj','obs',))
vz[:] = res[2]
vtemp = ncw.createVariable('temp', np.float32, ('lond','latd','traj','obs',))
vtemp[:] = res[3]
vsalt = ncw.createVariable('salin', np.float32, ('lond','latd','traj','obs',))
vsalt[:] = res[4]
vage = ncw.createVariable('age', np.float32, ('lond','latd','traj','obs',))
vage[:] = res[5]
vtime = ncw.createVariable('time', np.float32, ('lond','latd','traj','obs',))
vtime[:] = res[6]
if(dd==-1):
    vML = ncw.createVariable('MLr', np.float32, ('lond','latd','traj','obs',))
    vtime[:] = res[7]

ncw.close()




