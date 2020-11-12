The .nc files contains the particle trajectory data. The parameter settings can be read from the file name:
- dd is the dwelling depth in m
- sp is the sinking speed in m/day
- 2pic means 2 times pre industrial carbon configuration

If dd equals -1, it is the configuration where particles are backtracked until they reach the mixed layer depth. Whether a particle trajectory has reached the mixed layer depth is given by the variable 'MLr' (1 if the particle has reached the mixed laer depth and is tracked along the mixed layer depth).

The files contain fields with 4 dimensions:
- 'lon' the longitude (deg E)
- 'lat' the latitude (deg N)
- 'traj' the particles
- 'obs' the observations of particles as they are advected by time

(so each particle 'traj' is advected and 'obs' are the observations of this particle at different time steps as it is advected)

The files contain these variables:
- 'lon' the longitudes of the particle trajectories
- 'lat' the latitudes  of the particle trajectories
- 'z' the depths (m)  of the particle trajectories
- 'time' the time (s with respect to reference time)  of the particles 
- 'age' the age of the particles. Starts counting if the particles reach their dwelling depth
- 'temp' the temperature (degC)
- 'salin' the salinity (psu) 
