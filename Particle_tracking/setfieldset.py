import numpy as np
from parcels import (FieldSet, Field)

def set_the_ML_fieldset(ufiles, bfile):
    filenames = { 'U': {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'V' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'W' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'S' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'T' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles} ,
                }

    variables = {'U': 'UVEL',
                 'V': 'VVEL',
                 'W': 'WVEL',
                 'T': 'TEMP',
                 'S': 'SALT'}
                 #'B':'Bathymetry'}

    dimensions = {
                  'U':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep','time':'time'},#
                  'V':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep','time':'time'},#
                  'W':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep','time':'time'},#
                  'T':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time':'time'},#
                  'S':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time':'time'}}

    indices = {'lat':range(1300,2280), 'lon':range(0,1000)}

    timestamps = np.expand_dims(np.array([np.datetime64('2009-12-31') - np.timedelta64(x,'D') for x in range(len(ufiles))])[::-1], axis=1)
    fieldset = FieldSet.from_pop(filenames, variables, dimensions, allow_time_extrapolation=False, field_chunksize=False, indices=indices,
                                     timestamps=timestamps)

    # add the bathymetry field (independent of time and depth)
    bfiles = {'lon': bfile, 'lat': bfile, 'data': [bfile, ]}
    bvariables = ('B', 'bathymetry')
    bdimensions = {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D'}
    Bfield = Field.from_netcdf(bfiles, bvariables, bdimensions, allow_time_extrapolation=True, interp_method='bgrid_tracer', field_chunksize=False)
    fieldset.add_field(Bfield, 'B')

    # add the mixed layer depth field (independent of depth)
    MLfiles = {'lon': bfile, 'lat': bfile, 'data': ufiles}
    MLvariables = ('ML', 'HMXL')
    MLdimensions = {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D'}
    MLfield = Field.from_netcdf(MLfiles, MLvariables, MLdimensions, interp_method='bgrid_tracer', field_chunksize=False, indices=indices,
                                     timestamps=timestamps)
    fieldset.add_field(MLfield, 'ML')

    fieldset.U.vmax = 10
    fieldset.V.vmax = 10
    fieldset.W.vmax = 10
    return fieldset


def set_the_fieldset(ufiles, bfile):
    filenames = { 'U': {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'V' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'W' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'S' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles},
                'T' : {'lon': bfile,
                        'lat': bfile,
                        'depth': bfile,
                        'data':ufiles} ,
                }

    variables = {'U': 'UVEL',
                 'V': 'VVEL',
                 'W': 'WVEL',
                 'T': 'TEMP',
                 'S': 'SALT'}
                 #'B':'Bathymetry'}

    dimensions = {
                  'U':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep','time':'time'},#
                  'V':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep','time':'time'},#
                  'W':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep','time':'time'},#
                  'T':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time':'time'},#
                  'S':{'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time':'time'}}

    indices = {'lat':range(1300,2280), 'lon':range(0,1000)}

    bfiles = {'lon': bfile, 'lat': bfile, 'data': [bfile, ]}
    bvariables = ('B', 'bathymetry')
    bdimensions = {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D'}

    timestamps = np.expand_dims(np.array([np.datetime64('2009-12-31') - np.timedelta64(x,'D') for x in range(len(ufiles))])[::-1], axis=1)
    fieldset = FieldSet.from_pop(filenames, variables, dimensions, allow_time_extrapolation=False, field_chunksize=False, indices=indices,
                                     timestamps=timestamps)
    Bfield = Field.from_netcdf(bfiles, bvariables, bdimensions, allow_time_extrapolation=True, interp_method='bgrid_tracer', field_chunksize=False)
    fieldset.add_field(Bfield, 'B')
    fieldset.U.vmax = 10
    fieldset.V.vmax = 10
    fieldset.W.vmax = 10
    return fieldset

