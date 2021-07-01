import xarray as xr
import numpy as np  # needed for functions!

# source
exec(open("./src/penman_monteith_FAO56_FT.py").read())

# main variables for 1981-01-01
Z = xr.open_dataset("data/example_FT/DEM.nc").DEM
U2 = xr.open_dataset("data/example_FT/ws_01.nc").layer
LAT = xr.open_dataset("data/example_FT/Y.nc").Y
Tmin = xr.open_dataset("data/example_FT/tmin_1981-01-01.nc").layer
Tmax = xr.open_dataset("data/example_FT/tmax_1981-01-01.nc").layer
Hs = xr.open_dataset("data/example_FT/hs_1981-01-01.nc").layer
Td = xr.open_dataset("data/example_FT/td_1981-01-01.nc").layer


pm_evap = penman_monteith_FAO56_FT(time_i=1,
                                   tmax_i=Tmax.values,
                                   tmin_i=Tmin.values,
                                   sd_i=Hs.values,
                                   td_i=Td.values,
                                   lat_i=LAT.values,
                                   z_i=Z.values,
                                   u2_i=U2.values)





exp = xr.open_dataset("data/example_FT/DEM.nc")
exp["pm"] = (('latitude', 'longitude'), pm_evap)
exp["pm"].plot()

exp["pm"].to_netcdf("/home/adrian/Documents/Repos/Evapotranspiration/data/example_FT/output.nc")


### parallel

import urllib.request
import pandas as pd
import numpy as np
import datetime
import xarray as xr
from joblib import Parallel, delayed

PM_function = urllib.request.urlopen("https://raw.githubusercontent.com/adrHuerta/Evaporation/master/src/penman_monteith_FAO56_FT.py")
exec(PM_function.read())

range_time = pd.date_range("1981-01-01", "1981-01-10", freq="d")

def apply_PM_function(time_step):

    time_i = int(time_step.strftime('%j'))
    Tmax = xr.open_dataset("data/example_FT/example/tmax/tmax" + "_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
    Tmin = xr.open_dataset("data/example_FT/example/tmin/tmin" + "_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
    Hs = xr.open_dataset("data/example_FT/example/hs/hs" + "_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
    Td = xr.open_dataset("data/example_FT/example/td/td" + "_" + time_step.strftime('%Y-%m-%d') + ".nc").layer
    U2 = xr.open_dataset("data/example_FT/example/ws/ws" + "_" + time_step.strftime('%m') + ".nc").layer
    LAT = xr.open_dataset("data/example_FT/Y.nc").Y
    Z = xr.open_dataset("data/example_FT/DEM.nc").DEM

    pm_evap = penman_monteith_FAO56_FT(time_i=time_i,
                                       tmax_i=Tmax.values,
                                       tmin_i=Tmin.values,
                                       sd_i=Hs.values,
                                       td_i=Td.values,
                                       lat_i=LAT.values,
                                       z_i=Z.values,
                                       u2_i=U2.values)

    to_save = xr.open_dataset("data/example_FT/DEM.nc")
    to_save["pm"] = (('latitude', 'longitude'), pm_evap)
    to_save["pm"].to_netcdf("data/example_FT/example/eo/eo" + "_" + time_step.strftime('%Y-%m-%d') + ".nc")

Parallel(n_jobs=2, verbose=50)(
    delayed(apply_PM_function)(i) for i in range_time
)

xr.open_dataset("/home/adrian/Documents/Repos/Evapotranspiration/data/example_FT/example/eo/eo_1981-01-01.nc").pm.plot()