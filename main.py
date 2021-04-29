import xarray as xr
import numpy as np  # needed for functions!

# source
exec(open("./src/hargreaves_samani.py").read())
exec(open("./src/penman_monteith_FAO56.py").read())

# main variables for 1981-01-01
Z = xr.open_dataset("data/example_00/Z.nc").Z
U2 = xr.open_dataset("data/example_00/U2.nc").U2
LAT = xr.open_dataset("data/example_00/LAT.nc").LAT
Tmin = xr.open_dataset("data/example_00/Tmin.nc").Tmin
Tmax = xr.open_dataset("data/example_00/Tmax.nc").Tmax
Rs = xr.open_dataset("data/example_00/Rs.nc").Rs

# hs
hs_evap = hargreaves_samani(time_i=1,
                            tmax_i=Tmax.values,
                            tmin_i=Tmin.values,
                            lat_i=LAT.values)
# pm
pm_evap = penman_monteith_FAO56(time_i=1,
                                tmax_i=Tmax.values,
                                tmin_i=Tmin.values,
                                rs_i=Rs.values,
                                lat_i=LAT.values,
                                z_i=Z.values,
                                u2_i=U2.values)

# are pm computed in R == python? Yes! (I do this as a benchmarking)
exp = xr.open_dataset("/home/adrian/Documents/Repos/PISCOpet/PISCO_PMFAO56/1981-01-01.nc")
exp["pm_python"] = (('latitude', 'longitude'), pm_evap)
np.round(exp["pm_python"]/exp["layer"], 2).plot()