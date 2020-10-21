"""
Hargreaves-Samani equation for estimation of
reference evapotranspiration (mm/day) based on Vanderlinden et al. (2004) and
Shuttleworth (1993).

Parameters:
time_i : number of the day (julian day) since 1 january
tmax_i : daily maximum temperature (°C)
tmin_i : daily minimum temperature (°C)
lat_i : latitude in decimal degree (°)

References:
- Shuttleworth, W. J. 1993. ‘‘Evaporation.’’ Handbook of hydrology, D.
R. Maidment, ed., McGraw-Hill, New York
- Vanderlinden K, Giráldez JV, Van Meirvenne M (2004)
Assessing reference evapotranspiration by the Hargreaves method in southern
Spain. J Irrig Drain Eng 130(3):184–191
"""

def hargreaves_samani(time_i, tmax_i, tmin_i, lat_i):

    # equation 4.4.3
    delta = 0.4093 * np.sin((2 * np.pi * time_i / 365) - 1.405)
    # lat to radians
    lat_i = np.radians(lat_i)
    # equation 4.4.2
    W_s = np.arccos(-np.tan(lat_i) * np.tan(delta))
    # equation 4.4.5
    d_r = 1 + 0.033 * np.cos(2 * np.pi * time_i / 365)
    # equation 4.4.4
    Re = 15.392 * d_r * (
                W_s * np.sin(lat_i) * np.sin(delta) +
                np.cos(lat_i) * np.cos(delta) * np.sin(W_s)
    )
    # equation 4.2.44
    response = 0.0023 * Re * np.sqrt((tmax_i - tmin_i)) * ((tmax_i + tmin_i) / 2 + 17.8)

    return np.round(response, 1)
