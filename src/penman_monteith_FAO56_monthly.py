"""
Penman Monteith FAO56 equation for estimation of
reference evapotranspiration based on Zotarelli et al. (2010) and http://www.fao.org/3/x0490e/x0490e08.htm
At ten-day/monthly time step the computation is quite similar to the daily time step, but some changes are
important (see comments)

Parameters:
time_i : number of the day (julian day) since 1 january (centroid value, eg if is january is 15, and so on)
tmax_i : ten-day/monthly maximum temperature (°C)
tmin_i : ten-day/monthly minimum temperature (°C)
rs_i : monthly solar radiation  MJ/m(-2).day
lat_i : latitude in decimal degree (°)
z_i : elevation above sea level (m)
u2_i : mean wind speed (m/s)
tmax_i_minus1 : ten-day/monthly maximum temperature of the previous month (°C)
tmin_i_minus1 : ten-day/monthly minimum temperature of the previous month (°C)

Comments:
Some simplifications are done:
i) wind speed is automatically 2 m/s
ii) actual vapor pressure is only computed from tmin_i
iii) as at ten-day/monthly time step the soil heat flux (G) can
no be simplified, G is computed from the mean monthly air temperatures of
the previous and next month (Tmean_i - Tmean_i_minus1)

References:
- Zotarelli, Lincoln, et al. "Step by step calculation of the Penman-Monteith
Evapotranspiration (FAO-56 Method)." Institute of Food and Agricultural Sciences.
University of Florida (2010).
- http://www.fao.org/3/x0490e/x0490e08.htm
"""

def penman_monteith_FAO56_monthly(time_i, tmax_i, tmin_i, rs_i, lat_i, z_i, u2_i, tmax_i_minus1, tmin_i_minus1):
    # step 1
    tmean_i = (tmax_i + tmin_i)/2
    tmean_i_minus1 = (tmax_i_minus1 + tmin_i_minus1)/2
    # step 2
    rs_i = rs_i
    # step 3
    u2_i = u2_i
    # step 4
    DELTA_SlopeSat_i = 4098 * (0.6108 * np.exp( (17.27 * tmean_i) / (tmean_i + 237.3) )) / (np.power(tmean_i + 237.3, 2))
    # step 5
    P_i = 101.3 * np.power((293 - 0.0065 * z_i) / 293, 5.26)
    # step 6
    GAMMA_Psychrometric_i = 0.000665 * P_i
    # step 7
    DT_i =  DELTA_SlopeSat_i / (DELTA_SlopeSat_i + GAMMA_Psychrometric_i * (1 + 0.34 * u2_i))
    # step 8
    PT_i = GAMMA_Psychrometric_i / (DELTA_SlopeSat_i + GAMMA_Psychrometric_i * (1 + 0.34 * u2_i))
    # step 9
    TT_i = (900 / (tmean_i + 273)) * u2_i
    # step 10
    e_tmax_i = 0.6108 * np.exp((17.27 * tmax_i) / (tmax_i + 237.3))
    e_tmin_i = 0.6108 * np.exp((17.27 * tmin_i) / (tmin_i + 237.3))
    es_i = (e_tmax_i + e_tmin_i) / 2
    # step 11
    ea_i = e_tmin_i
    # step 12
    dr_i = 1 + 0.033 * np.cos(2 * np.pi * time_i / 365)
    delta_SD_i = 0.409 * np.sin((2 * np.pi * time_i / 365) - 1.39)
    # step 13
    lat_rad_i = np.pi * lat_i / 180
    # step 14
    Ws_i = np.arccos(-np.tan(lat_rad_i)*np.tan(delta_SD_i))
    # step 15
    Ra_i = (24 * 60 / np.pi) * 0.0820 * dr_i * ((Ws_i * np.sin(lat_rad_i) * np.sin(delta_SD_i)) + (np.cos(lat_rad_i) * np.cos(delta_SD_i) * np.sin(Ws_i)))
    # step 16
    Rso_i = (0.75 + 2 * np.float_power(10, -5) * z_i) * Ra_i
    # step 17
    Rns_i = (1 - 0.23) * rs_i
    # step 18
    Rnl_i = 4.903 * np.float_power(10, -9) * ((np.power(tmax_i + 273.16, 4) + np.power(tmin_i + 273.16, 4))/2) * (0.34 - 0.14 * np.sqrt(ea_i)) * ((1.35 * rs_i/Rso_i) - 0.35)
    # step 19
    Rn_i = Rns_i - Rnl_i
    G_i = 0.14 * (tmean_i - tmean_i_minus1)
    Rng_i = 0.408 * (Rn_i - G_i)
    # step 20
    ET_rad_i = DT_i * Rng_i
    # step 21
    ET_wind_i = PT_i * TT_i * (es_i - ea_i)
    # step 22
    response = ET_rad_i + ET_wind_i

    return np.round(response, 1)