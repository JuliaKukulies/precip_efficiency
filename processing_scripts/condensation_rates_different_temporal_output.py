"""
This script derives the estimated and simulated condensation rates as well as accumulated precip for different temporal output from idealized MCS simulations.

kukulies@ucar.edu

"""

import numpy as np
import datetime
import pandas as pd 
import xarray as xr
import wrf
from netCDF4 import Dataset
from pathlib import Path
from microphysics import microphysics_functions as micro
import warnings
warnings.filterwarnings("ignore")

##### MCS cases #####
path = Path('/glade/derecho/scratch/kukulies/idealized_mcs/23_2007-06-19_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5/4000/')
year = '2007'
month = '06'
day  = '19'
fnames = dict()

# file lists with different temporal output 
fnames['5min'] = list(path.glob('wrfout_pr*'))
fnames['5min'].sort()
fnames['10min'] = fnames['5min'][::2]
fnames['20min'] = fnames['5min'][::4]
fnames['30min'] = fnames['5min'][::6]
fnames['40min'] = fnames['5min'][::8]
fnames['50min'] = fnames['5min'][::10]
fnames['1H'] = fnames['5min'][::12]

# get time series for case 
date = year + month + day
start = datetime.datetime(int(year), int(month), int(day), 0, 0)
end = datetime.datetime(int(year), int(month), int(day), 7, 0)

for key in fnames.keys():
    times = pd.date_range(start, end , freq = key)
    files = fnames[key]
    assert times.size == len(files)
    print(len(files),' files detected for output timestep: ', key ,  flush = True)
    
    # loop through output files
    for fname in files:
        #### Read in data for that timestep
        mcs_case = xr.open_dataset(fname).squeeze()
        wrfin = Dataset(fname)
        
        vertical_velocity = wrf.getvar(wrfin, 'wa')
        temp = wrf.getvar(wrfin, 'tk')
        qcloud = mcs_case.QCLOUD
        pressure = wrf.getvar(wrfin, 'pres')

        #### Derive condensation rate from hourly variables instead of minute-output (in kg/kg/s)
        condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temp, pressure)
        rho = micro.get_air_density(pressure, temp).squeeze()
        condensation_rate_s = condensation_rate_s
        
        ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
        condensation_cloud = condensation_rate_s.where(qcloud > 0, 0 )
        # we are only interested in positive values, as negative values are evaporation 
        condensation_masked = condensation_cloud.where(condensation_cloud > 0, 0 ).data
        ### integrate over pressure levels to get kg/m2/s
        condensation = micro.pressure_integration(condensation_masked,-pressure.data)
        
        # vertical integration of process rate
        prwvcd_integrated = micro.pressure_integration(mcs_case.PRW_VCD.squeeze().data, -pressure.data)
        
        # accumulate precip, take instaneous values for w,p,T 
        if fname == files[0]: 
            surface_precip = mcs_case.RAINNC
            prwvcd = prwvcd_integrated
            condensation_rate = condensation
        else:
            surface_precip = np.dstack((surface_precip,  mcs_case.RAINNC  )) 
            prwvcd = np.dstack((prwvcd,  prwvcd_integrated ))
            condensation_rate = np.dstack((condensation_rate, condensation ))
        
    #### Save data for the case to netCDF4
    data_vars = dict(condensation_rate=(["south_north", "west_east", "time"], condensation_rate),
                     surface_precip=(["south_north", "west_east", "time"], surface_precip),
                     prwvcd=(["south_north", "west_east", "time"], prwvcd),
                     lats=(["south_north", "west_east"], mcs_case.XLAT.values),
                     lons=(["south_north", "west_east"], mcs_case.XLONG.values),)

    coords = dict(south_north=mcs_case.south_north.values, west_east=mcs_case.west_east.values, time = times) 
    data = xr.Dataset(data_vars=data_vars, coords=coords)                                                                 
    data.to_netcdf('/glade/derecho/scratch/kukulies/idealized_mcs/23_2007-06-19_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5/4000/idealized_mcs_condensation_timestep_dependence_' + key + '.nc')


