"""
This script derives the condensation rate from 3D model output for a variety of idealized MCS cases. 

kukulies@ucar.edu

"""

import numpy as np
import datetime
import pandas as pd 
import xarray as xr
import wrf
from netCDF4 import Dataset
from pathlib import Path
import microphysics_functions as micro
import warnings
warnings.filterwarnings("ignore")

##### 10 MCS cases #####
path = Path('/glade/scratch/kukulies/idealized_storms/23_2007-06-19_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5/4000/')
caseIDs = ['23']
resolution = '4000' 

for caseID in caseIDs:
    # get times
    year = '2011'
    month = '07'
    day  = '13'
    date = year + month + day
    start = datetime.datetime(int(year), int(month), int(day), 0, 0)
    end = datetime.datetime(int(year), int(month), int(day), 7, 0)
    times = pd.date_range(start, end, freq= '1MIN') 
    print('deriving data for case '+ caseID + ' ' +date+ ' with' + resolution +  ' meter resolution.')
    files = list(path.glob('wrfout*'))
    files.sort()
    assert times.shape[0]  == len(files)
    for fname in files:
        print(fname)
        #### Read in data for timestep
        mcs_case = xr.open_dataset(fname).squeeze()
        wrfin = Dataset(fname)
        # get variables 
        iwc = mcs_case.QSNOW +  mcs_case.QICE +  mcs_case.QGRAUP
        lwc =mcs_case.QRAIN + mcs_case.QCLOUD
        precip = mcs_case.RAINNC
        # get vertical velocity on mass points instead of staggered
        vertical_velocity = wrf.getvar(wrfin, 'wa')
        temp = wrf.getvar(wrfin, 'tk')
        qcloud = mcs_case.QCLOUD
        pressure = wrf.getvar(wrfin, 'pres')
        base_pressure =  mcs_case.PB


        # integate iwc and lwp over pressure
        iwp = micro.pressure_integration(iwc.data, -pressure.data)
        lwp = micro.pressure_integration(lwc.data, -pressure.data)

        #### Derive condensation rate (in kg/kg/s)
        condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temp, pressure, base_pressure)
        rho = micro.get_air_density(pressure, temp).squeeze()
        condensation_rate_s = condensation_rate_s 

        ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
        condensation_cloud = condensation_rate_s.where(qcloud > 0, 0 )
        # we are only interested in positive values, as negative values are evaporation 
        condensation_masked = condensation_cloud.where(condensation_cloud > 0, 0 ).data
        ### integrate over pressure levels to get kg/m2/s
        condensation = micro.pressure_integration(condensation_masked,-pressure.data)

        #### Concatenate
        if fname == files[0]:
            tiwp = iwp
            tlwp = lwp
            condensation_rate = condensation
            surface_precip = precip
        else:
            surface_precip = np.dstack((surface_precip, precip)) 
            tiwp  = np.dstack((tiwp, iwp ))
            tlwp = np.dstack((tlwp, lwp))
            condensation_rate = np.dstack((condensation_rate, condensation))

    #### Save data for the case to netCDF4
    data_vars = dict(tiwp=(["south_north", "west_east", "time"], tiwp),
                     tlwp=(["south_north", "west_east", "time"], tlwp),
                     surface_precip=(["south_north", "west_east", "time"],surface_precip),
                     condensation_rate=(["south_north", "west_east", "time"], condensation_rate),
                     lats=(["south_north", "west_east"], mcs_case.XLAT.values),
                     lons=(["south_north", "west_east"], mcs_case.XLONG.values),)

    coords = dict(south_north=mcs_case.south_north.values, west_east=mcs_case.west_east.values, time = times)
    data = xr.Dataset(data_vars=data_vars, coords=coords)                                                                 
    data.to_netcdf('/glade/scratch/kukulies/idealized_storms/23_2007-06-19_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5/4000/Idealized_MCS_' + caseID + '_' + date + '_' + resolution + '.nc')       
