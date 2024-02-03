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
parent_path = Path('/glade/derecho/scratch/kukulies/idealized_mcs/')

cases =['19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '03_2011-07-16_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5',
        '23_2007-06-19_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5',
        '10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '13_2003-08-30_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '17_2011-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '18_2010-06-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '38_2007-08-04_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '46_2009-06-14_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '07_2011-07-04_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5',  ]

for case in cases:
    fnames = dict()
    path = parent_path / case / '4000'
    
    # file lists with different temporal output 
    fnames['5min'] = list(path.glob('wrfout_pr*'))
    fnames['5min'].sort()
    fnames['10min'] = fnames['5min'][::2]
    fnames['20min'] = fnames['5min'][::4]
    fnames['30min'] = fnames['5min'][::6]
    fnames['40min'] = fnames['5min'][::8]
    fnames['50min'] = fnames['5min'][::10]
    fnames['1H'] = fnames['5min'][::12]
                                                                
    year = case[3:7]
    month = case[8:10]
    day = case[11:13]
                                                        
    # get time series for case 
    date = year + month + day
    start = datetime.datetime(int(year), int(month), int(day), 0, 0)
    end = datetime.datetime(int(year), int(month), int(day), 7, 0)
    print(case, year, month, day, flush = True)

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

            ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
            condensation_cloud = condensation_rate_s.where(qcloud > 0, 0 )
            # we are only interested in positive values, as negative values are evaporation 
            condensation_masked = condensation_cloud.where(condensation_cloud > 0, 0 ).data
            ### integrate over pressure levels to get kg/m2/s
            condensation = micro.pressure_integration(condensation_masked,-pressure.data)

            # get air density
            rho = micro.get_air_density(pressure.squeeze().data, temp.squeeze().data)
            # vertical integration of process rate
            prw_vcd = mcs_case.PRW_VCD.where(mcs_case.PRW_VCD > 0, 0 ).data
            prwvcd_integrated = micro.pressure_integration(prw_vcd , -pressure.data)

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
        data.to_netcdf('/glade/derecho/scratch/kukulies/idealized_mcs/'+ case +'/4000/idealized_mcs_condensation_timestep_dependence_' + key + '_prw_vcd.nc')


