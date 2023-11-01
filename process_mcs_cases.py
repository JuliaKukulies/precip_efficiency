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
path = Path('/glade/campaign/mmm/c3we/prein/Idealized_MCSs/wrfout_files/WRF/')
caseIDs = ['03', '07', '10', '13', '17', '18', '19' ,'23', '38', '46']
resolution = '4000'

for caseID in caseIDs:
    if caseID in ['03', '07', '23']:
        loc = 'Loc2'
    else:
        loc = 'Loc1'
    subdir = list(path.glob(str(caseID + '*CTRL*' + loc +'*Storm*TH5')))[0]
    # get times
    year = subdir.name[3:7]
    month = subdir.name[8:10]
    day = subdir.name[11:13]
    date = year + month + day
    start = datetime.datetime(int(year), int(month), int(day), 0, 0)
    end = datetime.datetime(int(year), int(month), int(day), 7, 0)
    times = pd.date_range(start, end, freq= '5MIN')
    print('deriving data for case '+ caseID + ' ' +date+ ' with' + resolution +  ' meter resolution.')
    pattern = subdir / resolution
    files = list(pattern.glob('wrfout_d01*'))
    files.sort()
    assert times.shape[0]  == len(files)
    for fname in files:
        print(subdir.name, fname)
        #### Read in data for timestep
        mcs_case = xr.open_dataset(fname).squeeze()
        wrfin = Dataset(fname)
        # get variables 
        iwc = mcs_case.QSNOW +  mcs_case.QICE +  mcs_case.QGRAUP
        lwc =mcs_case.QRAIN + mcs_case.QCLOUD
        precip = mcs_case.RAINNC
        # get vertical velocity on mass points instead of staggered
        vertical_velocity = wrf.getvar(wrfin, 'wa')
        temp = mcs_case.TK
        qcloud = mcs_case.QCLOUD
        pressure = mcs_case.P + mcs_case.PB
        # integate iwc and lwp over pressure
        iwp = micro.pressure_integration(iwc.data, pressure.data)
        lwp = micro.pressure_integration(lwc.data, pressure.data)
        
        #### Derive condensation rate (in kg/kg/s)
        condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temp, pressure)
        ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
        condensation_masked = condensation_rate_s.where(qcloud > 0, 0 )
        # we are only interested in positive values, as negative values are evaporation 
        condensation_masked = condensation_masked.where(condensation_masked >= 0, 0 ).data
        ### integrate over pressure levels to get kg/m2/s
        condensation = micro.pressure_integration(condensation_masked,pressure.data)
        
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
    data.to_netcdf('mcs_cases' / ('Idealized_storm_caseID_' + caseID + '_' + date + '_' + resolution + '.nc') )       
