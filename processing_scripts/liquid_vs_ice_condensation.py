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
path = Path('/glade/scratch/kukulies/idealized_storms/10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/4000/')
caseIDs = ['10']
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
    
        # get vertical velocity on mass points instead of staggered
        vertical_velocity = wrf.getvar(wrfin, 'wa')
        temp = wrf.getvar(wrfin, 'tk')
        qcloud = mcs_case.QCLOUD
        pressure = wrf.getvar(wrfin, 'pres')
        base_pressure =  mcs_case.PB

        # simulated process rates
        prw_vcd = mcs_case.PRW_VCD.where(mcs_case.PRW_VCD> 0, 0).squeeze()
        prw_vcd_liquid = prw_vcd.where(iwc <= 0, 0)
        prw_vcd_liquid_integration = micro.pressure_integration(prw_vcd_liquid.data, -pressure.data)
        prw_vcd_ice = prw_vcd.where(iwc > 0, 0)
        prw_vcd_ice_integration = micro.pressure_integration(prw_vcd_ice.data, -pressure.data)
 
        #### Derive condensation rate (in kg/kg/s)
        condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temp, pressure, base_pressure)
        ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
        condensation_cloud = condensation_rate_s.where(qcloud > 0, 0 )
        
        condensation_liquid = condensation_cloud.where(iwc <=0,0)
        condensation_ice= condensation_cloud.where(iwc > 0,0)
        # we are only interested in positive values, as negative values are evaporation 
        condensation_liquid = condensation_liquid.where(condensation_liquid > 0, 0 ).data
        condensation_ice = condensation_ice.where(condensation_ice > 0, 0 ).data
     
        ### integrate over pressure levels to get kg/m2/s
        condensation_liquid_int = micro.pressure_integration(condensation_liquid,-pressure.data)
        condensation_ice_int = micro.pressure_integration(condensation_ice,-pressure.data)

        #### Concatenate
        if fname == files[0]:
            condensation_rate_liquid = condensation_liquid_int
            condensation_rate_ice = condensation_ice_int
            prwvcd_ice  = prw_vcd_ice_integration
            prwvcd_liquid = prw_vcd_liquid_integration
        else:
            condensation_rate_liquid = np.dstack((condensation_rate_liquid, condensation_liquid_int))
            condensation_rate_ice = np.dstack((condensation_rate_ice, condensation_ice_int))
            prwvcd_liquid = np.dstack((prwvcd_liquid, prw_vcd_liquid_integration))
            prwvcd_ice = np.dstack((prwvcd_ice, prw_vcd_ice_integration))

    #### Save data for the case to netCDF4
    data_vars = dict(condensation_rate_liquid=(["south_north", "west_east", "time"], condensation_rate_liquid),
                     condensation_rate_ice=(["south_north", "west_east", "time"], condensation_rate_ice),
                     prwvcd_liquid=(["south_north", "west_east", "time"], prwvcd_liquid),
                     prwvcd_ice=(["south_north", "west_east", "time"], prwvcd_ice),
                     lats=(["south_north", "west_east"], mcs_case.XLAT.values),
                     lons=(["south_north", "west_east"], mcs_case.XLONG.values),)

    coords = dict(south_north=mcs_case.south_north.values, west_east=mcs_case.west_east.values, time = times)
    data = xr.Dataset(data_vars=data_vars, coords=coords)                                                                 
    data.to_netcdf('/glade/scratch/kukulies/idealized_storms/10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/4000/Idealized_MCS_' + caseID + '_' + date + '_' + resolution + '_condensation_rates.nc')       
