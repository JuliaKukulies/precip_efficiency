"""
Total ice and liquid content as well as condensation and deposition rate in a vertical cross section for idealized MCS case. 

"""

import xarray as xr 
from pathlib import Path 
import numpy as np 
from datetime import datetime 
import pandas as pd 
import microphysics_functions as micro
import wrf 
from netCDF4 import Dataset
times = np.arange(421)

path = Path('/glade/scratch/kukulies/idealized_storms/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/4000/') 
files = list(path.glob('wrfout*'))
files.sort()
assert len(files) == times.size

for fname in files:
    ds= xr.open_dataset(fname).squeeze()
    wrfin = Dataset(fname)
    lons = ds.XLONG.mean('south_north').data
    pressure_t = wrf.getvar(wrfin, 'pres')
    temperature = wrf.getvar(wrfin, 'tk')
    rho= micro.get_air_density(pressure_t, temperature).squeeze()
    total_ice = (ds.QGRAUP + ds.QSNOW + ds.QICE ) * rho 
    total_liquid = (ds.QRAIN + ds.QCLOUD) * rho
    total_ice = total_ice.sum('south_north').values
    total_liquid = total_liquid.sum('south_north').values
    total_condensation = ds.PRW_VCD.where(ds.PRW_VCD > 0).sum('south_north').values
    total_deposition = ds.PRS_SDE.where(ds.PRS_SDE >= 0, 0) + ds.PRS_IDE.where(ds.PRS_IDE >= 0,0) + ds.PRI_IDE.where(ds.PRI_IDE >= 0,0) + ds.PRG_GDE.where(ds.PRG_GDE >= 0,0) + ds.PRI_INU.where(ds.PRI_INU >= 0,0) + ds.PRI_IHA.where(ds.PRI_IHA >= 0,0) 
    total_deposition = total_deposition.sum('south_north').values 

    if fname == files[0]:
        ice = total_ice 
        liquid = total_liquid
        condensation = total_condensation
        deposition = total_deposition
    else:
        ice = np.dstack((ice, total_ice))
        liquid = np.dstack( (liquid, total_liquid))
        condensation =  np.dstack( (condensation, total_condensation)  )
        deposition  = np.dstack( (deposition, total_deposition) )
        
        
#### Save data for the case to netCDF4
data_vars = {'total_ice_condensate': (["bottom_top", "west_east", "time"], ice),
             'total_liquid_condensate':(["bottom_top", "west_east", "time"], liquid),
             'total_condensation': (["bottom_top", "west_east", "time"], condensation),
             'total_deposition': (["bottom_top", "west_east", "time"], deposition),
             'lons':(["west_east"] , lons),}

coords = dict(time= times, bottom_top= ds.bottom_top.values, west_east=ds.west_east.values)
data = xr.Dataset(data_vars=data_vars, coords = coords)                                   
data.to_netcdf('/glade/scratch/kukulies/idealized_storms/Idealized_MCS_19_profile.nc') 

