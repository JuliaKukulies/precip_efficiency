import xarray as xr 
from pathlib import Path 
import numpy as np 
from datetime import datetime 
import pandas as pd
from precip_efficiency import microphysics_functions as mic
import wrf
from netCDF4 import Dataset 

start = datetime(2011,7,13,0,0)
end =  datetime(2011,7,13,7,0)
times = pd.date_range(start, end, freq='1MIN') 
print(times.shape)

path = Path('/glade/scratch/kukulies/idealized_storms/10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/4000/') 
files = list(path.glob('wrfout*'))
files.sort()
print(len(files)) 

for fname in files:
    print(fname)  
    ds= xr.open_dataset(fname)
    wrfin = Dataset(fname)
    pressure_t = ds.P + ds.PB
    # get temperature field in Kelvin 
    temperature = wrf.getvar(wrfin, 'tk')
    height_t = wrf.getvar(wrfin, 'height')
    lats = ds.XLAT.squeeze().values
    lons = ds.XLONG.squeeze().values
    rho_t = mic.get_air_density(pressure_t, temperature).squeeze()
    if fname == files[0]:
        pressure = pressure_t
        height = height_t
        rho    = rho_t
    else:
        pressure = xr.concat([pressure, pressure_t], dim  = 'Time')
        height   = xr.concat([height, height_t], dim  = 'Time')
        rho      = xr.concat([rho, rho_t], dim = 'Time')

#### Save data for the case to netCDF4
data_vars = dict(pressure=(["time", "bottom_top", "south_north", "west_east"], pressure.data),
                 height = (["time", "bottom_up", "south_north", "west_east"], height.squeeze().data),
                 rho = (["time", "bottom_top", "south_north", "west_east"], rho.squeeze().data),
                 lats=(["south_north", "west_east",], lats),
                 lons=(["south_north", "west_east"], lons),) 

coords = dict(time= times, bottom_top = ds.bottom_top.values, south_north=ds.south_north.values, west_east=ds.west_east.values)
data = xr.Dataset(data_vars=data_vars, coords = coords)                                                                 
data.to_netcdf('/glade/scratch/kukulies/idealized_storms/10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/4000/Idealized_MCS_em_quarter_10_air_pressure_heights_rho.nc')







