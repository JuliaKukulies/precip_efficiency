import xarray as xr 
from pathlib import Path 
import numpy as np 
from datetime import datetime 
import pandas as pd 
import microphysics_functions as micro 

path = Path('/glade/scratch/kukulies/idealized_storms/10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/4000/') 
fname = path / 'Idealized_MCS_em_quarter_10_air_pressure_heights_rho.nc'
ds = xr.open_dataset(fname)
pressure = ds.pressure
lats = ds.lats.data
lons = ds.lons.data
times = np.arange(421)
process_rates = ['PRW_VCD', 'PRS_SDE', 'PRS_IDE', 'PRI_IDE', 'PRG_GDE', 'PRI_INU', 'PRI_IHA'] 


path2 = Path('/glade/scratch/kukulies/idealized_storms/') 
for processname in process_rates:
    print(processname)
    ds= xr.open_dataset(list(path2.glob('Idealized*' + processname +  '*.nc' ))[0] ) 
    process = ds[processname].where(ds[processname] > 0, 0) 
    integration = micro.pressure_integration(process.data, -pressure.data, axis = 1)        
    #### Save data for the case to netCDF4
    data_vars = {processname:(["time", "south_north", "west_east"], integration.data),
                'lats':(["south_north", "west_east",], lats),
                'lons':(["south_north", "west_east"], lons),}
    coords = dict(time= times, south_north=ds.south_north.values, west_east=ds.west_east.values)
    data = xr.Dataset(data_vars=data_vars, coords = coords)                                   
    data.to_netcdf( path2 /  ('Idealized_MCS_em_quarter_10_' +processname+'_integrated.nc') )
