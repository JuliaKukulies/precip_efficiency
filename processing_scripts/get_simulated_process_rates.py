"""
Get vertically integrated process rates (condensation and deposition rates) as output by idealized WRF simulations of MCS cases.

"""

import xarray as xr 
from pathlib import Path 
import numpy as np 
from datetime import datetime 
import pandas as pd 
import microphysics_functions as micro

times = np.arange(421)
print(times.shape)

path = Path('/glade/scratch/kukulies/idealized_storms/10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/4000/') 
files = list(path.glob('wrfout*'))
files.sort()
print(len(files)) 

process_rates = ['PRW_VCD', 'PRS_SDE', 'PRS_IDE', 'PRI_IDE', 'PRG_GDE', 'PRI_INU', 'PRI_IHA']

for processname in process_rates:
    print(processname)
    for fname in files:
        ds= xr.open_dataset(fname)
        process_t = ds[processname].where(ds[processname] > 0, 0).squeeze()
        pressure = ds.P.squeeze() + ds.PB.squeeze()
        lats = ds.XLAT.squeeze().values
        lons = ds.XLONG.squeeze().values
        integration = micro.pressure_integration(process_t.data, -pressure.data)
        if fname == files[0]:
            process = integration
        else:
            process = np.dstack((process, integration))
        print(process.shape, lats.shape, lons.shape, times.shape)
          
    #### Save data for the case to netCDF4
    data_vars = {processname:(["south_north", "west_east", "time"], process),
                'lats':(["south_north", "west_east",], lats),
                'lons':(["south_north", "west_east"], lons),}

    coords = dict(time= times, south_north=ds.south_north.values, west_east=ds.west_east.values)
    data = xr.Dataset(data_vars=data_vars, coords = coords)                                   
    data.to_netcdf('/glade/scratch/kukulies/idealized_storms/Idealized_MCS_10_' +processname+'_integrated.nc') 

