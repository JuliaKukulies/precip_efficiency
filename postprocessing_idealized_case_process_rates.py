import xarray as xr 
from pathlib import Path 
import numpy as np 
from datetime import datetime 
import pandas as pd 

start = datetime(2011,7,13,0,0)
end =  datetime(2011,7,13,7,0)
times = pd.date_range(start, end, freq='1MIN') 
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
        process_t = ds[processname].where(ds[processname] > 0, 0)
        lats = ds.XLAT.squeeze().values
        lons = ds.XLONG.squeeze().values
        if fname == files[0]:
            process = process_t
        else:
            process = xr.concat([process, process_t], dim  = 'Time')
          
    #### Save data for the case to netCDF4
    data_vars = {processname:(["time", "bottom_top", "south_north", "west_east"], process.data),
                'lats':(["south_north", "west_east",], lats),
                'lons':(["south_north", "west_east"], lons),}

    coords = dict(time= times, bottom_top = ds.bottom_top.values, south_north=ds.south_north.values, west_east=ds.west_east.values)
    data = xr.Dataset(data_vars=data_vars, coords = coords)                                   
    data.to_netcdf('/glade/scratch/kukulies/idealized_storms/Idealized_MCS_em_quarter_10_' +processname+'.nc') 

