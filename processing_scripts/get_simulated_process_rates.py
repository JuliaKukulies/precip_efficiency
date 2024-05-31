"""
Get vertically integrated process rates (condensation and deposition rates) as output by idealized WRF simulations of MCS cases.

"""
from netCDF4 import Dataset
import xarray as xr 
from pathlib import Path 
import numpy as np 
import datetime 
import pandas as pd 
from microphysics import microphysics_functions as micro
import wrf

path = Path('/glade/derecho/scratch/kukulies/idealized_mcs/')


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

# PRW_VCD
process_rates = [ 'PRS_SDE', 'PRS_IDE', 'PRI_IDE', 'PRG_GDE', 'PRI_INU', 'PRI_IHA']

grid_spacings = ['12000', '4000', '2000', '1000','500'] 
grid_spacings = ['4000']

for case in cases:
    for grid_spacing in grid_spacings: 
        year = case[3:7]
        month = case[8:10]
        day = case[11:13]
        date = year + month + day
        start = datetime.datetime(int(year), int(month), int(day), 0, 0)
        end = datetime.datetime(int(year), int(month), int(day), 7, 0)
        times = pd.date_range(start, end, freq= '5MIN')

        files = list((path / case / grid_spacing).glob('wrfout*pr*'))
        files.sort()
        assert len(files) == times.size
        print(case, grid_spacing, len(files), flush = True)

        for processname in process_rates:
            for fname in files:
                wrfin = Dataset(fname)
                ds= xr.open_dataset(fname)
                print(processname, fname, flush = True)
                process_t = ds[processname].where(ds[processname] < 0, 0).squeeze()
                process_t = - process_t
                pressure = ds.P.squeeze() + ds.PB.squeeze()
                pressure = wrf.getvar(wrfin, 'pres')
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
            data.to_netcdf('/glade/derecho/scratch/kukulies/idealized_mcs/idealized_mcs_'+case+'_' + grid_spacing + '_E_'+  processname+'_integrated.nc') 

