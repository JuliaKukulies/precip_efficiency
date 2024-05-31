'''

This script joins the data tiles from the parallel output option in WRF (io_format = 102) and writes them into one netcdf file.

kukulies@ucar.edu

'''
from scipy import stats
from netCDF4 import Dataset
import glob
from tqdm import tqdm
import scipy.interpolate as spint
import numpy as np
import xarray as xr
import sys
from pathlib import Path
import math

################################################# define paths #######################################################################

wrfout_splitfiles = '/glade/campaign/mmm/c3we/mingge/WRF_DOE/1KM/Thomson_YSU/sgp_20130617_07:00:00_L1_cheyenne/'
rslout_file = '/glade/campaign/mmm/c3we/mingge/WRF_DOE/1KM/Thomson_YSU/sgp_20130617_07:00:00_L1_cheyenne/rsl.out.0000'
wrfout_splitfiles = Path('/glade/derecho/scratch/kukulies/idealized_mcs/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/250/')
wrf = Path('/glade/work/kukulies/WRF/TempScripts/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/250/')
rslout_file = wrf / 'rsl.out.0000'
savedir = Path('/glade/derecho/scratch/kukulies/idealized_mcs/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/250/combined/')

#####################################################################################################################################

# get all the necessary information from rsl output file 

file = open(rslout_file)
all_lines = file.readlines()

# get number of x and y tiles
ixy_composition = np.where(np.array([' Ntasks in X' in all_lines[ii] for ii in range(20)]) == True)[0][0]
xy_composition = all_lines[ixy_composition]
x_comp, y_comp  = [int(i) for i in xy_composition.split() if i.isdigit()]

# get domain size
idomain_size = np.where(np.array([" ids,ide,jds,jde" in all_lines[ii] for ii in range(30)]) == True)[0][0]
domain_size = all_lines[idomain_size]
_, x_size, _, y_size = [int(i) for i in domain_size.split() if i.isdigit()]

x_tile_grids = int(math.ceil(x_size/x_comp))
y_tile_grids = int(math.ceil(y_size/y_comp))
################################################################################################

# variables to write into combined output 
var_list = ['RAINNC', 'PRW_VCD', 'PRS_SDE', 'PRG_GDE', 'PRI_IDE', 'PRI_IHA', 'PRI_INU', 'QCLOUD', 'W', 'T', 'P', 'PB', 'QSNOW', 'QGRAUP', 'QICE', 'QVAPOR']

var_list = ['QRAIN']
method = 'nearest'

# get list with all time steps 
path = Path('/glade/derecho/scratch/kukulies/idealized_mcs/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/500') 
fnames = list(path.glob('wrfout*00'))
fnames.sort()
timesteps = [tt.name for tt in fnames]
timesteps.sort()

for i in np.arange(len(timesteps)):
        t = timesteps[i]
        timesteps[i]= t.replace('pr', 'process_rates')

for tstep in timesteps[32:]:
    print(tstep, flush = True)
    # write output file 
    outfile = savedir / (tstep + '_combined_tiles_qrain')
    
    if outfile.is_file() is False: 
        DATA0 = xr.open_dataset(wrfout_splitfiles / str(tstep +'_'+str(0).zfill(4) ) ) 
        DATA0 = DATA0.interp(south_north=range(y_size-1), method= method)
        DATA0 = DATA0.interp(south_north_stag=range(y_size), method= method)
        DATA0 = DATA0.interp(west_east=range(x_size-1), method= method)
        DATA0 = DATA0.interp(west_east_stag=range(y_size), method= method)
        DATA0.to_netcdf(outfile)
        print('wrote dummy output file', str(outfile), flush = True) 

        for va in tqdm(range(len( var_list))):
            data_act = DATA0[var_list[va]]
            data_act[:] = np.nan
            if len(data_act.shape) < 2:
                continue
            la_start = 0
            for la in range(y_comp):
                lo_start = 0
                for lo in range(x_comp):
                    ii = lo + int(la * x_comp)
                    # get the datatile 
                    ncid=Dataset(wrfout_splitfiles / str(tstep + '_'+str(ii).zfill(4) ), mode='r')
                    #var_list = list(DATA0.keys() ) 
                    data_tile=ncid.variables[var_list[va]][:]
                    ncid.close()
                    la_stop = la_start + data_tile.shape[-2]
                    lo_stop = lo_start + data_tile.shape[-1]

                    if len(data_act.shape) == 3:
                        data_act[:, la_start:la_stop, lo_start:lo_stop] = data_tile
                    if len(data_act.shape) == 4:
                        data_act[:, :, la_start:la_stop, lo_start:lo_stop] = data_tile
                    lo_start = lo_stop
                la_start = la_stop

            # write the data to the NetCDF
            ncid=Dataset(outfile, mode='r+')
            ncid.variables[var_list[va]][:] = data_act
            ncid.close()
        print('finished joiner program for timestep:', str(tstep), flush = True)
