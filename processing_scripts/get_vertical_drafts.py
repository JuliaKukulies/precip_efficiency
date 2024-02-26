"""
This script derives the condensation rate and other vertically integrated variables from 3D model output for idealized MCS cases at different grid spacings.

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

#### selected cases ####
#'19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'03_2011-07-16_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5'
#'23_2007-06-19_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5'
#'10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'13_2003-08-30_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'17_2011-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'18_2010-06-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'38_2007-08-04_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'46_2009-06-14_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'07_2011-07-04_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5'

#'64_2012-06-17_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'58_2009-06-08_PGW_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5'
#'41_2005-06-10_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'68_2013-07-07_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'31_2006-08-18_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'16_2002-06-11_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'34_2010-07-12_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'51_2003-06-23_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'56_2008-06-18_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5'
#'35_2004-07-02_PGW_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5')


##### 10 MCS cases #####
path = Path('/glade/campaign/mmm/c3we/prein/Idealized_MCSs/wrfout_files/WRF/')
out_path = Path('/glade/derecho/scratch/kukulies/idealized_mcs/')

present_cases = ['03', '07', '10', '13', '17', '18', '19' ,'23', '38', '46']
future_cases  = ['64', '58', '41', '68', '31', '16', '34' ,'51', '56', '35']
caseIDs = present_cases + future_cases

resolutions = [12000, 4000, 2000, 1000, 500, 250]
heights = np.arange(0,24250, 250)

for resolution in resolutions:
    resolution = str(resolution)
    for caseID in caseIDs:
        if caseID in ['03', '07', '23', '58']:
            loc = 'Loc2'
        else:
            loc = 'Loc1'

        if caseID in present_cases:
            tag = 'CTRL'
        elif caseID in future_cases:
            tag = 'PGW'

        subdir = list(path.glob(str(caseID + '*'+ tag +'*' + loc +'*Storm*TH5')))[0]
        # get times
        year = subdir.name[3:7]
        month = subdir.name[8:10]
        day = subdir.name[11:13]
        date = year + month + day
        start = datetime.datetime(int(year), int(month), int(day), 0, 0)
        end = datetime.datetime(int(year), int(month), int(day), 7, 0)
        times = pd.date_range(start, end, freq= '5MIN')
        out_fname = out_path / ('idealized_mcs_'+ tag.lower()  +'_' + caseID + '_' + date + '_' + resolution + '_vertical.nc')
        
        if not Path( out_fname ).exists():
            print('deriving data for case '+ caseID + ' ' +date+ ' with' + resolution +  ' meter resolution.', flush = True)
            pattern = subdir / resolution
            if resolution == '250':
                pattern = pattern / 'Combined'
            files = list(pattern.glob('wrfout_d01*00'))
            files.sort()
            if caseID == '58':
                del files[-1]
                times = times[:-1]]
                
            print(len(files), times.shape[0])
            assert times.shape[0]  == len(files)
            for idx, fname in enumerate(files):
                print(subdir.name, fname)
                #### Read in data for timestep
                mcs_case = xr.open_dataset(fname).squeeze()
                wrfin = Dataset(fname)
                # get variables 
                iwc = mcs_case.QSNOW +  mcs_case.QICE + mcs_case.QGRAUP
                lwc =mcs_case.QRAIN + mcs_case.QCLOUD 

                # mean over horizontal dimension
                iwc = iwc.mean(["south_north","west_east"]).data
                lwc = lwc.mean(["south_north", "west_east"]).data
                try:
                    radar_reflectivity = mcs_case.DBZ.squeeze()
                except:
                    radar_reflectivity = xr.open_dataset(files[idx-1]).DBZ.squeeze()

                # get vertical velocity on mass points instead of staggered
                vertical_velocity = wrf.getvar(wrfin, 'wa')
                qcloud = mcs_case.QCLOUD

                # consider convective mass flux and vertical velocities (only where we have hydrometeors) 
                updraft_intensity = vertical_velocity.where( (radar_reflectivity  > 0)     &   (vertical_velocity > 1.5)).mean(["south_north", "west_east"]).data
                downdraft_intensity = vertical_velocity.where( (radar_reflectivity  > 0)   &  (vertical_velocity < -1.5)).mean(["south_north", "west_east"]).data

                #### Derive condensation rate (in kg/kg/s) 
                condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temp, pressure)
                ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
                condensation_masked = condensation_rate_s.where( qcloud > 0, 0)
                # we are only interested in positive values, as negative values are evaporation  
                condensation = condensation_masked.where(condensation_masked > 0, 0 ).mean(["south_north", "west_east"]).data
                
                #### Concatenate timesteps: ####                
                if fname == files[0]: 
                    tiwc = iwc
                    tlwc = lwc
                    condensation_rate = condensation
                    updrafts = updraft_intensity
                    downdrafts = downdraft_intensity  
                else:
                    tiwc  = np.dstack((tiwc, iwc))
                    tlwc = np.dstack((tlwc, lwc))
                    condensation_rate = np.dstack((condensation_rate, condensation))
                    updrafts = np.dstack((updrafts, updraft_intensity ))
                    downdrafts = np.dstack((downdrafts, downdraft_intensity ))

            #### Save data for the case to netCDF4 #### 
            data_vars = dict(tiwc=(["bottom_top", "time"], tiwc),
                             tlwc=(["bottom_top", "time"], tlwc),
                             condensation_rate=(["bottom_up", "time"], condensation_rate),
                             updrafts = (["bottom_up", "time"], updrafts),
                             downdrafts = (["bottom_up", "time"], downdrafts),
                             levels=(["bottom_top"], heights),) 
            
            coords = dict(levels=mcs_case.bottom_top.values, time = times)
            data = xr.Dataset(data_vars=data_vars, coords=coords)
            data.to_netcdf(out_fname)
            data.close()
            mcs_case.close()
        else:
            print(caseID, resolution, flush = True)
            continue
