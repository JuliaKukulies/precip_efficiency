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
caseIDs = future_cases

resolutions = [12000, 4000, 2000, 1000, 500, 250]
resolutions  = ['12000nc']

for resolution in resolutions:
    resolution = str(resolution)
    for caseID in caseIDs:
        if caseID in ['03', '07', '23', '58']:
            loc = 'Loc2'
        else:
            loc = 'Loc1'
        if caseID == '03':
            resolution = '12000'
        else:
            resolution = '12000nc'

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
        print(start, end, flush = True)
        times = pd.date_range(start, end, freq= '5MIN')
        out_fname = out_path / ('idealized_mcs_'+ tag.lower()  +'_' + caseID + '_' + date + '_' + resolution + '_vars.nc')
        
        if not Path( out_fname ).exists():
            print('deriving data for case '+ caseID + ' ' +date+ ' with' + resolution +  ' meter resolution.', flush = True)
            pattern = subdir / resolution
            if resolution == '250':
                pattern = pattern / 'Combined'
            files = list(pattern.glob('wrfout*'))
            files.sort()
            if caseID == '58':
                del files[-1]
                times = times[:-1]          
                
            print(len(files), times.shape[0])
            assert times.shape[0]  == len(files)
            for idx, fname in enumerate(files):
                print(subdir.name, fname)
                #### Read in data for timestep
                mcs_case = xr.open_dataset(fname).squeeze()
                wrfin = Dataset(fname)
                # get variables 
                iwc = mcs_case.QSNOW +  mcs_case.QICE +  mcs_case.QGRAUP
                lwc =mcs_case.QRAIN + mcs_case.QCLOUD
                precip = mcs_case.RAINNC
                tprecip = mcs_case.TOTAL_PRECIP
                try:
                    radar_reflectivity = mcs_case.DBZ.squeeze()
                except:
                    radar_reflectivity = xr.open_dataset(files[idx-1]).DBZ.squeeze()
                # get vertical velocity on mass points instead of staggered
                vertical_velocity = wrf.getvar(wrfin, 'wa')
                temp = wrf.getvar(wrfin, 'temp')
                qcloud = mcs_case.QCLOUD
                pressure = mcs_case.P + mcs_case.PB
                # integate iwc and lwp over pressure
                iwp = micro.pressure_integration(iwc.data, -pressure.data)
                lwp = micro.pressure_integration(lwc.data, -pressure.data)
                
                # consider convective mass flux and vertical velocities (only where we have hydrometeors)
                updraft_intensity = vertical_velocity.where( (radar_reflectivity  > 0 )  &   ( vertical_velocity > 0 )).mean('bottom_top')
                downdraft_intensity = vertical_velocity.where( (radar_reflectivity  > 0 ) &  ( vertical_velocity < 0 )).mean('bottom_top')
                
                rho = micro.get_air_density(pressure, temp) 
                updrafts_convective = vertical_velocity.where( (radar_reflectivity  > 0 ) &   ( vertical_velocity > 2 )) * rho
                downdrafts_convective = vertical_velocity.where( (radar_reflectivity  > 0 ) & ( vertical_velocity < -2)) * rho
                mass_flux_up = np.nansum(updrafts_convective, axis = 0)
                mass_flux_down = np.nansum(downdrafts_convective, axis = 0)

                #### Derive condensation rate (in kg/kg/s)
                condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temp, pressure)
                ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
                condensation_masked = condensation_rate_s.where( qcloud > 0, 0)
                # we are only interested in positive values, as negative values are evaporation 
                condensation_positive = condensation_masked.where(condensation_masked > 0, 0 ).data
                ### integrate over pressure levels to get kg/m2/s
                condensation = micro.pressure_integration(condensation_positive,-pressure.data)
                condensation_compare = micro.pressure_integration(condensation_rate_s.where(condensation_rate_s > 0),-pressure.data)
                
                #### Concatenate timesteps:  
                if fname == files[0]:
                    tiwp = iwp
                    tlwp = lwp
                    condensation_rate = condensation
                    condensation_rate2 = condensation_compare                    
                    total_precip = tprecip
                    vertical_mass_flux_up = mass_flux_up
                    vertical_mass_flux_down = mass_flux_down
                    updrafts = updraft_intensity
                    downdrafts = downdraft_intensity  
                else:
                    total_precip = np.dstack((total_precip, tprecip)) 
                    tiwp  = np.dstack((tiwp, iwp ))
                    tlwp = np.dstack((tlwp, lwp))
                    condensation_rate = np.dstack((condensation_rate, condensation))
                    vertical_mass_flux_up = np.dstack((vertical_mass_flux_up, mass_flux_up ))
                    vertical_mass_flux_down = np.dstack((vertical_mass_flux_down, mass_flux_down ))
                    updrafts = np.dstack((updrafts, updraft_intensity ))
                    downdrafts = np.dstack((downdrafts, downdraft_intensity ))
                    condensation_rate2 = np.dstack(( condensation_rate2 , condensation_compare)) 

            #### Save data for the case to netCDF4
            data_vars = dict(tiwp=(["south_north", "west_east", "time"], tiwp),
                             tlwp=(["south_north", "west_east", "time"], tlwp),
                             total_precip=(["south_north", "west_east", "time"],total_precip),
                             condensation_rate=(["south_north", "west_east", "time"], condensation_rate),
                             condensation_rate2=(["south_north", "west_east", "time"], condensation_rate2),
                             vertical_mass_flux_up = (["south_north", "west_east", "time"], vertical_mass_flux_up),
                             vertical_mass_flux_down = (["south_north", "west_east", "time"], vertical_mass_flux_down),
                             updrafts = (["south_north", "west_east", "time"], updrafts),
                             downdrafts = (["south_north", "west_east", "time"], downdrafts),
                             lats=(["south_north", "west_east"], mcs_case.XLAT.values),
                             lons=(["south_north", "west_east"], mcs_case.XLONG.values),)
            
            coords = dict(south_north=mcs_case.south_north.values, west_east=mcs_case.west_east.values, time = times)
            data = xr.Dataset(data_vars=data_vars, coords=coords)
            data.to_netcdf(out_fname)
            data.close()
            mcs_case.close()
        else:
            print(caseID, resolution, flush = True )
            continue
