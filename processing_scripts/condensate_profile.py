"""
Total ice and liquid content as well as condensation and deposition rate in a vertical cross section for idealized MCS case. 

"""

import xarray as xr 
from pathlib import Path 
import numpy as np 
from datetime import datetime 
import pandas as pd 
from microphysics import microphysics_functions as micro
import wrf 
from netCDF4 import Dataset


scales = ['12000', '4000', '2000', '1000', '500']


for scale in scales:
    print('get profiles for ', scale)

    times = np.arange(85)
    path = Path(('/glade/derecho/scratch/kukulies/idealized_mcs/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/' + scale ) ) 
    files = list(path.glob('wrfout*pr*'))
    files.sort()
    assert len(files) == times.size
    
    for fname in files:
        ds= xr.open_dataset(fname).squeeze()
        wrfin = Dataset(fname)
        lons = ds.XLONG.mean('south_north').data
        pressure_t = wrf.getvar(wrfin, 'pres')
        temperature = wrf.getvar(wrfin, 'tk')
        vertical_velocity = wrf.getvar(wrfin, 'wa')
        rho= micro.get_air_density(pressure_t, temperature).squeeze()
        
        #### Derive condensation rate (in kg/kg/s)                                                                     
        condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temperature, pressure_t)
        ### apply qcloud mask, because equation is conditional for grid cells with actual condensate                   
        condensation_masked = condensation_rate_s.where( ds.QCLOUD > 0)
        # we are only interested in positive values, as negative values are evaporation                                
        condensation_positive = condensation_masked.where(condensation_masked > 0 )
        condensation_rate  = condensation_positive.mean('south_north').values
        
        total_ice = (ds.QGRAUP + ds.QSNOW + ds.QICE ) #* rho 
        total_liquid = (ds.QRAIN + ds.QCLOUD) #* rho
        total_ice_cloud = total_ice.mean('south_north').values 
        total_liquid_cloud = total_liquid.mean('south_north').values
        total_ice = total_ice.where(total_ice > 0 ).mean('south_north').values 
        total_liquid = total_liquid.where(total_liquid > 0 ).mean('south_north').values
        
        total_condensation = ds.PRW_VCD.where(ds.PRW_VCD > 0).mean('south_north').values
        total_evaporation = ds.PRW_VCD.where(ds.PRW_VCD < 0).mean('south_north').values
        total_deposition = np.nansum([ds.PRS_SDE.where(ds.PRS_SDE > 0).mean('south_north').values,  ds.PRS_IDE.where(ds.PRS_IDE >0).mean('south_north').values,  ds.PRI_IDE.where(ds.PRI_IDE > 0).mean('south_north').values,  ds.PRG_GDE.where(ds.PRG_GDE > 0).mean('south_north', skipna= True).values,  ds.PRI_INU.where(ds.PRI_INU > 0).mean('south_north').values,  ds.PRI_IHA.where(ds.PRI_IHA >0).mean('south_north').values], axis = 0 ) 
        total_sublimation = np.nansum([ds.PRS_SDE.where(ds.PRS_SDE < 0).mean('south_north', skipna= True).values,  ds.PRS_IDE.where(ds.PRS_IDE < 0).mean('south_north', skipna= True).values,  ds.PRI_IDE.where(ds.PRI_IDE < 0).mean('south_north', skipna= True).values,  ds.PRG_GDE.where(ds.PRG_GDE < 0).mean('south_north', skipna= True).values, ds.PRI_INU.where(ds.PRI_INU < 0).mean('south_north', skipna= True).values,  ds.PRI_IHA.where(ds.PRI_IHA < 0).mean('south_north', skipna= True).values], axis = 0 )

        # include temperature profile 
        temp_profile = temperature.mean('south_north', skipna= True).values 
        
        if fname == files[0]:
            ice = total_ice 
            liquid = total_liquid
            ice_cloud = total_ice_cloud
            liquid_cloud = total_liquid_cloud
            condensation = total_condensation
            evaporation = total_evaporation
            deposition = total_deposition
            sublimation = total_sublimation
            temp = temp_profile
            condensation_derived = condensation_rate
        else:
            ice = np.dstack((ice, total_ice))
            liquid = np.dstack( (liquid, total_liquid))
            ice_cloud = np.dstack((ice_cloud, total_ice_cloud))
            liquid_cloud = np.dstack( (liquid_cloud, total_liquid_cloud))
            condensation =  np.dstack( (condensation, total_condensation)  )
            condensation_derived =  np.dstack( (condensation_derived, condensation_rate)  )
            evaporation =  np.dstack( (evaporation, total_evaporation)  )
            deposition  = np.dstack( (deposition, total_deposition) )
            sublimation = np.dstack( (sublimation, total_sublimation))
            temp = np.dstack( (temp, temp_profile))
            
    #### Save data for the case to netCDF4
    data_vars = {'total_ice_condensate': (["bottom_top", "west_east", "time"], ice),
                 'total_liquid_condensate':(["bottom_top", "west_east", "time"], liquid),
                 'total_ice_cloud': (["bottom_top", "west_east", "time"], ice_cloud),
                 'total_liquid_cloud':(["bottom_top", "west_east", "time"], liquid_cloud),
                 'total_condensation': (["bottom_top", "west_east", "time"], condensation),
                 'condensation_derived': (["bottom_top", "west_east", "time"], condensation_derived),
                 'total_deposition': (["bottom_top", "west_east", "time"], deposition),
                 'total_sublimation': (["bottom_top", "west_east", "time"], sublimation),
                 'total_evaporation': (["bottom_top", "west_east", "time"], evaporation),
                 'temperature': (["bottom_top", "west_east", "time"], temp),                
                 'lons':(["west_east"] , lons),}

    coords = dict(time= times, bottom_top= ds.bottom_top.values, west_east=ds.west_east.values)
    data = xr.Dataset(data_vars=data_vars, coords = coords)                                   
    data.to_netcdf(('/glade/derecho/scratch/kukulies/idealized_mcs/idealized_mcs_19_vertical-profile_' + scale + '.nc'))

