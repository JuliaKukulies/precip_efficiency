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


def wstagger_to_mass(W):
    """                                                                                                         
    W are the data on the top and bottom of a grid box                                                          
    A simple conversion of the stagger grid to the mass points.                                                 
                                                                                                                
    (row_j1+row_j2)/2 = masspoint_inrow                                                                         
                                                                                                                
    Input:                                                                                                      
        Wgrid with size (##+1)                                                                                  
    Output:                                                                                                     
        W on mass points with size (##)                                                                         
    """
    # create the first column manually to initialize the array with correct dimensions                          
    W_masspoint = (W[0, :, :]+W[1, :, :])/2. # average of first and second column    
    W_masspoint = np.expand_dims(W_masspoint, 0)

    W_num_levels = int(W.shape[0])-1 # we want one less level than we have                                      

    # Loop through the rest of the rows                                                                         
    # We want the same number of rows as we have columns.                                                       
    # Take the first and second row, average them, and store in first row in V_masspoint                        
    for lev in range(1,W_num_levels):
        lev_avg = (W[lev, :, :]+W[lev+1, :, :])/2.
        lev_avg = np.expand_dims(lev_avg, 0)
        W_masspoint = np.vstack((W_masspoint,lev_avg))
    return W_masspoint


scales = ['12000', '4000', '2000', '1000', '500']
scales = ['250']

for scale in scales:
    print('get profiles for ', scale)
    times = np.arange(48)
    path = Path(('/glade/derecho/scratch/kukulies/idealized_mcs/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/' + scale + '/combined/' ) ) 
    files = list(path.glob('wrfout*tiles'))
    files.sort()


    path = Path('/glade/derecho/scratch/kukulies/idealized_mcs/19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5/250/combined/')
    liquid_files = list(path.glob('*qrain'))
    liquid_files.sort()

    print(len(files), len(liquid_files), times.size)
    #assert len(files) == times.size
    
    for ii, fname in enumerate(files):
        ds= xr.open_dataset(fname).squeeze()
        wrfin = Dataset(fname)
        lons = ds.XLONG.mean('south_north').data
        pressure_t = wrf.getvar(wrfin, 'pres')
        #pressure_t = ds.P + ds.PB
        temperature = ds.T
        #vertical_velocity = wstagger_to_mass(ds.W)
        temperature = wrf.getvar(wrfin, 'tk')
        vertical_velocity = wrf.getvar(wrfin, 'wa')
        rho= micro.get_air_density(pressure_t, temperature).squeeze()

        # get QRAIN from additional file
        ds2 = xr.open_dataset(liquid_files[ii])
        QRAIN = ds2.QRAIN.squeeze()

        #### Derive condensation rate (in kg/kg/s)                                  
        condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temperature, pressure_t)
        ### apply qcloud mask, because equation is conditional for grid cells with actual condensate                   
        condensation_masked = condensation_rate_s.where( ds.QCLOUD > 0)
        # we are only interested in positive values, as negative values are evaporation                                
        condensation_positive = condensation_masked.where(condensation_masked > 0 )
        condensation_rate  = condensation_positive.mean('south_north').values
        
        total_ice = (ds.QGRAUP + ds.QSNOW + ds.QICE ) #* rho 
        total_liquid = (QRAIN + ds.QCLOUD) #* rho
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
    data.to_netcdf(('/glade/derecho/scratch/kukulies/idealized_mcs/idealized_mcs_19_vertical-profile_' + scale + '_qrain.nc'))

