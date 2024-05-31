
"""
This script derives the estimated and simulated condensation rates as well as accumulated precip for different temporal output from idealized MCS simulations.

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



##### MCS cases #####
parent_path = Path('/glade/derecho/scratch/kukulies/idealized_mcs/')

cases =['19_2011-07-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '03_2011-07-16_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5',
        '23_2007-06-19_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5',
        '10_2009-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '13_2003-08-30_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '17_2011-06-27_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '18_2010-06-13_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '38_2007-08-04_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '46_2009-06-14_CTRL_Midwest_-Loc1_MCS_Storm-Nr_JJA-8-TH5',
        '07_2011-07-04_CTRL_Midwest_-Loc2_MCS_Storm-Nr_JJA-8-TH5', ]

grid_spacings = ['12000', '4000', '2000', '1000', '500']
grid_spacings = [ '250']

for case in [cases[0]]:
    for delta_x in grid_spacings:
        fnames = dict()
        path = parent_path / case / delta_x / 'combined' 

        # file lists with different temporal output 
        fnames['5min'] = list(path.glob('wrfout*tiles'))
        fnames['5min'].sort()
        #fnames['10min'] = fnames['5min'][::2]
        #fnames['15min'] =  fnames['5min'][::3]
        #fnames['25min'] = fnames['5min'][::5]
        #fnames['35min'] = fnames['5min'][::7]
        #fnames['45min'] = fnames['5min'][::9]
        #fnames['55min'] = fnames['5min'][::11]
        #fnames['1H'] = fnames['5min'][::12]

        year = case[3:7]
        month = case[8:10]
        day = case[11:13]

        # get time series for case 
        date = year + month + day
        print(case, year, month, day, flush = True)
        start = datetime.datetime(int(year), int(month), int(day), 3, 0)
        end = datetime.datetime(int(year), int(month), int(day), 6, 55)
        
        for key in list(fnames.keys()):
            print(key, flush = True)
            times = pd.date_range(start, end , freq = key)
            files = fnames[key]
            print(times.size, 'timesteps', flush = True)
            assert times.size == len(files)
            print(len(files),' files detected for output timestep: ', key ,  flush = True)

            # loop through output files
            for fname in files:
                print(fname, flush = True)
                #### Read in data for that timestep
                mcs_case = xr.open_dataset(fname).squeeze()
                print(mcs_case.QCLOUD.mean(), flush = True)
                wrfin = Dataset(fname)
                #vertical_velocity = wstagger_to_mass(mcs_case.W)
                #temp = mcs_case.T 
                vertical_velocity = wrf.getvar(wrfin, 'wa')
                temp = wrf.getvar(wrfin, 'tk')
                qcloud = mcs_case.QCLOUD
                #pressure = mcs_case.P + mcs_case.PB
                pressure = wrf.getvar(wrfin, 'pres')

                #### Derive condensation rate from hourly variables instead of minute-output (in kg/kg/s)
                condensation_rate_s= micro.get_condensation_rate(vertical_velocity, temp, pressure)
                rho = micro.get_air_density(pressure, temp).squeeze()

                ### apply qcloud mask, because equation is conditional for grid cells with actual condensate 
                condensation_cloud = condensation_rate_s.where(qcloud > 0, 0 )
                # we are only interested in positive values, as negative values are evaporation 
                condensation_masked = condensation_cloud.where(condensation_cloud > 0, 0 ).data
                ### integrate over pressure levels to get kg/m2/s
                condensation = micro.pressure_integration(condensation_masked,-pressure.data)

                # get air density
                rho = micro.get_air_density(pressure.squeeze().data, temp.squeeze().data)
                # vertical integration of process rate
                prw_vcd = mcs_case.PRW_VCD.where(mcs_case.PRW_VCD > 0, 0 ).data
                prwvcd_integrated = micro.pressure_integration(prw_vcd , -pressure.data)

                evapo= mcs_case.PRW_VCD.where(mcs_case.PRW_VCD < 0, 0 ).data
                evapo_integrated = micro.pressure_integration(evapo , -pressure.data)

                # getting the deposition rates
                total_deposition = mcs_case.PRS_SDE.where(mcs_case.PRS_SDE > 0, 0 ) + mcs_case.PRS_IDE.where(mcs_case.PRS_IDE >0, 0) + mcs_case.PRI_IDE.where(mcs_case.PRI_IDE > 0, 0) + mcs_case.PRG_GDE.where(mcs_case.PRG_GDE > 0, 0 ) + mcs_case.PRI_INU.where(mcs_case.PRI_INU > 0, 0) +  mcs_case.PRI_IHA.where(mcs_case.PRI_IHA > 0, 0)
                # integrate these
                deposition_integrated = micro.pressure_integration(total_deposition, -pressure.data)

                # accumulate precip, take instaneous values for w,p,T 
                if fname == files[0]: 
                    surface_precip = mcs_case.RAINNC
                    prwvcd = prwvcd_integrated
                    condensation_rate = condensation
                    evaporation_rate = evapo_integrated
                    deposition_rate = deposition_integrated
                else:
                    surface_precip = np.dstack((surface_precip,  mcs_case.RAINNC  )) 
                    prwvcd = np.dstack((prwvcd,  prwvcd_integrated ))
                    condensation_rate = np.dstack((condensation_rate, condensation ))
                    evaporation_rate = np.dstack((evaporation_rate, evapo_integrated))
                    deposition_rate =  np.dstack((deposition_rate, deposition_integrated ))

            #### Save data for the case to netCDF4
            data_vars = dict(evaporation_rate=(["south_north", "west_east", "time"], evaporation_rate),
                             surface_precip=(["south_north", "west_east", "time"], surface_precip),
                             deposition_rate=(["south_north", "west_east", "time"], deposition_rate),
                             condensation_rate=(["south_north", "west_east", "time"], condensation_rate),
                             prwvcd=(["south_north", "west_east", "time"], prwvcd),
                             lats=(["south_north", "west_east"], mcs_case.XLAT.values),
                             lons=(["south_north", "west_east"], mcs_case.XLONG.values),)
            
            coords = dict(south_north=mcs_case.south_north.values, west_east=mcs_case.west_east.values, time = times) 
            data = xr.Dataset(data_vars=data_vars, coords=coords)                                                                 
            data.to_netcdf('/glade/derecho/scratch/kukulies/idealized_mcs/'+ case + '/' + delta_x + '/idealized_mcs_condensation_temporal_dependence_' + key + '_error.nc')


