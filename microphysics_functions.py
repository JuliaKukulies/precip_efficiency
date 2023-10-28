"""
A collection of functions to handle the conversion and derivation of different microphyiscal quantities. 

"""
import numpy as np 

##### Physical constants #####

# specific gas constant air in J/(kg*K)  
#alternatively: R = 8.31446261815324 in J/(K*mol)
R = 287

# gravitational acceleration in m/s^2
g = 9.81
    
# specific heat at constant pressure in J/kg/K                                                                                                                                                                                
cp = 1.005

# molar mass of air in kg/mol 
#(only needed for he ideal gas law if gas constant is given in J/(mol*K) 
molar_mass_air = 28.97 / 1000

#############################

def get_air_density(pressure, temperature): 
    """
    Get air density based on the ideal gas law 

    Args: 
      pressure(np.array or xr.DataArray): air pressure field in P
      temperature(np.array or xr.DataArray): temperature field in K 

    Returns:
      rho_dry: air density for a dry air mass in kg/m3.
    
    """
    rho_dry  = pressure / (R*temperature)
    return rho_dry 


def mixing_ratio_to_density(mixing_ratio, air_density): 
    """
    Converts the mixing ratio of a substance (e.g. cloud ice) to density.
    
    Args:
       mixing_ratio(np.array or xr.DataArray): field with mixing ratio of a given hydrometer in kg/kg
       air_density(np.array or xr.DataArray): field with air density in kg/m3 at a given temperature and pressure

    Returns:
    
       density(np.array or xr.DataArray): field with air density in kg/m3 at a given temperature and pressure

    Returns:
    
       density: density field of the hydrometer in kg/m3 (same dimensions as input).
    
    """
    density= mixing_ratio * air_density
    return density




def pressure_integration(mixing_ratio):
    """
    Integrates the mixing ratio of a hydrometeor over pressure which results in kg/m2. 

    Args:
      mixing ratio(np.array): 3D or 4D field of mixing ratio where one dimensions are pressure levels

    Returns:
      integrated_mass(np.array): 2D or 3D (if time dimension) of integrated mixing ratio in kg/m2
   
    """


    
def get_saturation_vapor_pressure(temperature): 
    """
    Estimates the saturation vapor pressure for a given temperature 
    using the August-Roche-Magnus approximation. 

    Args:
       temperature: temperature or temperature field in K

    Returns:
       es: saturation vapor pressure in Pa 
    
    """
    es = 6.1094 * exp(17.625 * temperature/ (T+243.04) )
    
    return es


def get_dp_dt():
    """
    Derives the change rate in saturation vapor pressure with temperature using the Clausius-Clapeyron equation. 

    
    """


def get_condensation_rate(vertical_velocity, temperature, pressure, rho):
    """

    Estimates the condensation rate from standard model output based on saturation adjustment. 
    
    """

    

    # get change in super-saturation with temperature                                                                        
    #dqs_dT = 0
    #  latent heat of vaporization                                                                                                                                                                                           
    Lv =  2.5*10^6
    # get saturation vapor pressure


    
    
    condensation_rate = q*vertical_velocity *(dqs_dT*cp**(-1) - qs* rho/ (pressure-vapor_pressure_sat) ) * (1 + dqs_dT * Lv/cp)**(-1)

    return condensation_rate









