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
      pressure: air pressure field in P
      temperature: temperature field in K 

    """
    rho_dry  = pressure / (R*temperature)
    return rho_dry 


def mixing_ratio_to_density(mixing_ratio): 
    """
    Converts the mixing ratio of a substance (e.g. cloud ice) in kg/kg to density in kg/m3. 
 
    """

def pressure_integration(mixing_ratio):
    """
    Integrates the mixing ratio of a hydrometeor over pressure which results in kg/m2. 

    """

def get_saturation_vapor_pressure(temperature): 
    """
    Estimates the saturation vapor pressure for a given temperature 
    using the August-Roche-Magnus approximation. 

    """















