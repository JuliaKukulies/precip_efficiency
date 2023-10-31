"""
A collection of functions to handle the conversion and derivation of different microphyiscal quantities. 

"""
import numpy as np 
import math

##### Physical constants #####

# specific gas constant air in J/(kg*K)  
#alternatively: R = 8.31446261815324 in J/(K*mol)
R = 287.053
# specific gas constant for water vapor in J/(kg*K)
Rv = 461.5 

# gravitational acceleration in m/s^2
g = 9.81
    
# specific heat at constant pressure in J/kg/K
cp = 1.005

# molar mass of air in kg/mol 
#(only needed for he ideal gas law if gas constant is given in J/(mol*K) 
molar_mass_air = 28.97 / 1000
molar_mass_q =  0.01801588 # kg/mol

# latent heat of vaporization in J/kg
Lv = 2.25* 10**6 

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
    Converts the mixing ratio of a substance (e.g. cloud ice) to a density.
    
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


def vapor_pressure_to_mr(vapor_pressure):
    """
    Converts the vapor pressure of water in Pa to a mixing ratio q in kg/kg. 
    """
    mixing_ratio = R/Rv * (vapor_pressure/ (pressure-vapor_pressure))
    return mixing_ratio


def pressure_integration(mixing_ratio, pressure ):
    """
    Integrates the mixing ratio of a hydrometeor over pressure which results in kg/m2. 

    Args:
      mixing ratio(np.array): 3D or 4D field of mixing ratio where one dimensions are pressure levels
      pressure(np.array): 1D or multidimensional field with pressure levels in Pa (make sure pressure data is in right direction!) 
    
    Returns:
      integrated_mass(np.array): 2D or 3D (if time dimension) of integrated mixing ratio in kg/m2
   
    """
    return  np.trapz(mixing_ratio, pressure.squeeze, axis = 0) * 1/g


    
def get_saturation_vapor_pressure(temperature): 
    """
    Estimates the saturation vapor pressure for a given temperature 
    using the August-Roche-Magnus approximation. 

    Args:
       temperature: temperature or temperature field in K

    Returns:
       es: saturation vapor pressure in Pa 
    
    """
    # convert temperature to Celsius degrees for this equation
    temperature -= 273.15
    
    es = 6.1094 * math.exp(17.625 * temperature/ (temperature+243.04) )
    # convert from hPa to Pa
    es/= 10

    return es


def get_dp_dT(temperature):
    """
    Derives the change rate in saturation vapor pressure with temperature using the Clausius-Clapeyron equation. 

    Args:
       temperature: temperature value or field
    Returns:
       dp_dT: change rate of saturation vapor pressure with temperature in Pa/K
    
    """
    es = get_saturation_vapor_pressure(temperature)

    dp_dT = Lv/ (R*temperature**2 / es)
    return dp_dT 


def get_condensation_rate(vertical_velocity, temperature, pressure):
    """
    Estimates the condensation rate from standard model output based on saturation adjustment.

    Args:
      vertical_velocity(np.array or xr.DataArray): field with vertical velocity in m/s
      temperature(np.array or xr.DataArray): temperature field in K
      pressure(np.array or xr.DataArray): air pressure field in Pa


    Returns:
      condensation_rate: field with condensation rates for every grid point in kg/m2/s
    
    """
    # get change rate of saturation mixing ratio with temperature
    dp_dT = get_dp_dT(temperature)
    # convert saturation vapor pressure to mixing ratio
    qqs_dT = convert_vapor_pressure_to_mr(dp_dT, pressure)
    
    # get saturation vapor pressure
    es = get_saturation_pressur(temperature)
    
    # get air density
    rho = get_air_density(pressure, temperature)
    
    condensation_rate_s = g*vertical_velocity *(dqs_dT*cp**(-1) - (qs* rho)/ (pressure-es) ) * (1 + dqs_dT * Lv/cp)**(-1)

    # convert condensation rate in kg/kg per second to kg/m2/s:
    condensation_rate = convert_mixing_ratio_to_density(condensation_rate_s)

    return condensation_rate







