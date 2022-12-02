'''CryoFM(TM) is a library of functions useful for cryogenic fluid management 
   
   Default units are SI. The reference report can be accessed at: 
   https://drive.google.com/file/d/1sTNNPRgGdC4JrDt5UGz7pkyhoysBNkRZ/view
   Copyright 2022 Moran Innovation LLC. Licensed under the Apache License, 
   Version 2.0. 
   '''

import math

# DIMENSIONLESS NUMBERS

def bond(accel, radius_fs, dens_liq, dens_vap, surf_tens):
    """Bond number (ratio of acceleration to capillary forces)
    
    Keyword arguments:
    accel -- local acceleration, m/s^2
    radius_fs -- free surface radius, m^2
    dens_liq -- density of saturated liquid, kg/m^3
    dens_vap -- density of saturated vapor, kg/m^3
    surf_tens -- fluid surface tension, N/m
    """
    return (dens_liq - dens_vap) * accel * (radius_fs)**2 / surf_tens

def reynolds(velocity, dim_char, density, visc_dyn):
    """Reynolds number (ratio of inertia to viscous forces in a flowing fluid)
    
    Keyword arguments:
    velocity -- fluid velocity, m/s
    dim_char -- characteristic dimension, m
                distance from leading edge for external flow, or
                hydraulic diameter for internal flow (4*area/wetted perimeter)
    density -- fluid density, kg/m^3
    visc_dyn -- dynamic viscosity, Pa-s
    """
    return density * velocity * dim_char / visc_dyn

def grashof(accel, cte, temp_surf, temp_bulk, dim_char, visc_kin):
    """Grashof number (ratio of buoyancy to viscous force)
    
    Keyword arguments:
    accel -- local acceleration, m/s^2
    cte -- fluid coefficient of thermal expansion, 1/K (~1/T for ideal gas)
    temp_surf -- surface temperature, K
    temp_bulk -- bulk fluid temperature, K
    dim_char -- characteristic dimension, m
                vertical length for vertical flat plate, or
                hydraulic diameter for internal flow (4*area/wetted perimeter)
    visc_dyn -- dynamic viscosity, Pa-s
    """
    return (accel * cte * abs(temp_surf - temp_bulk) * dim_char**3
            / visc_kin**2)

def jakob(cp_final, temp_final, temp_init, dh_final):
    """Jakob number (ratio of sensible heat to latent heat during phase change)
    
    Keyword arguments:
    cp_final -- specific heat at final conditions, J/kg-K
    temp_final -- final fluid temperature, K
    temp_init -- initial fluid temperature, K
    dh_final -- latent heat of vaporization at final conditions, J/kg
    """
    return cp_final * (temp_final- temp_init) / dh_final

def froude(veloc, accel, length_char):
    """Froude number (ratio of the flow inertia to the external field)
    
    Keyword arguments:
    veloc -- fluid velocity, m/s
    accel -- local acceleration, m/s^2
    length_char -- characteristic length, m
    """
    return veloc**2 / (accel * length_char)

def froude_beta(veloc, accel, length_char, beta, dT):
    """Thermal expansion Froude number
    
    Keyword arguments:
    veloc -- fluid velocity, m/s
    accel -- local acceleration, m/s^2
    length_char -- characteristic length, m
    beta -- coefficient of thermal expansion, 1/k
    dT -- temperature differential, K
    """
    return froude(veloc, accel, length_char) / (beta * dT)

# PRESSURIZATION

def press_mass_JaFr(int_longdim, height_ullage, jakob, froude, k=1, m_ref=1):
    """Pressurant mass prediction based on Ja & Fr correlation
    
    Keyword arguments:
    tank_shape -- 'vertical cylinder', 'spherical', or 'horizontal cylinder'
    jakob -- Jakob number
    froude -- thermal expansion froude number
    k -- pressurant mass equation constant (default 1.0)
    mass_ref -- reference mass, kg (default 1.0)
    """
    # Validated ranges: 0.07<Ja<0.17 and 0.001<Fr<1.15
    # Ref: Ludwig and Dreyer, Cryogenics 63 (2014) 1-16
    if(int_longdim/height_ullage < 2):
        m_cn = k * jakob**(3/2) / froude**(1/3)
    else:  # modified equation for larger surface area
        m_cn = k * jakob / froude**(1/3)
    return m_cn * m_ref

# SLOSHING

def init_grad_thick(diffus_sat, time):
    """Initial thermal boundary layer gradient thickness (liquid)
    
    Keyword arguments:
    diffus_sat -- saturated liquid diffusivity
    time -- elapsed time from start of ramp pressurization to slosh start
    """
    # Ref: Ludwig, et al., Intl J of Heat and Mass Transfer 66 (2013) 223-234
    return math.sqrt(math.pi * diffus_sat * time)

def slosh_nusselt(slosh_reynolds, reynolds_crit=4000.):
    """Sloshing Nusselt number
    
    Keyword arguments:
    slosh_reynolds -- sloshing Reynolds number
    reynolds_crit -- critical sloshing Reynolds number (default 4000.)
    """
    # Ref: Ludwig, et al., Intl J of Heat and Mass Transfer 66 (2013) 223-234
    return (slosh_reynolds/reynolds_crit)**0.69

def slosh_jakob(dens_liq, cp_liq, temp_sat, temp_liq, dens_vapsat, 
                latent_heat):
    """Sloshing Jakob number
    
    Keyword arguments:
    dens_liq -- density at bulk liquid temp and press
    cp_liq -- specific heat at constant pressure at bulk liquid temp and press
    temp_sat -- saturation temperature at vapor pressure
    temp_liq -- temperature of bulk liquid
    dens_vapsat -- density of saturated vapor
    latent_heat -- latent heat of vaporization
    """
    # Ref: Ludwig, et al., Intl J of Heat and Mass Transfer 66 (2013) 223-234
    return (dens_liq * cp_liq * (temp_sat - temp_liq)/ dens_vapsat
            / latent_heat)

def slosh_reynolds(ang_wave_freq, wave_height, visc_kinliq):
    """Sloshing Reynolds number
    
    Keyword arguments:
    ang_wave_freq -- angular wave frequency (2*pi*frequency)
    wave_height -- wave height (amplitude) during sloshing
    visc_kinliq -- kinematic viscosity at bulk liquid temp and press
    """
    # Ref: Ludwig, et al., Intl J of Heat and Mass Transfer 66 (2013) 223-234
    return ang_wave_freq * wave_height**2 / visc_kinliq

def tslosh_delta(grad_init,nusselt_slosh,diffus_sat):
    """Time prior to slosh when initial thermal gradient is set up in liquid
    
    Keyword arguments:
    grad_init -- initial thermal boundary layer gradient thickness (liquid)
    nusselt_slosh -- sloshing nusselt number
    diffus_sat -- saturated liquid diffusivity"""
    # Ref: Ludwig, et al., Intl J of Heat and Mass Transfer 66 (2013) 223-234
    return grad_init**2 / (math.pi * nusselt_slosh * diffus_sat)

# TANK GEOMETRIES

def sphere_area_inter(radius, height_liq):
    """ Interface area based on liquid height in spherical tank

    Keyword arguments:
    radius -- spherical tank internal radius
    height_liq -- liquid height from tank bottom
    """
    return math.pi * (2 * radius - height_liq) * height_liq
