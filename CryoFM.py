'''CryoFM(TM) is a library of functions useful for cryogenic fluid management 
   
   Default units are SI. The reference report can be accessed at: 
   https://drive.google.com/file/d/1sTNNPRgGdC4JrDt5UGz7pkyhoysBNkRZ/view
   Copyright 2022 Moran Innovation LLC. Licensed under the Apache License, 
   Version 2.0. 
   '''

# Dimensionless numbers

def bond(accel, diam_fs, dens_liq, dens_vap, surf_tens):
    """Bond number (ratio of acceleration to capillary forces)
    
    Keyword arguments:
    accel -- local acceleration, m/s^2
    diam_fs -- free surface diameter, m^2
    dens_liq -- density of saturated liquid, kg/m^3
    dens_vap -- density of saturated vapor, kg/m^3
    surf_tens -- fluid surface tension, N/m
    """
    return (dens_liq - dens_vap) * accel * (diam_fs)**2 / surf_tens

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

def jakob(cp_final, temp_final, temp_init, dh_final)
    """Jakob number (ratio of sensible heat to latent heat during phase change)
    
    Keyword arguments:
    cp_final -- specific heat at final conditions, J/kg-K
    temp_final -- final fluid temperature, K
    temp_init -- initial fluid temperature, K
    dh_final -- latent heat of vaporization at final conditions, J/kg
    """
return cp_final * (temp_final- temp_init) / dh_final