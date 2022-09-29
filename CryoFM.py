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

# Raleigh number