'''CryoFM(TM) is a library of functions useful for cryogenic fluid management 
   
   Default units are SI. The reference report can be accessed at: 
   https://drive.google.com/file/d/1sTNNPRgGdC4JrDt5UGz7pkyhoysBNkRZ/view
   Copyright 2022 Moran Innovation LLC. Licensed under the Apache License, 
   Version 2.0. 
   '''


import CoolProp.CoolProp as cp


# Fluid properties

def list_fluids():
    """Commmon cryogenic fluids list used by CoolProp
    
    Keyword Arguments: none
    Ref: http://coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
    """
    return ["Common cryogenic fluids:", "ParaHydrogen", "Hydrogen", 
            "OrthoHydrogen", "Oxygen", "Methane", "Nitrogen", "Helium", "Neon"]


# Dimensionless numbers

def bond(fluid, accel, diam_fs, press):
    """Bond number (ratio of acceleration to capillary forces)
    
    Keyword arguments:
    fluid -- fluid type supported by CoolProp
    accel -- local acceleration, m/s^2
    diam_fs -- free surface diameter, m
    press -- vapor pressure, Pa
    """
    # saturated liquid density, kg/m^3
    density_liquid = cp.PropsSI('D', 'P', press, 'Q', 0, fluid)
    # saturated vapor density, kg/m^3
    density_vapor = cp.PropsSI('D', 'P', press, 'Q', 1, fluid)
    # saturated surface tension,  N/m  
    surf_tension = cp.PropsSI('surface_tension', 'P', press, 'Q', 0, fluid)
    # Bond number, dimensionless 
    return (density_liquid - density_vapor) * accel * (diam_fs)**2 / surf_tension

def reynolds(fluid, velocity, dim_char, temp, press):
    """Reynolds number (ratio of inertia to viscous forces in a flowing fluid)
    
    Keyword arguments:
    fluid -- fluid type supported by CoolProp
    velocity -- fluid velocity, m/s
    dim_char -- characteristic dimension, m
                distance from leading edge for external flow, or
                hydraulic diameter for internal flow (4*area/wetted perimeter)
    temp -- fluid temperature (K)
    press -- fluid pressure (Pa)
    """
    # fluid density, kg/m^3
    density = cp.PropsSI('D', 'T', temp, 'P', press, fluid)
    # dynamic viscosity, Pa-s
    visc_dynamic = cp.PropsSI('V', 'T', temp, 'P', press, fluid)  
    # Reynolds number, dimensionless
    return density * velocity * dim_char / visc_dynamic

# Raleigh number