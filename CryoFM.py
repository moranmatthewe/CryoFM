'''CryoFM(TM) is a library of functions useful for cryogenic fluid management. 
   Default units are SI. Copyright 2022 Moran Innovation LLC. Licensed under 
   the Apache License, Version 2.0. The reference report can be accessed at: 
   https://drive.google.com/file/d/1sTNNPRgGdC4JrDt5UGz7pkyhoysBNkRZ/view'''


import CoolProp.CoolProp as CP


# Practice setting up a new class
class DimLess:
    def__init__(self, name, inputs, properties, value)
    self.name = name
    self.inputs = inputs
    self.properties = properties
    self.value = value

Bo.name = DimLess.name('Bond')


# Fluid properties

# List of commmon cryogenic fluids supported in CoolProp; full list at
# http://coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
def list_fluids():
    return ["Common cryogenic fluids:", "ParaHydrogen", "Hydrogen", 
            "OrthoHydrogen", "Oxygen", "Methane", "Nitrogen", "Helium", "Neon"]


# Dimensionless numbers

# Bond number based on fluid, acceleration (m/s^2), free surface diameter (m), 
# and pressure (Pa)
def bond(fluid, accel, diam_fs, press):
    # saturated liquid density, kg/m^3
    density_liquid = CP.PropsSI('D', 'P', press, 'Q', 0, fluid)
    # saturated vapor density, kg/m^3
    density_vapor = CP.PropsSI('D', 'P', press, 'Q', 1, fluid)
    # saturated surface tension,  N/m  
    surf_tension = CP.PropsSI('surface_tension', 'P', press, 'Q', 0, fluid)
    # Bond number, dimensionless 
    return (density_liquid - density_vapor) * accel * (diam_fs)**2 / surf_tension

# Reynolds number based on fluid, velocity (m/s), characteristic length or
# hydraulic diameter = 4*area/wetted perimeter (m), fluid temperature (K), and
# fluid pressure (Pa)
def reynolds(fluid, velocity, length, temp, press):
    # fluid density, kg/m^3
    density = CP.PropsSI('D', 'T', temp, 'P', press, fluid)
    # dynamic viscosity, Pa-s
    visc_dynamic = CP.PropsSI('V', 'T', temp, 'P', press, fluid)  
    # Reynolds number, dimensionless
    return density * velocity * length / visc_dynamic

# Raleigh number