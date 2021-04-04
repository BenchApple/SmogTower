# Ben Chappell

# This file stores the functions that calculate the various forces with different parameters.

import math

# Calculates the gravitational force. This force is directed downwards
# d_part is the diameter of the particle
def grav_force(rho_part, d_part):
    g = 9.81998
    return (math.pi / 6) * rho_part * g * (d_part ** 3)

def bouyant_force(rho_part, d_part):
    g = 9.81998
    return (math.pi / 6) * rho_part * g * (d_part ** 3)

# rho_fluid is density of fluid 
# diameter is diameter of object, v_apparent is apparent velocity of object.
# mu_fluid is dynamic viscosity of fluid.
def reynolds(rho_fluid, diameter, v_apparent, mu_fluid):
    return (diameter * rho_fluid * abs(v_apparent)) / mu_fluid

# Returns the drag coefficient of air assuming stokes flow.
def drag_coeff(rho_fluid, diameter, v_apparent, mu_fluid):
    return 24 / reynolds(rho_fluid, diameter, v_apparent, mu_fluid)

# Returns the drag force given the given parameters.
# Drag force changes with time until terminal velocity because 
# v_apparent = - v_vert(t)
def drag_force(rho_fluid, diameter, v_apparent, mu_fluid):
    C_d = drag_coeff(rho_fluid, diameter, v_apparent, mu_fluid)
    return (1/2) * rho_fluid * C_d * ((math.pi / 4) * diameter ** 2) * (v_apparent ** 2)

# Describes the force on the particle due to the electric collector field
# Foce acts downwards.
# q is the charge of the particle, sigma is the charge of the plate
def electric_collector(q, sigma):
    # This is the permittivity constant
    e_0 = 8.854187 * math.pow(10, -12)
    return (q * sigma) / (4 * math.pi * e_0)

# This equation depends on the current concentration of particles c(t)
# at current height D_vert(t). This equation will assume that the given values
# are at the given instant of time desired.
# Height H is maximum height where the concentration begins to wear off.
def electric_other(q, c_oft, D_vert_oft, H):
    e_0 = 8.854187 * math.pow(10, -12)
    return ((q ** 2 * c_oft) / 2 * e_0) * (2 * D_vert_oft - H)

# Returns the mass of a particle given its density and diameter
def mass(rho_part, d_part):
    return (math.pi / 6) * rho_part * (d_part ** 3)

    

