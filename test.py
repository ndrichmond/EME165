import sympy as sym
import numpy as np
from pint import UnitRegistry
u = UnitRegistry()

T_air = (20 + 273.15) *u.kelvin #kelvin
L_1 = 1.20 * u.meter #meters, core length
L_2 = 0.25 * u.meter #meters, handle length
r_core = 0.02 * u.meter #meters, radius of core (and handle)
k_core = 0.35 *u.watt / (u.meter * u.kelvin) #W/mK
k_handle = 0.05 *u.watt / (u.meter * u.kelvin) #W/mK
R_tc = 10**-5 #m^2K/W, thermal resistance of handle-core boundary
R_tc_luke = 10**-4 #m^2K/W, thermal resistance of luke-handle boundary
T_luke = 37.1 + 273.15 #kelvin, temperature of Luke's hand
qdot_r0 = 10**6 * u.watt/(u.meter**3)
h = 50 #W/m^2k
pi = np.pi

def qgen(r1,r2,dx): #calcluate the heat generated in a volume
    return qdot_r0*pi*(dx)*((r2**2-r1**2) + ((1/(2*r_core**2))*(r1**4-r2**4)))

print(qgen(0 * u.meter, 0.02*u.meter,1*u.meter))