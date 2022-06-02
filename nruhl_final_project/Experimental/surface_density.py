# Author: Nathaniel Ruhl

# In some cases, surface pressure and temp is given for planets. This notebook is used to convert to surface density using the ideal gas law, p = rho R_spec T.

R = 8.3124    # Ideal gas constant, [J/K mol]

# M0: Molar Mass [kg/mol]
# p0: Surface pressure [Pa]
# T0: Surface temp [K]
def rho0(p0, T0, M0):
    rh0 = (M0*p0)/(R*T0)
    return rh0
