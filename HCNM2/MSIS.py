# Author: Nathaniel Ruhl
# The function in this script can be used to extract an atmospheric density profile from the pymsis python library

import numpy as np
from pymsis import msis

# This function gets msis density profile directly from pymsis
def get_pymsis_density(datetime, lon, lat, f107, ap, version, alts=np.arange(0, 1000, 1.0)):
    f107a = f107
    aps = [[ap] * 7]
    output = msis.run(datetime, lon, lat, list(alts), f107, f107a, aps, version=version)

    densities = output[0, 0, 0, :, 0]  # kg/m^3
    densities = densities / 1000  # g/cm^3
    return alts, densities
