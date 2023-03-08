# Author: Seamus Flannery
# Numbers from NASA Jupiter Fact Sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
# and

Jupiter2 = {
    "Mass": 1.89813 * 10 ** 27,  # kg
    "Radius": 69911,   # km, equatorial radius
    "surface_gravity": 24.79,  # m/s^2, acceleration due to gravity at surface level

    # Density Profile
    "surface_density": 0.065,  # g/cm^3
    "scale_height": 15.9,  # km, scale height for density exponential

    # Volumetric mix of the atmosphere, (included an extra 0.05% O)
    # Molecular hydrogen (H2) - 89.8% +/- 2.0%; Helium (He) - 10.2% +/- 2.0%
    # Methane (CH4) - 0.3% +/- 0.1%; Ammonia (NH3) - 0.026% +/-0.004%; Hydrogen Deuteride (HD) - 0.0028% +/- 0.001%;
    # Ethane (C2H6) - 0.00058% +/- 0.00015%; Water (H2O) - 0.0004% (varies with pressure)
    "mix_N": 0.0,
    "mix_O": 0.0,
    "mix_Ar": 0.0,
    "mix_C": 0.0,
    "mix_H": 0.87,
    "mix_He": 0.13

}

if __name__ == '__main__':
    print(Jupiter["mix_N"]+Jupiter["mix_O"]+Jupiter["mix_Ar"]+Jupiter["mix_C"])
