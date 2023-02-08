# Author: Nathaniel Ruhl
# Cross Section formulas created by Balucinska-Church and McCammon (1992) (xsects in units of cm^2/g)

import numpy as np

# Valid energy range is 0.03 keV to 10 keV
# Cross Sections are for elemental Oxygen, Nitrogen, and Argon
# I put these functions in this class mainly for namespacing purposes
class BCM:
    def __init__(self):
        return

    # This function can either takes the energy as an array of values or just a single value. Each of the other functions only take in single values for energy at a time
    @staticmethod
    def get_total_xsect(mean_energy_kev, mix_N, mix_O, mix_Ar, mix_C, mix_H, mix_He):
        energy_ev = mean_energy_kev * 1000
        if isinstance(energy_ev, np.ndarray):
            O_xsect = np.array(list(map(BCM.oxygen_xsect, energy_ev)))
            N_xsect = np.array(list(map(BCM.nitrogen_xsect, energy_ev)))
            Ar_xsect = np.array(list(map(BCM.argon_xsect, energy_ev)))
            C_xsect = np.array(list(map(BCM.carbon_xsect, energy_ev)))
            H_xsect = np.array(list(map(BCM.hydrog_xsect(), energy_ev)))
            He_xsect = np.array(list(map(BCM.helium_xsect(), energy_ev)))
            xsect_total = mix_O * O_xsect + mix_N * N_xsect + mix_Ar * Ar_xsect + mix_C * C_xsect + mix_H * H_xsect + mix_He * He_xsect
        else:
            xsect_total = mix_O * BCM.oxygen_xsect(energy_ev) + mix_N * BCM.nitrogen_xsect(energy_ev) \
                + mix_Ar * BCM.argon_xsect(energy_ev) + mix_C * BCM.carbon_xsect(energy_ev) + mix_H *  \
                BCM.hydrog_xsect(energy_ev) + mix_He * BCM.helium_xsect(energy_ev)
        return xsect_total

    @staticmethod
    def oxygen_xsect(energy_ev):
        Elog = np.log(energy_ev)
        if energy_ev < 531.7:
            X = 2.57264 + (10.9321 * Elog) + \
                (-1.79383*Elog**2) + (0.102619*Elog**3)
        else:
            X = 16.53869 + (0.6428144 + 3.)*Elog - 0.3177744 * \
                Elog**2 + 7.9471897e-3 * (Elog ** 3)
        oxygen_xsect = np.exp(X)/(energy_ev**3)
        return oxygen_xsect

    @staticmethod
    def nitrogen_xsect(energy_ev):
        Elog = np.log(energy_ev)
        if energy_ev < 401.0:
            X = 9.24058 + (7.02985 * Elog) + (-1.08849 * Elog *
                                              Elog) + (0.0611007 * Elog * Elog * Elog)
        else:
            X = -13.0353 + (15.4851 * Elog) + (-1.89502 *
                                               Elog**2) + (0.0769412*Elog**3)
        nitrogen_xsect = np.exp(X)/(energy_ev**3)
        return nitrogen_xsect

    @staticmethod
    def argon_xsect(energy_ev):
        Elog = np.log(energy_ev)
        if energy_ev < 3202.9:
            X = -330.3509 + (267.7433 + 3.) * Elog - 78.90498 * Elog**2 \
                + 10.35983 * (Elog ** 3) - 0.5140201 * (Elog ** 4)
        elif energy_ev < 3202.9:
            X = -5.71870 + (8.85812 * Elog) + (-0.307357 * Elog * Elog) \
                + (0.00169351 * (Elog ** 3)) + (-0.0138134 * (Elog ** 4)) \
                + (0.00120451 * (Elog ** 5))
        else:
            X = 19.1905 + (2.74276 * Elog) + (-0.164603 *
                                              Elog * Elog) + (0.00165895*Elog**3)
        argon_xsect = np.exp(X)/(energy_ev**3)
        return argon_xsect

    @staticmethod
    def carbon_xsect(energy_ev):
        Elog = np.log(energy_ev)
        if energy_ev < 284.0:
            X = 8.74161 + (7.13348*Elog) + (-1.14604*Elog*Elog) + (0.0677044*Elog**3)
        else:
            X = 3.81334 + (8.93626*Elog) + (-1.06905*Elog*Elog) + (0.0422195*Elog**3)

        carbon_xsect = np.exp(X)/(energy_ev**3)
        return carbon_xsect

    @staticmethod
    def hydrog_xsect(energy_ev):
        Elog = np.log(energy_ev)
        X = 21.46941 + (3. - 2.060152) * Elog - 0.1492932 * (Elog**2) + 5.4634294e-3 * (Elog ** 3)
        hydrog_xsect = np.exp(X)/(energy_ev**3)
        return hydrog_xsect

    @staticmethod
    def helium_xsect(energy_ev):
        #data sets from Balucinska & Church 1992 fortran file
        c1 = [-2.953607e1, 7.083061e0, 8.678646e-1, -1.221932e0, 4.052997e-2, 1.317109e-1, -3.265795e-2, 2.500933e-3]
        c2 = [-2.465188e1, 4.354679e0, -3.553024e0, 5.573040e0, -5.872938e0, 3.720797e0, -1.226919e0, 1.576657e-1]
        Q = [2.81e0, 2.51e0, 2.45e0, 2.44e0]
        Nu = [1.610e0, 2.795e0, 3.817e0, 4.824e0]
        Gamma = [2.64061e-3, 6.20116e-4, 2.56061e-4, 1.320159e-4]
        lambDa = 12398.54e0 / energy_ev  # originally named lambda, shadows a python convention, hence capitalized D
        X = np.log10(lambDa)
        if lambDa > 503.97:
            helium_xsect = 0
        elif lambDa < 46:
            y = 0
            for i in range(8):
                y += c2[i]*(X**[i-1])
        else:
            y = 0
            for i in range(8):
                y += c1[i]*(X**[i-1])
            for i in range(4):
                # This section adapted from FANO(A, B, C, lambda) in B&C(1998)
                eps = 911.2671 / lambDa  # energy in Rydbergs
                epsi = 3 - 1 / (Nu[i] ** 2) + 1.807317
                X = 2 * (eps - epsi) / Gamma[i]
                y += np.log10(((X - Q[i]) ** 2) / (1 + X ** 2))
        sigma = 10**y
        helium_xsect = sigma * 6.022045e23 / 4.0026
        return helium_xsect

    @staticmethod
    def neon_xsect(energy_ev):
        Elog = np.log(energy_ev)
        if energy_ev < 867:
            X = -3.04041 + (13.0071 * Elog) + (-1.93205 * Elog**2) + (0.0977639 * Elog**3)
        else:
            X = 17.6007 + (3.29278 * Elog) + (-0.263065 * Elog**2) + (5.68290E-3 * Elog**3)
        neon_xsect = np.exp(X)/(energy_ev**3)
        return neon_xsect


def main():
    import matplotlib.pyplot as plt
    E_kev = np.linspace(1, 10, 1000)
    plt.figure(1)
    plt.plot(E_kev, BCM.get_total_xsect(E_kev, 1, 0, 0, 0, 0, 0), label = "N")
    plt.plot(E_kev, BCM.get_total_xsect(E_kev, 0, 1, 0, 0, 0, 0), label="O")
    plt.plot(E_kev, BCM.get_total_xsect(E_kev, 0, 0, 1, 0, 0, 0), label="Ar")
    plt.plot(E_kev, BCM.get_total_xsect(E_kev, 0, 0, 0, 1, 0, 0), label="C")
    plt.legend()

    plt.figure(2)
    E_kev = np.linspace(0.4, 1.0, 1000)
    plt.title("Photoelectric Cross Sections below 1 keV")
    plt.plot(E_kev, BCM.get_total_xsect(E_kev, 0, 1, 0, 0, 0, 0), label="Oxygen")
    plt.plot(E_kev, BCM.get_total_xsect(E_kev, 1, 0, 0, 0, 0, 0), label="Nitrogen")
    plt.plot(E_kev, BCM.get_total_xsect(E_kev, 0, 0, 1, 0, 0, 0), label="Argon")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Cross Section ($cm^2/g$)")
    plt.legend()
    plt.show()
    return 0

if __name__ == '__main__':
    main()
