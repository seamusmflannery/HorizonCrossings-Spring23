# Author: Nathaniel Ruhl
# This file explores the dimensionless parameters identified with the Buckingham Pi theorem

import numpy as np
import matplotlib.pyplot as plt

from AnalyzeCrossing import AnalyzeCrossing

def main(SAT):
    time_array = np.arange(0, SAT.time_final+1, 1)
    transmit = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit[i] = np.exp(-SAT.tau_gauss(t, 50))

    # Dimensionless parameters
    h_star = SAT.tan_alt(time_array)/SAT.scale_height
    lambda_star = SAT.rho0*SAT.sigma*SAT.tan_alt(time_array)
    Rstar = SAT.rho0*SAT.sigma*SAT.scale_height
    plt.figure(1)
    plt.plot(h_star, lambda_star, label=fr"{SAT.cb} satellite, $R^*=${Rstar}")
    plt.ylabel(r"$\lambda^*$")
    plt.xlabel(r"$h^*$")
    plt.legend()
    plt.figure(2)
    plt.plot(h_star, transmit, label=f"{SAT.cb} satellite")
    plt.ylabel("Transmittance")
    plt.xlabel("$h^*$")
    plt.legend()
    plt.figure(3)
    plt.plot(lambda_star, transmit, label=f"{SAT.cb} satellite")
    plt.legend()
    plt.figure(4)
    plt.plot(time_array, transmit, label=f"{SAT.cb} satellite")
    plt.legend()
    return 0

if __name__ == '__main__':
    ES = AnalyzeCrossing(cb="Earth", H=420)
    MS = AnalyzeCrossing(cb="Mars", H=420)
    P1S = AnalyzeCrossing(cb="P1", H=420)
    # ES2 = AnalyzeCrossing(cb="Earth", H=500)
    main(ES)
    main(P1S)
    plt.show()

