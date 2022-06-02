# Author: Nathaniel Ruhl

# This script creates the chi-squared plot to show that Newton's method is well-behaved when solving for scale height

import numpy as np

from AnalyzeCrossing import AnalyzeCrossing
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Note that these are used
STD = 0.05  # Standard deviation of normal distribution from which noise is generated for generating the transmittance 'data'

# range in which noise is added to the model and in which rho0 is solved
COMP_RANGE = [0.01, 0.9]

# This function generates noisy data for a horizon crossing
def generate_crossing(SAT, plot_bool):
    # I copied these global variables here so I can call this function from another script
    STD = 0.05
    COMP_RANGE = [0.01, 0.9]
    time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)
    transmit_model = np.zeros_like(time_array)
    transmit_data = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_model[i] = np.exp(-SAT.tau_gauss(t, N=10))

    for i, t in enumerate(time_array):
        if (transmit_model[i] > COMP_RANGE[0]) & (transmit_model[i] < COMP_RANGE[1]):
            transmit_data[i] = transmit_model[i]*np.random.normal(1, STD)
        else:
            transmit_data[i] = transmit_model[i]

    if plot_bool is True:
        plt.figure()
        plt.title("Simulated horizon crossing data")
        plt.ylabel("Transmittance")
        plt.xlabel("Time (s)")
        plt.plot(time_array, transmit_data, ".")

    return transmit_data

# Function to calculate a transmittance array for a given SAT/orbital parameters
def calc_transmit(SAT, time_array):
    transmit_array = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_array[i] = np.exp(-SAT.tau_gauss(t, N=10))
    return transmit_array

SAT = AnalyzeCrossing(cb="Earth", H=420, E_kev=4.0)
transmit_data = generate_crossing(SAT, plot_bool=False)
comp_indices = np.where((transmit_data > COMP_RANGE[0]) & (
    transmit_data < COMP_RANGE[1]))[0]  # indices for comparison

time_array = np.arange(0, SAT.time_final + 1, 1, dtype=float)

# list of scaling factors for scale_height
alpha_list = np.arange(0.9, 1.1, 0.01)
chisq_list = []
for alpha in alpha_list:
    SAT.scale_height = alpha*SAT.scale_height
    transmit_array = calc_transmit(SAT, time_array)

    # Calculate chisq
    
    # Creat interpolating function for model
    transmit_model_vs_time = interp1d(x=time_array, y=transmit_array, kind='cubic')
    chisq = np.sum(((transmit_data[comp_indices]-transmit_model_vs_time(time_array[comp_indices]))/STD)**2)
    chisq_list.append(chisq)

plt.figure()
plt.plot(alpha_list, chisq_list, '.')
plt.show()



