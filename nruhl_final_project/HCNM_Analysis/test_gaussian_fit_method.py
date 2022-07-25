# Author:Seamus Flannery
# This script tests a gaussian fit method

import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import numpy as np
import sys
import os
import matplotlib.animation as manimation

# add working directory, str(Path(__file__).parents[1])
sys.path.append("/Users/nathanielruhl/Documents/HorizonCrossings-Summer22/nruhl_final_project/")
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project/")

# import local modules
from AnalyzeCrossing import AnalyzeCrossing

# Global parameters to be used in the analysis
cb_str = "Earth"  # P2: size of Earth, but a thinner atmosphere
hc_type = "rising"
N0 = 5000  # average number of unattenuated counts in data
E_kev = 1.5
H = 500  # km, orbital altitude
bin_size = 1.0
# range of transmittance in which to compare the curves
comp_range = [0.01, 0.99]
animate = True # animates sliding and analysis process for education and debugging
animate_res = 8 # fractional number of points to be rendered in animation
# Parameters involved in generating hc data
std = 0.05  # standard deviation of normally-distributed noise


def gaussian_fit(t_in, t_min_guess, range_list, chisq_list):
    t_min_guess_index = np.where(range_list >= t_min_guess)[0][0]
    start = t_min_guess_index - int(len(range_list) / 25)
    end = t_min_guess_index + int(len(range_list) / 25)
    low_b = range_list[t_min_guess_index - int(len(range_list) / 75)]
    high_b = range_list[t_min_guess_index + int(len(range_list) / 75)]
    chisq_min_guess = chisq_list[t_min_guess_index]
    low_k = 0.8 * chisq_min_guess
    high_k = 1.2 * chisq_min_guess
    low_bounds = [-500, low_b, 0, low_k]
    high_bounds = [0, high_b, 10, high_k]
    gaussian = lambda x, a, b, c, k: k + a * np.exp(-((x - b) ** 2) / (2 * c ** 2))
    popt, pcov = curve_fit(gaussian, range_list[start:end], chisq_list[start:end], bounds=(low_bounds, high_bounds))
    t_out = gaussian(t_in, *popt)
    return t_out


range = np.load("gaussian_test_data_dim.npy")[0]
chisq = np.load("gaussian_test_data_dim.npy")[1]

print(gaussian_fit(1999.73, 1999.72, range, chisq))

t_min = np.argmin(chisq)
print(range[t_min])
start = t_min - int(len(range)/25)
end = t_min + int(len(range)/25)
gaussian = lambda x, a, b, c, k : k + a * np.exp(-((x - b) ** 2) / (2 * c ** 2))
# gaussian = lambda x, a, b, k: k + a * np.exp(-1 * b * x ** 2)
# gaussian = lambda x, a, b, c: a*x**2+b*x+c
# FOR BRIGHT SOURCES, N0=5000:
# popt, pcov = curve_fit(gaussian, range[start:end], chisq[start:end], bounds=([-500, 1999.7, 0, 0], [0, 1999.75, 10, 2000]))
# FOR DIM SOURCES, N0=200:
popt, pcov = curve_fit(gaussian, range[start:end], chisq[start:end], bounds=([-500, 1999.7, 0, 0], [0, 1999.75, 10, 200]))
# popt, pcov = curve_fit(gaussian, range, chisq, bounds=([-500, 1999.68, 0, 0], [0, 1999.78, 10, 10000]))
print(popt)
plt.figure(1)
plt.plot(range, chisq, color="blue")
# plt.plot(range, gaussian(range, *popt), color="brown")
# plt.plot(range[start:end], gaussian(range[start:end], *popt), color="red")
smoothrange = np.linspace(range[0], range[len(range)-1], 1000)
plt.plot(smoothrange, gaussian(smoothrange, *popt), color="red")
plt.plot(range[start:end], chisq[start:end], color="green")
#plt.figure(2)
#plt.plot(range, chisq, color="blue")
#plt.plot(range, gaussian(range, *popt), color="red")
plt.show()


