# Author: Nathaniel Ruhl & Seamus Flannery
# This script makes and animates the curve comparison. This script should be used whenever working to improve the algorithm

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import numpy as np
import sys

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
animate = True  # animates sliding and analysis process for education and debugging
animate_res = 8  # fractional number of points to be rendered in animation
# Parameters involved in generating hc data
std = 0.05  # standard deviation of normally-distributed noise
np.random.seed(3)


# This function generates the horizon crossing times and transmittance arrays for both model and data
def generate_crossings(sat, hc_type):
    # Define time_model differently if setting or rising HC
    time_model = np.arange(0, sat.time_final, bin_size)
    if hc_type == "setting":
        time_model = np.flip(time_model)

    transmit_model = sat.calc_transmit_curve(time_model)

    transmit_data = transmit_model.copy()

    if hc_type == "rising":
        transmit_data = np.insert(transmit_data, 0, np.zeros(int(100 / bin_size)))
        transmit_data = np.append(transmit_data, np.ones(int(100 / bin_size)))
        noise_range = np.where(transmit_data > 0.05)[0]
        transmit_data[noise_range] += np.random.normal(0, std, len(transmit_data[noise_range]))
    elif hc_type == "setting":
        transmit_data = np.insert(transmit_data, 0, np.ones(int(100 / bin_size)))
        transmit_data = np.append(transmit_data, np.zeros(int(100 / bin_size)))
        noise_range = np.where(transmit_data > 0.05)[0]
        transmit_data[noise_range] += np.random.normal(0, std, len(transmit_data[noise_range]))

    time_data = np.arange(2000 - 100, 2000 + sat.time_final + 100, bin_size)
    return time_model, transmit_model, time_data, transmit_data


# This class executes the CurveComparison
def gaussian(x, a, b, c, k):
    return k + a * np.exp(-((x - b) ** 2) / (2 * c ** 2))


class CurveComparison:
    def __init__(self, sat, hc_type, N0, interp_type="cubic"):
        self.interp_type = interp_type
        self.sat = sat
        self.hc_type = hc_type  # want this to be autonomous
        self.N0 = N0
        self.time_model, self.transmit_model, self.time_data, self.transmit_data = generate_crossings(self.sat,
                                                                                                      self.hc_type)
        # First step to identify t0
        self.t0_1 = self.locate_t0_step1()
        self.t0_2, self.t0_guess_list, self.chisq_list = self.locate_t0_step2()
        self.dt_2, self.chisq_upper_t2, self.chisq_lower_t2 = self.analyze_chisq()
        self.t0_3, self.t0_range_list, self.chisq_list_fine = self.locate_t0_fine()
        self.dt_e, self.t0_e = self.analyze_chisq_fine()

    # This function is used to identify the model time (from t0)
    def time_vs_transmit_model(self, transmit_approx_50):
        frange = np.where((self.transmit_model > 0.01) & (self.transmit_model < 0.99))[
            0]  # indices for which transmit vs time is 1-1
        time_vs_transmit = interp1d(x=self.transmit_model[frange], y=self.time_model[frange], kind=self.interp_type)
        return time_vs_transmit(transmit_approx_50)

    # first guess of t0 by comparing 50% points
    def locate_t0_step1(self):
        if self.hc_type == "setting":
            index50_approx_data = np.where(self.transmit_data < 0.51)[0][0]
        elif self.hc_type == "rising":
            index50_approx_data = np.where(self.transmit_data > 0.49)[0][0]
        transmit50_approx_data = self.transmit_data[index50_approx_data]
        t50_approx_data = self.time_data[index50_approx_data]

        dt50_approx_model = self.time_vs_transmit_model(transmit50_approx_data)

        if self.hc_type == "rising":
            t0_1 = t50_approx_data - dt50_approx_model
        elif self.hc_type == "setting":
            t0_1 = t50_approx_data + dt50_approx_model
        return t0_1

    # This function slides the model past the data at time intervals of desired_precision and calculates chi_sq
    def locate_t0_step2(self):
        desired_precision = 0.01
        t0_1_index = np.where(self.time_data >= self.t0_1)[0][0]

        # Define the data points in the full crossing time range
        if self.hc_type == "setting":
            time_crossing_data = self.time_data[t0_1_index - len(self.time_model):t0_1_index]
            rate_data = self.N0 * self.transmit_data[t0_1_index - len(self.time_model):t0_1_index]
            transmit_data = self.transmit_data[t0_1_index -
                                               len(self.time_model):t0_1_index]
        elif self.hc_type == "rising":
            time_crossing_data = self.time_data[t0_1_index:t0_1_index + len(self.time_model)]
            rate_data = self.N0 * self.transmit_data[t0_1_index:t0_1_index + len(self.time_model)]
            transmit_data = self.transmit_data[t0_1_index -
                                               len(self.time_model):t0_1_index]
        # set ranges for testing
        weight_range = np.where((self.transmit_model >= comp_range[0]) & (self.transmit_model <= comp_range[1]))[0]
        side_range = 1 / 25 * len(weight_range)
        t_start_list = np.arange(self.t0_1 - side_range,
                                 self.t0_1 + side_range,
                                 desired_precision)
        if animate:
            plt.figure(1)
            fig, ax = plt.subplots(1, 2)
            fig.suptitle("Horizon Crossing Curve Comparison")
            fig.supxlabel("Time (sec)")
        chisq_list = np.zeros(len(t_start_list))
        for indx, t0_guess in enumerate(t_start_list):
            if animate:
                ax[0].clear()
            # define interpolating function and array for the model
            if self.hc_type == "rising":
                time_crossing_model = np.arange(t0_guess, t0_guess + self.sat.time_final, bin_size)
            elif self.hc_type == "setting":
                time_crossing_model = np.flip(np.arange(t0_guess, t0_guess - self.sat.time_final, -bin_size))
            else:
                continue

            # Note that however this interpolation is done, the model and data times need to be in the same order
            model_rate_vs_time = interp1d(time_crossing_model, self.N0 * self.transmit_model, kind=self.interp_type)
            model_rate_interp = model_rate_vs_time(time_crossing_data[weight_range])
            # List of model values at times where data points are

            if any(model_rate_interp <= 0):
                print(self.interp_type + " spline went negative")

            # Chi-squared test in weight_range of full curve

            chisq = np.sum(
                (rate_data[weight_range] - model_rate_interp) ** 2 / model_rate_interp)
            chisq_list[indx] = chisq

            if animate and indx % animate_res == 0:
                ax[0].plot(time_crossing_data[weight_range], model_rate_interp, 'b-')
                ax[0].plot(time_crossing_data[weight_range], rate_data[weight_range], 'kx')
                ax[1].plot(t0_guess, chisq, 'r.')
                ax[0].set_xlim([min(time_crossing_data[weight_range] - 1), max(time_crossing_data[weight_range]) + 1])
                ax[0].set_ylim([0, max(rate_data[weight_range]) + 10])
                ax[0].set_ylabel("Counts/sec")
                ax[1].set_ylabel(r"$\chi^2$")
                plt.pause(0.01)

        t0_2 = t_start_list[np.argmin(chisq_list)]
        # np.save("gaussian_test_data_dim.npy", [t_start_list, chisq_list])
        return t0_2, t_start_list, chisq_list

    # Methods to calculate the chisq+1 error
    def analyze_chisq(self):
        chisq_func = lambda t: self.chisq_vs_time(t, self.t0_guess_list, self.chisq_list)
        t_min_chisq = float(least_squares(chisq_func, x0=self.t0_2).x)
        print("min chisq step 2:" + str(t_min_chisq))
        plus_1 = lambda t: chisq_func(t) - (chisq_func(t_min_chisq) + 1)
        upper_t0 = root_scalar(plus_1, method="secant", x0=t_min_chisq + 0.1, x1=t_min_chisq + 0.3).root
        lower_t0 = root_scalar(plus_1, method="secant", x0=t_min_chisq - 0.1, x1=t_min_chisq - 0.3).root
        if animate:
            print("+1 Root found: " + str(upper_t0))
            print("+1 Root found: " + str(lower_t0))
        # return the larger of the two
        if self.chisq_vs_time(upper_t0, self.t0_guess_list, self.chisq_list) > self.chisq_vs_time(lower_t0,
                                                                                                  self.t0_guess_list,
                                                                                                  self.chisq_list):
            chisq_error = upper_t0 - self.t0_2
        else:
            chisq_error = abs(lower_t0 - self.t0_2)

        plus_2 = lambda t: chisq_func(t) - (chisq_func(t_min_chisq) + 2)
        upper_t0_2 = root_scalar(plus_2, method="secant", x0=t_min_chisq + 0.1, x1=t_min_chisq + 0.3).root
        lower_t0_2 = root_scalar(plus_2, method="secant", x0=t_min_chisq - 0.1, x1=t_min_chisq - 0.3).root
        if animate:
            print("+2 Root found: " + str(upper_t0_2))
            print("+2 root y: " + str(chisq_func(upper_t0_2)))
            print("+2 Root found: " + str(lower_t0_2))
            print("+2 root y: " + str(chisq_func(lower_t0_2)))

        return chisq_error, upper_t0_2, lower_t0_2

    def chisq_vs_time(self, t0, range_list, chisq_list):
        func = interp1d(range_list, chisq_list, kind='cubic', fill_value="extrapolate")
        return func(t0)

    def gaussian_fit(self, t_min_guess, range_list, chisq_list):
        t_min_guess_index = np.where(range_list >= t_min_guess)[0][0]
        start = t_min_guess_index - int(len(range_list) / 25)
        end = t_min_guess_index + int(len(range_list) / 25)
        low_b = range_list[t_min_guess_index - int(len(range_list) / 75)]
        upper_b = range_list[t_min_guess_index + int(len(range_list) / 75)]
        chisq_min_guess = chisq_list[t_min_guess_index]
        low_k = 0.8 * chisq_min_guess
        upper_k = 1.5 * chisq_min_guess
        low_bounds = [-500, low_b, 0, low_k]
        upper_bounds = [0, upper_b, 10, upper_k]
        popt, pcov = curve_fit(gaussian, range_list[start:end], chisq_list[start:end], bounds=(low_bounds, upper_bounds))
        if animate:
            print("a: " + str(popt[0]) + ", b: " + str(popt[1]) + ", c: " + str(popt[2]) + ", k: " + str(popt[3]))
            print("lower fit bounds: ")
            print(low_bounds)
            print("upper fit bounds")
            print(upper_bounds)
            plt.figure(4)
            plt.plot(range_list[start:end], chisq_list[start:end], color="brown")
        return popt

    def locate_t0_fine(self):
        zoom_factor = 10
        resolution_factor = 10
        desired_precision = 0.01 / resolution_factor
        t0_2_index = np.where(self.time_data >= self.t0_2)[0][0]

        # Define the data points in the full crossing time range
        if self.hc_type == "setting":
            time_crossing_data = self.time_data[t0_2_index - len(self.time_model):t0_2_index]
            rate_data = self.N0 * self.transmit_data[t0_2_index - len(self.time_model):t0_2_index]
            transmit_data = self.transmit_data[t0_2_index -
                                               len(self.time_model):t0_2_index]
        elif self.hc_type == "rising":
            time_crossing_data = self.time_data[t0_2_index:t0_2_index + len(self.time_model)]
            rate_data = self.N0 * self.transmit_data[t0_2_index:t0_2_index + len(self.time_model)]
            transmit_data = self.transmit_data[t0_2_index -
                                               len(self.time_model):t0_2_index]
        # set up ranges for fine testing
        weight_range = np.where((self.transmit_model >= comp_range[0]) & (self.transmit_model <= comp_range[1]))[0]
        # side_range = 1 / (25 * zoom_factor) * len(weight_range)
        t0_range_list = np.arange(self.chisq_lower_t2,
                                  self.chisq_upper_t2,
                                  desired_precision)
        if animate:
            plt.figure(2)
            fig, ax = plt.subplots(1, 2)
            fig.suptitle("Horizon Crossing Curve Comparison Fine Zoom")
            fig.supxlabel("Time (sec)")
        chisq_list_fine = np.zeros(len(t0_range_list))
        for indx, t0_guess in enumerate(t0_range_list):
            if animate:
                ax[0].clear()
            # define interpolating function and array for the model
            if self.hc_type == "rising":
                time_crossing_model = np.arange(t0_guess, t0_guess + self.sat.time_final, bin_size)
            elif self.hc_type == "setting":
                time_crossing_model = np.flip(np.arange(t0_guess, t0_guess - self.sat.time_final, -bin_size))

            # Note that however this interpolation is done, the model and data times need to be in the same order
            model_rate_vs_time = interp1d(time_crossing_model, self.N0 * self.transmit_model, kind=self.interp_type)
            model_rate_interp = model_rate_vs_time(time_crossing_data[weight_range])
            # List of model values at times where data points are

            if any(model_rate_interp <= 0):
                print(self.interp_type + " spline went negative")

            # Chi-squared test in weight_range of full curve
            chisq = np.sum(
                (rate_data[weight_range] - model_rate_interp) ** 2 / model_rate_interp)
            chisq_list_fine[indx] = chisq

            if animate and indx % animate_res == 0:
                ax[0].plot(time_crossing_data[weight_range], model_rate_interp, 'b-')
                ax[0].plot(time_crossing_data[weight_range], rate_data[weight_range], 'kx')
                ax[1].plot(t0_guess, chisq, 'r.')
                ax[0].set_xlim([min(time_crossing_data[weight_range] - 1), max(time_crossing_data[weight_range]) + 1])
                ax[0].set_ylim([0, max(rate_data[weight_range]) + 10])
                ax[0].set_ylabel("Counts/sec")
                ax[1].set_ylabel(r"$\chi^2$")
                plt.pause(0.01)

        t0_3 = t0_range_list[np.argmin(chisq_list_fine)]
        if animate:
            plt.figure(3)
            interp_range = np.linspace(min(t0_range_list), max(t0_range_list), 1000)
            interp_data = self.chisq_vs_time(np.linspace(min(t0_range_list), max(t0_range_list), 1000), t0_range_list,
                                             chisq_list_fine)
            ax[1].plot(interp_range, interp_data, color="red")
            plt.show()
        return t0_3, t0_range_list, chisq_list_fine

    def analyze_chisq_fine(self):
        popt = self.gaussian_fit(self.t0_3, self.t0_range_list, self.chisq_list_fine)
        chisq_func = lambda t: gaussian(t, *popt)
        t_min_chisq = float(least_squares(chisq_func, x0=self.t0_3).x)
        print("min chisq step 3: " + str(t_min_chisq))
        plus_1 = lambda t: chisq_func(t) - (chisq_func(t_min_chisq) + 1)
        upper_t0 = root_scalar(plus_1, method="secant", x0=t_min_chisq + 0.1, x1=t_min_chisq + 0.3).root
        lower_t0 = root_scalar(plus_1, method="secant", x0=t_min_chisq - 0.1, x1=t_min_chisq - 0.3).root
        if animate:
            smoothrange = np.linspace(self.t0_range_list[0], self.t0_range_list[len(self.t0_range_list) - 1], 1000)
            plt.figure(4)
            plus_one_line = np.full(1000, chisq_func(t_min_chisq)+1)
            plt.title("gaussian fit check")
            plt.plot(smoothrange, gaussian(smoothrange, *popt), color="red")
            plt.plot(self.t0_range_list, self.chisq_list_fine, "bo", markersize="2")
            plt.plot(smoothrange, plus_one_line, color="green")
            plt.show()
            print("Final +1 Root found: " + str(upper_t0))
            print("Final +1 Root found: " + str(lower_t0))
        # return the larger of the two
        if chisq_func(upper_t0) > chisq_func(lower_t0):
            chisq_error = upper_t0 - self.t0_3
        else:
            chisq_error = abs(lower_t0 - self.t0_3)

        return chisq_error, t_min_chisq


if __name__ == "__main__":
    import time
    start_time = time.time()
    sat = AnalyzeCrossing(cb_str, H, E_kev)
    comp_obj = CurveComparison(sat, hc_type, N0=N0)
    print("dt_e: " + str(comp_obj.dt_e))
    print("t0_e: " + str(comp_obj.t0_e))
    print("--- %s seconds ---" % (time.time() - start_time))
