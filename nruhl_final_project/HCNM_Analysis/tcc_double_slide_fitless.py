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
N0 = 5378  # average number of unattenuated counts in data
E_kev = 1.5
H = 2200  # km, orbital altitude
bin_size = 1.0
# range of transmittance in which to compare the curves
comp_range = [0.01, 0.99]
animate = True  # animates sliding and analysis process for education and debugging
plot_toggle = True
animate_res = 8  # fractional number of points to be rendered in animation
# Parameters involved in generating hc data
std = 0.05  # standard deviation of normally-distributed noise
np.random.seed(3)
initial_precision = 0.01


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
class CurveComparison:
    def __init__(self, sat, hc_type, N0, interp_type="linear"):
        self.interp_type = interp_type
        self.sat = sat
        self.hc_type = hc_type  # want this to be autonomous
        self.N0 = N0
        self.time_model, self.transmit_model, self.time_data, self.transmit_data = generate_crossings(self.sat,
                                                                                                      self.hc_type)
        # First step to identify t0
        self.t0_1 = self.locate_t0_step1()
        self.t0_1_index = np.where(self.time_data >= self.t0_1)[0][0]
        # Locate t0_2 and narrow range, calculate chisq for narrowed range
        self.t0_2, self.t0_guess_list, self.chisq_list = self.locate_t0_step2()
        # centered around t0_2, narrows range further to +3 chisq,
        self.refined_search_range = self.refine_search_range()
        # Repeat step 2 on narrower range, calculate chisq for narrower range, higher resolution
        self.t0_e, self.t0_range_list, self.chisq_list_fine, self.dt_e = self.locate_t0_fine()

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
        # Define the data points in the full crossing time range
        if self.hc_type == "setting":
            time_crossing_data = self.time_data[self.t0_1_index - len(self.time_model):self.t0_1_index]
            rate_data = self.N0 * self.transmit_data[self.t0_1_index - len(self.time_model):self.t0_1_index]
            transmit_data = self.transmit_data[self.t0_1_index -
                                               len(self.time_model):self.t0_1_index]
        elif self.hc_type == "rising":
            time_crossing_data = self.time_data[self.t0_1_index:self.t0_1_index + len(self.time_model)]
            rate_data = self.N0 * self.transmit_data[self.t0_1_index:self.t0_1_index + len(self.time_model)]
            transmit_data = self.transmit_data[self.t0_1_index -
                                               len(self.time_model):self.t0_1_index]
        # set ranges for testing
        weight_range = np.where((self.transmit_model >= comp_range[0]) & (self.transmit_model <= comp_range[1]))[0]
        side_range = 1 / 10 * len(weight_range)
        if self.sat.cb != "Earth":
            t_start_list = np.arange(self.t0_1 - 2, self.t0_1 + 2, initial_precision)
        else:
            t_start_list = np.arange(self.t0_1 - side_range, self.t0_1 + side_range, initial_precision)
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
        return t0_2, t_start_list, chisq_list

    # Return the range of +3 chisq around the lowest chisq value for a later, higher granularity search
    def refine_search_range(self):
        minchisq = min(self.chisq_list)
        min_index = np.argmin(self.chisq_list)
        left_half = self.chisq_list[:min_index]  # data to the left of the min
        right_half = self.chisq_list[min_index:]  # data to the right of the min
        # finds left bound on +3 chisq
        newmin_index = max(np.where(left_half > minchisq + 3)[0])
        # finds right bound on +3 chisq - convoluted but takes the right half of the data, finds where its bigger than
        # plus three from minchisq, finds the exact value of the chisq there in order to be able to find that value in
        # the full chisq list (we don't want the right half array indexing)
        newmax_index = np.where(self.chisq_list == right_half[min(np.where(right_half > minchisq + 3)[0])])[0][0]
        refined_search_range = np.array(self.t0_guess_list[newmin_index:newmax_index])
        if plot_toggle:
            smoothrange = np.linspace(self.t0_guess_list[0], self.t0_guess_list[len(self.t0_guess_list) - 1], 1000)
            plt.figure(3)
            plus_one_line = np.full(1000, minchisq + 1)
            plus_three_line = np.full(1000, minchisq + 3)
            plt.plot(self.t0_guess_list, self.chisq_list, "bo", markersize="2", label="chisq data")
            plt.plot(smoothrange, plus_one_line, color="green", label="chisq minimum +1")
            plt.plot(smoothrange, plus_three_line, color="green", label="chisq minimum +3")
            plt.legend()
            plt.show()
        return refined_search_range

    def locate_t0_fine(self):
        refinement_factor = 25
        global initial_precision
        refined_precision = initial_precision / refinement_factor

        # Define the data points in the full crossing time range
        if self.hc_type == "setting":
            time_crossing_data = self.time_data[self.t0_1_index - len(self.time_model):self.t0_1_index]
            rate_data = self.N0 * self.transmit_data[self.t0_1_index - len(self.time_model):self.t0_1_index]
            transmit_data = self.transmit_data[self.t0_1_index -
                                               len(self.time_model):self.t0_1_index]
        elif self.hc_type == "rising":
            time_crossing_data = self.time_data[self.t0_1_index:self.t0_1_index + len(self.time_model)]
            rate_data = self.N0 * self.transmit_data[self.t0_1_index:self.t0_1_index + len(self.time_model)]
            transmit_data = self.transmit_data[self.t0_1_index -
                                               len(self.time_model):self.t0_1_index]
        # set up ranges for fine testing
        weight_range = np.where((self.transmit_model >= comp_range[0]) & (self.transmit_model <= comp_range[1]))[0]
        t0_range_list = np.arange(min(self.refined_search_range),
                                  max(self.refined_search_range),
                                  refined_precision)
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

        t0_4 = t0_range_list[np.argmin(chisq_list_fine)]
        t0_4_index = np.where(t0_range_list == t0_4)[0][0]
        minchisq = min(chisq_list_fine)
        left_half = chisq_list_fine[:t0_4_index]  # data to the left of the min
        right_half = chisq_list_fine[t0_4_index:]  # data to the right of the min
        # finds left bound on +1 chisq
        lower_chisq_plus_one_index = max(np.where(left_half > minchisq + 1)[0])
        # finds right bound on +1 chisq - convoluted but takes the right half of the data, finds where it's bigger than
        # plus one from minchisq, finds the exact value of the chisq there in order to be able to find that value in
        # the full chisq list (we don't want the right half array indexing)
        chisq_plus_one = right_half[min(np.where(right_half > minchisq + 1)[0])]
        upper_chisq_plus_one_index = np.where(chisq_list_fine == right_half[min(np.where(right_half > minchisq + 1)[0])])[0][0]
        if abs(t0_range_list[lower_chisq_plus_one_index] - t0_4) > abs(t0_range_list[upper_chisq_plus_one_index] - t0_4):
            final_uncertainty = abs(t0_range_list[lower_chisq_plus_one_index] - t0_4)
        else:
            final_uncertainty = abs(t0_range_list[upper_chisq_plus_one_index] - t0_4)
        if plot_toggle:
            smoothrange = np.linspace(t0_range_list[0], t0_range_list[len(t0_range_list) - 1], 1000)
            plt.figure(4)
            plus_one_line = np.full(1000, minchisq + 1)
            plus_two_line = np.full(1000, minchisq + 2)
            plt.plot(t0_range_list, chisq_list_fine, "bo", markersize="2", label="chisq data")
            plt.plot(smoothrange, plus_one_line, color="green", label="chisq minimum +1, dt_e = "+str(final_uncertainty))
            plt.plot(smoothrange, plus_two_line, color="green", label="chisq minimum +2")
            plt.legend()
            plt.show()
        return t0_4, t0_range_list, chisq_list_fine, final_uncertainty


if __name__ == "__main__":
    import time

    start_time = time.time()
    sat = AnalyzeCrossing(cb_str, H, E_kev)
    comp_obj = CurveComparison(sat, hc_type, N0=N0)
    print("dt_e: " + str(comp_obj.dt_e))
    print("t0_e: " + str(comp_obj.t0_e))
    print("--- %s seconds ---" % (time.time() - start_time))
