# Author: Nathaniel Ruhl

# import local modules
from OrbitModel import OrbitModel
from LocateR0hc import LocateR0hc

# import observation dictionaries
from ObservationDictionaries.v4641 import v4641

# Process of HCNM to be used in every script

# Splitting the process up means that we can loop over eband_derived_inputs and store orbit_derived_inputs to analyze a single hc without re-calculating r0_hc each time

obs_dict = v4641
e_band = [1.0, 2.0] # keV
bin_size = 1.0 # sec

####### ORBIT-DERIVED INPUTS ########
# Future directions:
# class variables for keplerian elements
# class OrbitDerivedParameters(obs_dict, orbit_model_type)
#1)
r_array, t_array = OrbitModel.define_orbit_model(obs_dict, "mkf", time_step=0.01)

#2)
r0_obj = LocateR0hc(obs_dict, r_array, t_array)
t0_model_index, lat_gp, lon_gp = r0_obj.return_r0_hc_params()
print(r_array[t0_model_index])
# Remember that lat_gp, lon_gp can be defined from g_unit

orbit_derived_inputs = (r_array, t_array, t0_model_index)

'''

######## EBAND-DERIVED INPUTS #######
# class EbandDerivedParameters(obs_dict, e_band, bin_size, dE_eff)
# def return_eband_derived_outputs()
#3)
normalized_spectrum = get_normalized_spectrum(obs_dict, e_band)

# 4) NICER or RXTE ("evt" (unbinned) or "lc" (binned))? It will be done differently in each case
rate_data, time_data, unattenuated_rate = bin_data(obs_dict, e_band, bin_size)

# 5) class variables like mix_NOAr, pymsis_args, etc..
transmit_model, transmit_times = TransmitOrbitModel(obs_dict, e_band, dE_eff)
rate_model = unattenuated_rate*transmit_model

eband_derived_inputs = (eband, bin_size, rate_data, time_data, rate_model, time_crossing_model)

# 6)
t0_e, dt_e = CurveComparison(obs_dict, orbit_derived_inputs, eband_derived_inputs)

# 7) weighted_mean_HC.py

'''
