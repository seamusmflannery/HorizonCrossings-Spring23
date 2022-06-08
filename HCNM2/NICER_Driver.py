# Author: Nathaniel Ruhl

# Driver for HCNM process

import numpy as np

# import local modules
from OrbitModel import OrbitModel
from LocateR0hc import LocateR0hc
from ReadEVT import ReadEVT
from NormalizeSpectrumNICER import NormalizeSpectrumNICER
from TransmitModel import TransmitModel
from CurveComparison import CurveComparison

# import observation dictionaries
from ObservationDictionaries.v4641 import v4641

# Process of HCNM to be used in every script


def main():
    obs_dict = v4641
    e_band = [2.0, 3.0] # keV
    bin_size = 1.0  # sec

    #1) Bin the data. Can identify hc_type in this step
    evt_obj = ReadEVT(obs_dict)
    rate_data, time_data, unattenuated_rate = evt_obj.return_crossing_data(e_band, bin_size)
    print(unattenuated_rate)
    # 2) Define orbit model
    r_array, t_array = OrbitModel.define_orbit_model(obs_dict, "mkf", time_step=0.01)

    # 3) LocateRohc (must know hc_type here (from obs_dict or step 1))
    r0_obj = LocateR0hc(obs_dict, r_array, t_array)
    t0_model_index, lat_gp, lon_gp = r0_obj.t0_model_index, r0_obj.lat_gp, r0_obj.lon_gp
    del r0_obj
    print(f'r0_hc = {r_array[t0_model_index]}')

    orbit_derived_inputs = (r_array, t_array, t0_model_index, lat_gp, lon_gp)

    # 4)
    spec_obj = NormalizeSpectrumNICER(evt_obj, e_band)
    normalized_amplitudes, bin_centers = spec_obj.normalized_amplitudes, spec_obj.bin_centers
    del evt_obj
    del spec_obj

    eband_derived_inputs = (e_band, bin_size, normalized_amplitudes, bin_centers)

    # 5) RXTE Analysis starts here. Not currently working
    model_obj = TransmitModel(obs_dict, orbit_derived_inputs, eband_derived_inputs)
    transmit_model, time_crossing_model = model_obj.calculate_transmit_model()

    # GENERATE DATA IF NEEDED
    # time_crossing_model = np.arange(0, 300, bin_size)
    # transmit_model = np.exp(-np.exp(-time_crossing_model+20))

    # 6)
    model_and_data_tuple = (time_crossing_model, transmit_model, time_data, rate_data, unattenuated_rate)

    comp_obj = CurveComparison(obs_dict, model_and_data_tuple)
    print(comp_obj.t0_1)
    t0_e, dt_e = comp_obj.t0_e, comp_obj.dt_e
    del comp_obj

    print(f"Time at the position r0_hc = {r_array[t0_model_index]}:")
    print(f"Crossing: t0_e = {t0_e} +/- {dt_e} sec")
    print(f"Input Orbit Model: t0 = {t_array[t0_model_index]} sec")

    #7 weighted_mean_HC.py


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
