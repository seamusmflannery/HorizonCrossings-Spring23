# Author: Nathaniel Ruhl

# Driver for HCNM process for a NICER observation

# import local modules for HCNM Analysis
from Modules.OrbitModel import OrbitModel
from Modules.LocateR0hc import LocateR0hc
from Modules.ReadEVT import ReadEVT
from Modules.NormalizeSpectrumNICER import NormalizeSpectrumNICER
from Modules.TransmitModel import TransmitModel
from Modules.CurveComparison import CurveComparison

# import observation dictionaries
from ObservationDictionaries.v4641NICER import v4641NICER

# Process of HCNM to be used in every script


def main():
    obs_dict = v4641NICER
    e_band = [1.0, 2.0]  # keV
    bin_size = 1.0   # sec

    # 1) Define orbit model
    r_array, t_array = OrbitModel.define_orbit_model(obs_dict, "mkf", time_step=0.01)

    # 2) LocateR0hc (must know hc_type here (from obs_dict or step 1))
    r0_obj = LocateR0hc(obs_dict, r_array, t_array)
    t0_model_index, lat_gp, lon_gp = r0_obj.return_orbit_data()
    del r0_obj

    orbit_derived_inputs = (r_array, t_array, t0_model_index, lat_gp, lon_gp)

    # 3) Bin the data. If repeating for a different energy band, can start analysis at step 3.
    # Can also identify hc_type in this step.
    evt_obj = ReadEVT(obs_dict)
    rate_data, time_data, unattenuated_rate = evt_obj.return_crossing_data(e_band, bin_size)

    spec_obj = NormalizeSpectrumNICER(evt_obj, e_band)
    normalized_amplitudes, bin_centers = spec_obj.return_spectrum_data()
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
    t0_e, dt_e = comp_obj.t0_e, comp_obj.dt_e
    del comp_obj

    print(f"Time at the position r0_hc = {r_array[t0_model_index]}:")
    print(f"Crossing: t0_e = {t0_e} +/- {dt_e} sec")
    print(f"Input Orbit Model: t0 = {t_array[t0_model_index]} sec")

    #7) weighted_mean_HC.py


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))