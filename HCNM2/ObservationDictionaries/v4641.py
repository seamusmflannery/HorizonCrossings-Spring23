import numpy as np

# import local modules
import tools as tools

v4641 = {   # BASIC OBSERVATION INFO
            "detector": "NICER",
            "source_name": "V4641 Sgr.",
            "source_RA": 274.839,  # deg
            "source_DEC": -25.407, # deg
            "starECI": tools.celestial_to_geocentric(np.deg2rad(274.839), np.deg2rad(-25.407)),
            "hc_type": "rising",

            # 2 fields below will need to be made autonomous eventually
            # Used in LocateR0hc.py
            "h_unit": np.array([-0.72023941, -0.30720814, 0.62199545]),
            "R_orbit": 6796.83588, # km

            # USER-ENTERED INFORMATION
            "crossing_time_range": np.array([300 + 1.92224e8, 760 + 1.92224e8]),  # seconds in MET
            "spectrum_time_range": np.array([550 + 1.92224e8, 690 + 1.92224e8]),  # only for NICER
            "t50_datetime": 1,  # approximate
            "f107": 100,
            "Ap": 7.0,

            # PATHS TO DATA FILES (from cwd, HCNM2/)
            "spectrum_path": None,

            "evt_path": "Data/NICER/2-3-20-v4641/NICER_events.evt", # NICER events file
            "mkf_path": "Data/NICER/2-3-20-v4641/ISS_orbit.mkf", # NICER orbital solution

            "lc_path": None, # RXTE binned data
            "orb_path": None, # RXTE orbital solution

            "aster_path": None # ASTER Labs orbital solution
}
