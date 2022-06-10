import numpy as np
from pathlib import Path
cwd = str(Path(__file__).parents[1])  # HCNM2/ is cwd

from Modules import tools as tools

crabNICER = {   # BASIC OBSERVATION INFO
            "detector": "NICER",
            "source_name": "Crab Nebula",
            "obsID": 4522010103,
            "source_RA": 83.63317,  # deg
            "source_DEC": 22.01453,  # deg
            "starECI": tools.celestial_to_geocentric(np.deg2rad(83.63317), np.deg2rad(22.01453)),
            "hc_type": "rising",

            # USER-ENTERED INFORMATION

            "crossing_time_range": np.array([240165000, 240165500]),  # seconds in MET
            "spectrum_time_range": np.array([240165300, 240165400]),  # only for NICER
            "f107": 75.2,
            "ap": 2,

            # 2 FIELDS FOR USER TO DETERMINE
            # Used in LocateR0hc.py
            "h_unit": np.array([0.7279706, -0.29417425, 0.6192902]),
            "R_orbit": 6801,  # km (approximate)

            # PATHS TO DATA FILES (from cwd, HCNM2/)

            "evt_path": cwd + "/Data/NICER/obs4522010103/ni4522010103_0mpu7_cl.evt",  # NICER events file
            "mkf_path": cwd + "/Data/NICER/obs4522010103/ni4522010103.mkf",  # NICER orbital solution

            "lc_path": None,  # RXTE binned data
            "orb_path": None,  # RXTE orbital solution
            "spectrum_path": None,  # RXTE only

            "aster_path": None  # ASTER Labs orbital solution
}
