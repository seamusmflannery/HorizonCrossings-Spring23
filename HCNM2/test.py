# This script shows an example of how to manually change the working directory if needed

from pathlib import Path
import os

os.chdir(str(Path(__file__).parents[1]))
from ObservationDictionaries.v4641NICER import v4641NICER
from Modules.OrbitModel import OrbitModel

r_array, t_array = OrbitModel.define_orbit_model(v4641NICER, "mkf", time_step=0.01)

print(r_array)