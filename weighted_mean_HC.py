# This script calculates the "weighted mean" for the HC results table

import numpy as np

def calc_weighted_mean(Delta_e_list, delta_e_list):
    w_unnorm = 1/delta_e_list**2  # unnormalized weights list
    w_norm = w_unnorm/np.sum(w_unnorm)  # normalized weights list
    delta_tot = np.sqrt(np.sum((w_norm*delta_e_list)**2))   # weighted uncertainty
    Delta_tot = np.sum(w_norm*Delta_e_list)  # weighted error
    return Delta_tot, delta_tot


# Values of v0_mkf for the two crossings
v0_v4641 = 7.658770226280014
v0_crab = 7.650021684935115

# Ellipsoid results table
print("Ellipsoid:")
Delta_t, delta_t = calc_weighted_mean(
np.array([-0.12, 0.43, 0.44, 0.16]),
np.array([0.13, 0.18, 0.27, 0.38])
)

print(f"Time: {Delta_t} +/- {delta_t}")

Delta_r = v0_v4641*Delta_t
delta_r = v0_v4641*delta_t

print(f"Position: {Delta_r} +/- {delta_r}")

# Sphere Results table
print("Sphere:")
# Ellipsoid results table
Delta_t, delta_t = calc_weighted_mean(
np.array([1.11, 1.64, 1.66, 1.39]),
np.array([0.13, 0.18, 0.25, 0.39])
)

print(f"Time: {Delta_t} +/- {delta_t}")

Delta_r = v0_v4641*Delta_t
delta_r = v0_v4641*delta_t

print(f"Position: {Delta_r} +/- {delta_r}")
print(Delta_r/Delta_t)

# Results for average density profile
print("Average density (V4641):")

Delta_t, delta_t = calc_weighted_mean(
np.array([0.73, 0.33, 0.33, 0.09]),
np.array([0.17, 0.20, 0.27, 0.39])
)

print(f"Time: {Delta_t} +/- {delta_t}")

Delta_r = v0_v4641*Delta_t
delta_r = v0_v4641*delta_t

print(f"Position: {Delta_r} +/- {delta_r}")

# Results for average density profile (REDO)
print("Average density (Crab):")

Delta_t, delta_t = calc_weighted_mean(
np.array([0.56, 0.29, 0.04, 0.64]),
np.array([0.03, 0.05, 0.10, 0.49])
)

print(f"Time: {Delta_t} +/- {delta_t}")

Delta_r = v0_crab*Delta_t
delta_r = v0_crab*delta_t

print(f"Position: {Delta_r} +/- {delta_r}")

