# def do_one():
#     altitude_list = np.arange(400, 2100, 100)
#     dt_list = np.zeros_like(altitude_list, float)
#     dr_list = np.zeros_like(dt_list)  # lists containing uncertainties corresponding to altitude_list
#
#     for i, alt in enumerate(altitude_list):
#         sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
#         comp_obj = CurveComparison(sat, hc_type, N)
#         dt_list[i] = comp_obj.dt_e
#         dr_list[i] = comp_obj.dt_e * sat.R_orbit * sat.omega
#
#     # Plot results
#     plt.figure(1)
#     plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dt_list, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.legend()
#
#     plt.figure(2)
#     plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dr_list,
#              label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.show()


# def do_a_bunch(min_alt, max_alt, alt_interval, how_many):
#     altitude_list = np.arange(min_alt, max_alt, alt_interval)
#     dt_list = np.zeros_like(altitude_list, float)
#     dr_list = np.zeros_like(dt_list)
#     for j in range(how_many):
#         #np.random.seed(j)
#         for i, alt in enumerate(altitude_list):
#             sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
#             comp_obj = CurveComparison(sat, hc_type, N)
#             dt_list[i] += comp_obj.dt_e / how_many
#             dr_list[i] += comp_obj.dt_e * sat.R_orbit * sat.omega / how_many
#             print("seed: " + str(j) + " altitude: " + str(alt))
#     plt.figure(1)
#     # plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dt_list, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.legend()
#     plt.savefig("bunches/dt_v_alt.png")
#
#     plt.figure(2)
#     # plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dr_list,
#              label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.savefig("bunches/dr_v_alt.png")
#     plt.show()
#
#     dt = open("bunches/dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     dt.write(str(dt_list))
#     dt.close()
#     dr = open("bunches/dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     dr.write(str(dr_list))
#     dr.close()
#     altlist = open("bunches/alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     altlist.write(str(altitude_list))
#     altlist.close()
#     return 0


# def do_a_bunch_median(min_alt, max_alt, alt_interval, how_many):
#     altitude_list = np.arange(min_alt, max_alt, alt_interval)
#     dt_list = np.zeros((len(altitude_list),how_many), dtype="float")
#     dr_list = np.zeros((len(altitude_list),how_many), dtype="float")
#     for j in range(how_many):
#         for i, alt in enumerate(altitude_list):
#             sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
#             comp_obj = CurveComparison(sat, hc_type, N)
#             dt_list[i][j] = comp_obj.dt_e / how_many
#             dr_list[i][j] = comp_obj.dt_e * sat.R_orbit * sat.omega / how_many
#             print("iteration: " + str(j) + " altitude: " + str(alt))
#
#     dt_median = np.median(dt_list, axis=1)
#     dr_median = np.median(dr_list, axis=1)
#     plt.figure(1)
#     # plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dt_median, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.legend()
#     plt.savefig("bunches/med_dt_v_alt.png")
#
#     plt.figure(2)
#     # plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dr_median, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.savefig("bunches/med_dr_v_alt.png")
#     plt.show()
#
#     dt_data = open("bunches/med_dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     dt_data.write(str(dt_list))
#     dt_data.close()
#     dr_data = open("bunches/med_dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     dr_data.write(str(dr_list))
#     dr_data.close()
#     altlist = open("bunches/med_alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     altlist.write(str(altitude_list))
#     altlist.close()
#     return 0

# def plot_points(min_alt, max_alt, alt_interval, how_many):
#     altitude_list = np.arange(min_alt, max_alt, alt_interval)
#     dt_list = np.zeros((len(altitude_list), how_many), dtype="float")
#     dr_list = np.zeros((len(altitude_list), how_many), dtype="float")
#     for j in range(how_many):
#         for i, alt in enumerate(altitude_list):
#             sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
#             comp_obj = CurveComparison(sat, hc_type, N)
#             dt_list[i][j] = comp_obj.dt_e
#             dr_list[i][j] = comp_obj.dt_e * sat.R_orbit * sat.omega
#             print("iteration: " + str(j) + "/" + str(how_many) + " altitude: " + str(alt) +
#                   ", " + str(
#                 round((j * len(altitude_list) + i + 1) * 100 / (how_many * len(altitude_list)), 2)) + "% complete")
#
#     plt.figure(1)
#     # plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dt_list, '.', label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.legend()
#     plt.savefig("bunches/points_dt_v_alt.png")
#
#     plt.figure(2)
#     # plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
#     plt.plot(altitude_list, dr_list, '.', label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
#     plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
#     plt.xlabel("Orbital altitude (km)")
#     plt.savefig("bunches/points_dr_v_alt.png")
#     plt.show()
#
#     dt_data = open("bunches/points_dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     dt_data.write(str(dt_list))
#     dt_data.close()
#     dr_data = open("bunches/points_dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     dr_data.write(str(dr_list))
#     dr_data.close()
#     altlist = open("bunches/points_alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
#     altlist.write(str(altitude_list))
#     altlist.close()
#     return 0


# plt.plot(altitude_list, dt_66, label=fr"upper 66 data")
# plt.plot(altitude_list, dt_33, label=fr"lower 33 data")
# plt.plot(altitude_list, dt_median, label=fr"median data")
# plt.plot(altitude_list, dt_66fit, label=fr"upper 66 poly fit")
# plt.plot(altitude_list, dt_33fit, label=fr"lower 33 poly fit")
# plt.plot(altitude_list, dt_medianfit, label=fr"median poly fit")
# plt.plot(altitude_list, dt_66sqfit, label=fr"upper 66  sq-1 fit")
# plt.plot(altitude_list, dt_33sqfit, label=fr"lower 33 sq-1 fit")
# plt.plot(altitude_list, dt_mediansqfit, label=fr"median sq-1 fit")
# plt.plot(altitude_list, dt_66_invlog_fit, label=fr"66 invlog fit")
# plt.plot(altitude_list, dt_33_invlog_fit, label=fr"33 invlog fit")
# error = np.ones(len(altitude_list))*(dt_66_invlog_fit-dt_33_invlog_fit) / 2
# plt.errorbar(altitude_list, dt_median_invlog_fit, yerr=error, fmt="+")
# plt.plot(altitude_list, dt_median_exponential_fit, label=fr"median exponential fit")


# plt.plot(altitude_list, dr_66, label=fr"upper 66 data")
# plt.plot(altitude_list, dr_33, label=fr"lower 33 data")
# plt.plot(altitude_list, dr_66fit, label=fr"upper 66 poly fit")
# plt.plot(altitude_list, dr_33fit, label=fr"lower 33 poly fit")
# plt.plot(altitude_list, dr_medianfit, label=fr"median poly fit")
# plt.plot(altitude_list, dr_66sqfit, label=fr"upper 66  sq-1 fit")
# plt.plot(altitude_list, dr_33sqfit, label=fr"lower 33 sq-1 fit")
# plt.plot(altitude_list, dr_mediansqfit, label=fr"median sq-1 fit")
# plt.plot(altitude_list, dr_median, label=fr"median data")
# plt.plot(altitude_list, dr_66_invlog_fit, label=fr"std. dev", color="red")
# plt.plot(altitude_list, dr_33_invlog_fit, color="red")
# dt_median_exponential_fit = plot_exponential_fit(altitude_list, dt_median, "true")