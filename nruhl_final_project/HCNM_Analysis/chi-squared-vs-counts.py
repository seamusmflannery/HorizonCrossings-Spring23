# This script plots chisq+1 stdErr vs unnatenuated count rate for the results in Breck paper table

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def sqrt_inverse(x, a):
    return a/np.sqrt(x)

counts = [19, 40, 89, 244, 275, 584, 1406, 5378]
stdErr = [0.38, 0.27, 0.20, 0.14, 0.09, 0.07, 0.05, 0.03]  # delta t_{hc}

# seperate lists so we can plot with different colors
counts_v4 = [19, 40, 89, 244]
stdErr_v4 = [0.38, 0.27, 0.20, 0.14]
counts_crab = [275, 584, 1406, 5378]
stdErr_crab = [0.09, 0.07, 0.05, 0.03]

xdata = np.linspace(counts[0],counts[-1], 300)
b, b_err = curve_fit(sqrt_inverse,counts, stdErr)

plt.plot(xdata, sqrt_inverse(xdata, b), label=r"$b/ \sqrt{N_0}$ trendline, b= %.2f" % b)
# plt.plot(counts, stdErr, ".", label=r"$\delta t_{hc}$ from Crab and V4641 Sgr horizon crossings")
plt.plot(counts_v4, stdErr_v4, ".", label=r"$\delta t_{hc}$ from V4641 Sgr horizon crossing")
plt.plot(counts_crab, stdErr_crab, ".", label=r"$\delta t_{hc}$ from Crab Nebula horizon crossing")
# plt.title(r"Standard Error vs Unnatenuated Counts")
plt.ylabel(r"$\chi^2+1$ standard error, $\delta t_{hc}$ (sec)")
plt.xlabel(r"Number of unnatenuated counts, $N_0$")
plt.legend()

plt.show()
