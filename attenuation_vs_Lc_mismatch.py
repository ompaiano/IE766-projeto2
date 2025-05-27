import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

import coupled

# Waveguide parameters:
n0 = 1.45 # silica refractive index.
n1 = 3.48 # silicon refractive index.
width = 500e-9  # meters.
height = 220e-9 # meters.
# gap = 75e-9
# Optical parameters:
lambida = 1580e-9  # meters.
p = 1
q = 1
# Devive parameters (Max attenuation requirement):
loss_threshold_dB = -0.3
loss_threshold = 10**(loss_threshold_dB/10)
print(loss_threshold)

num_z = int(1e5)
z = np.linspace(0.0, 10e-6, num_z)
num_gaps = 200
gaps = np.linspace(10e-9, 100e-9, num_gaps)
dist_mismatch = []
for gap in gaps:
    Lc, kappa = coupled.solve_coupled_waveguides(n0, n1, width, height, gap, lambida, p, q)
    Pa, Pb = coupled.calc_optical_power_evolution(z, kappa, delta=0)
    dist_mismatch.append(Lc - z[np.argmax(Pb>loss_threshold)])

dist_mismatch = np.array(dist_mismatch)
plt.plot(gaps*1e9, dist_mismatch*1e9, "k-")
plt.xlim(gaps[0]*1e9, gaps[-1]*1e9)
plt.xlabel("gap [nm]")
plt.ylabel("$\\Delta L$ (-0.3 dB) [nm]")
plt.grid(which="both")
plt.grid(which="minor", alpha=.25, color="k", linestyle="--", linewidth=.5)
plt.minorticks_on()
plt.savefig("max_acceptable_dist_mismatch.pdf", bbox_inches="tight")
plt.show()
plt.close()