import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

import coupled
import supermodes


# Waveguide parameters:
n0 = 1.45 # silica refractive index.
n1 = 3.48 # silicon refractive index.
width = 500e-9  # meters.
height = 220e-9 # meters.
gap = 75e-9

# Optical parameters:
lambida = 1560e-9  # meters.
p = 1
q = 1

num_z = int(1e3)
z = np.linspace(0.0, 12e-6, num_z)

Lc, kappa = coupled.solve_coupled_waveguides(n0, n1, width, height, gap, lambida, p, q)
Pa, Pb = coupled.calc_optical_power_evolution(z, kappa, delta=0)

plt.plot(z*1e6, Pa, label="$P_a(z)$")
plt.plot(z*1e6, Pb, "--", label="$P_b(z)$")
plt.xlim(z[0]*1e6, z[-1]*1e6)
plt.legend()
plt.xlabel("z [$\\mu$m]")
plt.ylabel("Normalized power")
plt.grid(which="both")
plt.grid(which="minor", alpha=.25, color="k", linestyle="--", linewidth=.5)
plt.minorticks_on()
plt.savefig("power_vs_z.pdf", bbox_inches="tight")
plt.show()
plt.close()