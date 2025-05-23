import numpy as np
import matplotlib.pyplot as plt

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

z = np.linspace(0.0, 10e-6, 1000)

Lc, kappa = coupled.solve_coupled_waveguides(n0, n1, width, height, gap, lambida, p, q)
Pa, Pb = coupled.calc_optical_power_evolution(z, kappa, delta=0)

plt.plot(z*1e6, Pa)
plt.plot(z*1e6, Pb, "--")
plt.xlabel("z [$\\mu$m]")
plt.ylabel("Normalized power [dimensionless]")
plt.grid(which="both")
plt.grid(which="minor", alpha=.25, color="k", linestyle="--", linewidth=.5)
plt.minorticks_on()
plt.savefig("power_vs_z.pdf", bbox_inches="tight")
plt.show()
plt.close()