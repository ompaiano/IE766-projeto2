import numpy as np
from numpy import pi
from scipy.constants import lambda2nu, nu2lambda
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

import coupled


coupled.sanity_check_compare_against_okamotos_results()

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

# Sweep parameters:
gaps = np.linspace(20e-9, 400e-9, 50)
Lc = []
for gap in gaps:
    coupling_length = coupled.solve_coupled_waveguides(n0, n1, width, height, gap, lambida, p, q)[0]
    Lc.append(coupling_length)
Lc = np.array(Lc)

plt.plot(gaps*1e9, Lc*1e6)
plt.xlabel("gap [nm]")
plt.ylabel("$L_c$ [$\\mu$m]")
plt.grid(which="both")
plt.grid(which="minor", alpha=.25, color="k", linestyle="--", linewidth=.5)
plt.minorticks_on()
plt.savefig("Lc_vs_gap.pdf", bbox_inches="tight")
plt.show()
plt.close()

# More sweeps, now also on frequency:

# Sweep params:
freq_band = 5e12 # 5 THz.
lambida_i = nu2lambda(lambda2nu(lambida) - freq_band/2)
lambida_f = nu2lambda(lambda2nu(lambida) + freq_band/2)
lambidas = np.linspace(lambida_i, lambida_f, 3, endpoint=True)
gaps = np.linspace(10e-9, 100e-9, 50)
Lcs = [[] for _ in range(len(lambidas))]
linestyles = [":", "-", "--"]
for lc, lambida, lnstyle in zip(Lcs, lambidas[::-1], linestyles):
    for gap in gaps:
        coupling_length = coupled.solve_coupled_waveguides(n0, n1, width, height, gap, lambida, p, q)[0]
        lc.append(coupling_length)
    lc = np.array(lc)
    plt.plot(gaps*1e9, lc*1e6, 
        label=f"$\\lambda_0=$ {lambida*1e6:.3f} $\\mu$m $\\Rightarrow\\nu=$ {lambda2nu(lambida)*1e-12:.2f} THz",
        linestyle=lnstyle)

plt.xlim(gaps[0]*1e9, gaps[-1]*1e9)
plt.ylim(3, 10)
plt.legend()
plt.xlabel("gap [nm]")
plt.ylabel("$L_c$ [$\\mu$m]")
plt.grid(which="both")
plt.grid(which="minor", alpha=.25, color="k", linestyle="--", linewidth=.5)
plt.minorticks_on()
plt.gcf().set_size_inches(6, 5)
plt.savefig("Lc_vs_gap_sweep_freq.pdf", bbox_inches="tight")
plt.show()
plt.close()
    
## assumir um erro de fabricação na diferenca de largura dos guias e ver como isso muda o max power coupling efficiency.