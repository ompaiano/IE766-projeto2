import numpy as np
from numpy import pi
from scipy.constants import lambda2nu, nu2lambda
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

import coupled
import supermodes as sm


# Plot Parameters:
colors = ["g", "b", "k", "c", "m"]


# Coupled mode theory stuff:
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
plt.savefig("Lc_vs_gap__cpt_vs_sm.pdf", bbox_inches="tight")
plt.show()
plt.close()

# More sweeps, now also on frequency:

# Sweep params:
freq_band = 10e12 # 5 THz.
lambida_i = nu2lambda(lambda2nu(lambida) - 2*freq_band/2)
lambida_f = nu2lambda(lambda2nu(lambida) + 2*freq_band/2)
num_freqs = 5
lambidas = np.linspace(lambida_i, lambida_f, num_freqs, endpoint=True)
lambidas = 1e-9*np.array([1520, 1540, 1560, 1580, 1600])
gaps = np.linspace(5e-9, 105e-9, 50)
Lcs = [[] for _ in range(len(lambidas))]
linestyles = [":", (0, (3, 5, 1, 5, 1, 5)), "-", "-.", "--"] # for more personalized linestyles, check: https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html

for lc, lambida, c, lnstyle in zip(Lcs, lambidas[::-1], colors, linestyles):
    for gap in gaps:
        coupling_length = coupled.solve_coupled_waveguides(n0, n1, width, height, gap, lambida, p, q)[0]
        lc.append(coupling_length)
    lc = np.array(lc)
    plt.plot(gaps*1e9, lc*1e6, 
        label=f"CMT, $\\lambda_0=$ {lambida*1e6:.3f} $\\mu$m", #  $\\Rightarrow\\nu=$ {lambda2nu(lambida)*1e-12:.2f} THz",
        linestyle=lnstyle, c=c)

# Now we do supermodes stuff:
markers = ["o", "x", "*", "s", "P"]
wls = [1520, 1540, 1560, 1580, 1600]
sm_data = sm.get_data()

for wl, c, mk in zip(wls[::-1], colors, markers):
    aux_data = sm_data[ sm_data.wavelength == wl ]
    print(type(aux_data.wavelength))
    print(aux_data.neff_even.values)
    Lc = sm.calc_coupling_length(
        aux_data.neff_even.values.real, 
        aux_data.neff_odd.values.real, 
        aux_data.wavelength * 1e-9)
    plt.plot(aux_data.gap * 1e9, Lc*1e6, c+mk,
        label=f"SM, $\\lambda_0=$ {wl/1e3:.3f} $\\mu$m")


plt.xlim(gaps[0]*1e9, .98*gaps[-1]*1e9)
# plt.ylim(3, 10)
plt.legend(fontsize=11) # loc="lower right")
plt.xlabel("gap [nm]")
plt.ylabel("$L_c$ [$\\mu$m]")
plt.grid(which="both")
plt.grid(which="minor", alpha=.25, color="k", linestyle="--", linewidth=.5)
plt.minorticks_on()
plt.gcf().set_size_inches(6, 6)
plt.savefig("Lc_vs_gap__sweep_freq__cmt_vs_sm.pdf", bbox_inches="tight")
plt.show()
plt.close()
    
## assumir um erro de fabricação na diferenca de largura dos guias e ver como isso muda o max power coupling efficiency.