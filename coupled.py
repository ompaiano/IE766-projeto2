"""Module to solve rectangular waveguides coupling by means of coupled mode 
theory.

See Okamoto, chapter 4.

"""

import numpy as np
from numpy import pi
from scipy import optimize
import matplotlib.pyplot as plt


def kx_system(ks, p, a, n0, n1, k):
    eq1 = (p-1)*pi/2.0 - ks[0]*a + np.atan( (n1**2 *ks[1])/(n0**2 *ks[0]) )
    eq2 = k**2  * (n1**2 - n0**2) - ks[0]**2 - ks[1]**2
    return [eq1, eq2] 

def ky_system(ks, p, a, n0, n1, k):
    eq1 = (p-1)*pi/2.0 - ks[0]*a + np.atan(ks[1]/ks[0])
    eq2 = k**2  * (n1**2 - n0**2) - ks[0]**2 - ks[1]**2
    return [eq1, eq2] 

def calc_kappa(k, kx, gx, a, n0, n1, gap):
    """Calculate kappa (codirectional mode coupling in coupled mode theory) 
    for the particular case of eq. 4.91 of Okamoto's 'Fundamentals of Optical 
    Waveguides'.

    This method is approximate and applies for waveguides and mode satisfying 
    the following assumptions: 
        1. $E^x_{11}$ mode.
        2. rectangular waveguide.

    Parameters
    ----------
    k: float
        Wavenumber (in vacuum).
    kx: float
        Componente of the effective wave-vector in the x direction.
    gx: float
        Exponential decay wavenumber (check the actual name for this parameter).
        for the x direction.
    a: float
        Half of the waveguide width in the x direction.
    n0: float
        Cladding refractive index.
    n1: float
        Core refractive index.
    gap: float
        Gap between the two waveguides.

    Returns
    -------
    float
    
    """
    delta = 0.5 * (n1**2 - n0**2) / n1**2 
    normalized_frequency = k * n1 * a * np.sqrt(2*delta)
    kappa = np.sqrt(2*delta) * a * (a* kx * gx)**2 * np.exp(-gx*gap)
    kappa /= (1 + gx*a) * normalized_frequency**3
    return kappa, delta, normalized_frequency

def calc_coupling_length(kappa, delta):
    """Calculate the length for maximum power transfer bewteen waveguides.

    """
    Lc = 0.5*pi / np.sqrt(kappa**2 + delta**2)
    return Lc

def solve_coupled_waveguides(
    n0, n1, width, height, gap, lambida, p, q, verbose=True):
    """Solve co-directionally coupled waveguides through coupled waveguide 
    theory.

    TO DO: RETURN ALL PARAMETERS!!

    """
    k = 2*pi / lambida # wavenumber (for vacuum).
    n_avg = np.mean((n0, n1))
    ksx0 = (n_avg*k, k)
    ksy0 = (n_avg*k, k)

    solx = optimize.root(kx_system, ksx0, args=(p, width/2, n0, n1, k))
    soly = optimize.root(ky_system, ksy0, args=(q, height/2, n0, n1, k))
    kx, gamma_x = solx.x
    ky, gamma_y = soly.x
    beta = np.sqrt( (k*n1)**2 - kx**2 - ky**2 )
    kappa, deltao, norm_freq = calc_kappa(k, kx, gamma_x, width/2, n0, n1, gap)
    coupling_length = calc_coupling_length(kappa, delta=0)

    if verbose:
        # PRINTAR OS PARAMETROS TB!!
        print("Status da solução na direção x: ", solx.success, "\n\tMessage: ", solx.message)
        print("Status da solução na direção y: ", soly.success, "\n\tMessage: ", soly.message)
        print("kx/k, gamma_x/k: ", kx/k, gamma_x/k)
        print("ky/k, gamma_y/k: ", ky/k, gamma_y/k)
        print("kx*a, gamma_x*a: ", kx*width/2, gamma_x*width/2)
        print("effective index (beta/k): ", beta/k)
        print("kappa, deltao, norm_freq: ", kappa, deltao, norm_freq)
        print("Coupling length [micrometers]: ", coupling_length*1e6)

    return coupling_length, kappa


def sanity_check_compare_against_okamotos_results():
    print("Okamoto parameters: ")
    # Okamoto parameters (sanity check, it works!):
    n0 = 1.45
    n1 = 1.4544
    width = 8e-6
    height = width
    D = 3*width/2
    gap = D - width
    lambida = 1550e-9  # meters.
    p = 1
    q = 1
    solve_coupled_waveguides(n0, n1, width, height, gap, lambida, p, q, verbose=True)

def calc_max_power_coupling_efficiency(kappa, delta):
    """Calculate the maximum power coupling efficiency.
    """
    F = 1 / (1 + (delta/kappa)**2)
    return F

def calc_optical_power_evolution(z, kappa, delta):
    F = calc_max_power_coupling_efficiency(kappa, delta)
    q = np.sqrt( kappa**2 + delta**2 )
    Pa = 1 - F * np.sin(q*z)**2
    Pb = F * np.sin(q*z)**2
    return Pa, Pb

if __name__ == "__main__":
    # Waveguide parameters:
    n0 = 1.45 # silica refractive index.
    n1 = 3.48 # silicon refractive index.
    width = 500e-9  # meters.
    height = 220e-9 # meters.
    gap = 75e-9

    # Optical parameters:
    lambida = 1560e-9  # meters.
    k = 2*pi / lambida # wavenumber (for vacuum).
    p = 1
    q = 1

    a_avg = np.mean((width, height)) / 2
    n_avg = np.mean((n0, n1))
    ksx0 = (n_avg*k, 1/a_avg)
    ksy0 = (n_avg*k, 1/a_avg)

    ksx0 = (n_avg*k, k)
    ksy0 = (n_avg*k, k)

    solx = optimize.root(kx_system, ksx0, args=(p, width/2, n0, n1, k))
    soly = optimize.root(ky_system, ksy0, args=(q, height/2, n0, n1, k))
    kx, gamma_x = solx.x
    ky, gamma_y = soly.x
    print("Status da solução na direção x: ", solx.success, "\n\tMessage: ", solx.message)
    print("Status da solução na direção y: ", soly.success, "\n\tMessage: ", soly.message)

    # Print the solution
    print("kx/k, gamma_x/k: ", kx/k, gamma_x/k)
    print("ky/k, gamma_y/k: ", ky/k, gamma_y/k)
    print(kx*width/2, gamma_x*width/2)

    beta = np.sqrt( (k*n1)**2 - kx**2 - ky**2 )
    print("effective index (beta/k): ", beta/k)

    # Coupled mode theory stuff
    kappa, deltao, norm_freq = calc_kappa(k, kx, gamma_x, width/2, n0, n1, gap)
    print("kappa, deltao, norm_freq: ", kappa, deltao, norm_freq)

    coupling_length = calc_coupling_length(kappa, delta=0)
    print("Coupling length [micrometers]: ", coupling_length*1e6)