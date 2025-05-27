import numpy as np

def calc_coupling_length(neff_even, neff_odd, lambida):
    """

    """
    neff_diff = np.abs(neff_even - neff_odd)
    coupling_length = 0.5 * lambida / neff_diff
    return coupling_length

def calc_optical_power_evolution(z, kappa, delta):
    F = calc_max_power_coupling_efficiency(kappa, delta)
    q = np.sqrt( kappa**2 + delta**2 )
    Pa = 1 - F * np.sin(q*z)**2
    Pb = F * np.sin(q*z)**2
    return Pa, Pb