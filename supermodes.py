import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


def get_data(lambida=None, gap=None):
    """Get supermodes data.

    Parameters
    ----------
    lambida: int
        Wavelength in nanometers.
        Defaults to None, which results in retrieving information for all 
        wavelengths. 
    gap: int
        Gap between waveguides in nanometers.
        Defaults to None, which results in retrieving information for all gap 
        sizes. 

    Returns
    -------
    pandas.dataframe
        The keys for the dataframe are:
            - "wavelength";
            - "gap";
            - "neff_even";
            - "neff_odd".
    """
    data = _read_supermodes_files(path=None)
    if lambida is not None:
        data = data[ data["lambda"] == lambida ]
    if gap is not None:
        data = data[ data["gap"] == gap ]
    return data

def _read_supermodes_files(path=None):
    """

    Parameters
    ----------
    path: str
        Path to root of the supermodes data.

    """
    if path is None:
        path = "./lumerical/supermodes/gap-sweep-and-frequency/"
    wavelength = []
    gap = []
    neff_even = []
    neff_odd = []
    wavelength_paths = [f.path for f in os.scandir(path) if f.is_dir()]
    for wl_path in wavelength_paths:
        g = np.loadtxt(wl_path + "/gaps.txt")
        ne = np.loadtxt(wl_path + "/neffs-TE0_even.txt", dtype=str)
        no = np.loadtxt(wl_path + "/neffs-TE0_odd.txt", dtype=str)
        wl = [ int(wl_path.split("/")[-1]) for _ in range(len(g)) ]
        gap.extend(g)
        neff_even.extend(ne)
        neff_odd.extend(no)
        wavelength.extend(wl)
    neff_even = [np.complex128(n.replace("i", "j")) for n in neff_even]
    neff_odd  = [np.complex128(n.replace("i", "j")) for n in neff_odd]
    data = {
        "wavelength": wavelength,
        "gap": gap,
        "neff_even": neff_even,
        "neff_odd": neff_odd
    }
    return pd.DataFrame.from_dict(data)


if __name__ == "__main__":

    data = get_data()

    data_1560 = data[ data.wavelength == 1560 ]

    print(type(data_1560.wavelength))

    print(data_1560.neff_even.values)

    Lc = calc_coupling_length(
        data_1560.neff_even.values.real, 
        data_1560.neff_odd.values.real, 
        data_1560.wavelength * 1e-9)

    print(Lc)

    plt.plot(data_1560.gap * 1e9, Lc*1e6, "ko")
    plt.xlabel("gap [nm]")
    plt.ylabel("$L_c$ [$\\mu$m]")
    plt.grid(which="both")
    plt.grid(which="minor", alpha=.25, color="k", linestyle="--", linewidth=.5)
    plt.minorticks_on()
    plt.gcf().set_size_inches(6, 5)
    plt.show()

