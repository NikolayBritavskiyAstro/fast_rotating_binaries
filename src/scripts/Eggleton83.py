#!/usr/bin/python

"""
This uses the approximation from Eggleton 1983 to the Roche size
https://ui.adsabs.harvard.edu/abs/1983ApJ...268..368E/abstract
"""
import numpy as np


def func(q):
    """fitting formula by Eggleton 1983"""
    return 0.49 * q ** (2.0 / 3) / (0.6 * q ** (2.0 / 3) + np.log(1 + q ** (2.0 / 3)))


def Roche_ratios(m1, m2):
    q_donor = m1 / m2
    rl_donor = func(q_donor)
    q_accretor = m2 / m1
    rl_accretor = func(q_accretor)
    return rl_donor, rl_accretor


def Roche_lobes(m1, m2, a, ecc=0):
    """
    Parameters:
    ----------
    m1: `float` stellar mass arbitrary units
    m2: `float` stellar mass same unit as above
    a:  `float` separation arbitrary units
    e: ` float` eccentricity, if provided returns rl at periastron
    Returns:
    -------
    rl1: `float` Roche lobe size of the star of mass m1
    rl2: `float` Roche lobe size of the star of mass m2
    the value is in the same units as the input separation
    """
    q_donor = m1 / m2
    rl_donor = func(q_donor) * a * (1 - ecc)
    q_accretor = m2 / m1
    rl_accretor = func(q_accretor) * a * (1 - ecc)
    return rl_donor, rl_accretor


if __name__ == "__main__":
    M1 = float(input("M1 in solar masses? "))
    M2 = float(input("M2 in solar masses? "))
    a = float(input("M2 in solar radii? "))
    rl1, rl2 = Roche_lobes(M1, M2, a)
    print("Roche lobe star 1:", rl1, "Rsun")
    print("Roche lobe star 2:", rl2, "Rsun")
