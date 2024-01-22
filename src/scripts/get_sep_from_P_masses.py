#!/usr/local/bin/python

from math import pi

import astropy.units as u
from astropy.constants import G


def get_sep_from_P_masses(M1, M2, P):
    """ M1, M2 in solar masses, P in days """
    # quantities in cgs units
    # convert to cgs
    try:
        M1 = ((M1.value) * u.Msun).to(u.g)
        M2 = ((M2.value) * u.Msun).to(u.g)
        P = ((P.value) * u.day).to(u.s)
    except AttributeError:
        M1 = (M1 * u.Msun).to(u.g)
        M2 = (M2 * u.Msun).to(u.g)
        P = (P * u.day).to(u.s)
    a = (G.to(u.cm ** 3 / (u.g * u.s ** 2)) * (M1 + M2) * P * P / (4 * pi ** 2)) ** (1 / 3.0)
    a = a.to(u.Rsun)
    return a


if __name__ == "__main__":
    M1 = 24
    M2 = 1
    P = 7
    a_HD45485 = get_sep_from_P_masses(M1, M2, P)
    M1 = 15
    M2 = 1.5
    P = 3.5
    a2 = get_sep_from_P_masses(M1, M2, P)
    M1 = 7
    M2 = 1
    P = 5
    a3 = get_sep_from_P_masses(M1, M2, P)
    # print(a_HD45485, a2, a3)
