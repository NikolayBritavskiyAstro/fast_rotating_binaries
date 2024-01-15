import numpy as np
from astropy import units as u
from get_sep_from_P_masses import get_sep_from_P_masses
from math import pi
import paths
from Eggleton83 import Roche_lobes, Roche_ratios
from classify_trip import mlp_classifier

def get_Rzams_from_mass_radius_rel(M):
    """Rough mass-radius R = Rsun (M/Msun)^exp relation for main sequence massive stars
    Note that this assumes polytropic homogeneous stars.
    Typically if M <= 1.4: exp = 0.9 else: exp = 0.6

    Parameters:
    ----------
    M: `float` or np.array, mass of the star in Msun units
    exp: `float` exponent, default to approximate value for massive stars

    Returns:
    -------
    radius: `float` or np.array in Rsun units
    """
    try:
        mval = M.value
    except AttributeError:
        mval = M
    try:
        # mval is a scalar
        if mval<= 1.4:
            exp = 0.9
        else:
            exp = 0.6
    except ValueError:
        # mval is an array
        i_low_mass = mval <= 1.4
        exp = np.asarray([0.6]*len(mval), dtype=float)
        exp[i_low_mass] = 0.9
    return (mval ** exp) * u.Rsun

def get_ZAMS_periastron_RLOF_min_a(q, mtot, e=0):
    """Provide the minimum orbital separation so that neither stars
    fill their Roche lobe at ZAMS. If a non-zero eccentricity value is
    provided uses the Roche radii at periastron (i.e. uses a(1-e)
    instead of the separation a in the Eggleton 1983 formula)

    Parameters:
    -----------
    q: `float` or `np.array` mass ratio
    mtot: `float` or `np.array` total mass of the binary
    e: `float` or `np.array` eccentricity of the binary

    Returns:
    --------
    min_a: `float` or `np.array` minimum orbital separation such as a> RL1+RL2 at periastron, in solar radii
    """
    # solve for masses
    m1zams = mtot/(1+q)
    m2zams = m1zams*q
    # assign units
    try:
        m1zams = m1zams.value*u.Msun
        m2zams = m2zams.value*u.Msun
    except AttributeError:
        m1zams = m1zams*u.Msun
        m2zams = m2zams*u.Msun
    # get ZAMS radii from analytic mass-radius relation
    rzams1 = get_Rzams_from_mass_radius_rel(m1zams)
    rzams2 = get_Rzams_from_mass_radius_rel(m2zams)
    # get fi = R_RLi/(a_i(1-e_in)) or in other words: R_RLi = fi*a*(1-e)
    # TheRoche geometry only applies to circular binaries but we consider
    # here the separation at periastron to compensate for this
    # (approximate the ellipse with the circle of periastron radius)
    f1, f2 = Roche_ratios(m1zams.value, m2zams.value)
    # condition we want:
    # periastron passage does not cause RLOF
    # a_in such that rzams1 < f1*a_in*(1-e_in) and rzams2 < f2*a(1-e_in)
    min_a_in = np.minimum(rzams1/(f1*(1.-e)), rzams2/(f2*(1.-e)))
    return min_a_in # same units as radii, should be Rsun

def binary(M1, M2, P, e=0, name=None):
    """ Assign units and return values"""
    try:
        binary_system =  {"M1": M1.value*u.Msun,
                          "M2": M2.value*u.Msun,
                          "P" : P.value*u.day,
                          "a" : get_sep_from_P_masses(M1.value*u.Msun, M2.value*u.Msun, P.value*u.day).to(u.Rsun),
                          "e" : e,
                          "name": name}
    except AttributeError:
        binary_system =  {"M1": M1*u.Msun,
                          "M2": M2*u.Msun,
                          "P" : P*u.day,
                          "a" : get_sep_from_P_masses(M1*u.Msun, M2*u.Msun, P*u.day).to(u.Rsun),
                          "e" : e,
                          "name": name}
    return binary_system


def draw_triples(observed_binary, f_dm=0.0, N_triples=10):
    """ produce a triple system that can be accepted as initial conditions
    Parameters:
    ----------
    observed_binary: `dict`, dictionary containing the properties of the binary
                             observed we want to produce with a triple merger
    f_dm: `float`, fraction of mass lost at the merger moment, defaults to zero

    """
    m_post_merger = observed_binary["M1"].value # post-merger mass  == present day primary mass (neglect wind)
    m1zams_plus_m2zams =  m_post_merger/(1-f_dm) # again neglect pre-merger wind
    q_out = [observed_binary["M2"].value/m1zams_plus_m2zams]*np.ones(N_triples)
    max_a_out = observed_binary["a"].value
    # draw inner binary eccentricity
    e_in = np.random.random(N_triples)
    # draw inner binary mass ratio
    q_in = np.random.random(N_triples)
    min_a_in = get_ZAMS_periastron_RLOF_min_a(q_in, m1zams_plus_m2zams, e_in).value
    # if min_a_in > max_a_out:
    # Count all these systems as unstable
    # because they can't even form. They will be assigned a P_unstable=-5
    i_too_large_min_a_in = min_a_in >= max_a_out
    min_a_in[i_too_large_min_a_in] = max_a_out
    a_in = np.random.uniform(min_a_in, max_a_out, N_triples)
    a_out = np.random.uniform(a_in, max_a_out, N_triples)
    inclination = np.random.uniform(0,pi, N_triples)
    e_out = observed_binary["e"]*np.ones(N_triples)
    triple_list = np.asarray((q_in, q_out, a_in/a_out, e_in, e_out, inclination)).T
    return triple_list, i_too_large_min_a_in


def ML_stability_test(triple, ML_model="./mlp_model_trip_ghost.pkl"):
    """ Uses Pavan Vynatheya et al. 2023 "ghost orbits" classifier (latest version
    has of Sep 9th, 2023) to assign a probability of being dynamically
    unstable."""
    # print(triple)
    q_in = triple[0]
    q_out = triple[1]
    a_ratio = triple[2]
    e_in = triple[3]
    e_out = triple[4]
    inclination = triple[5]
    stability = mlp_classifier(ML_model, q_in, q_out, a_ratio, e_in, e_out, inclination)
    # print(stability)
    return stability

def wrapper(system):
    print(system["name"])
    file_name = str(paths.static)+"/"+system["name"]+"_"+"triples"+".txt"
    with open(file_name, "w") as output_file:
        output_file.writelines("q_in\t\t\tq_out\t\t\ta_ratio\t\t\te_in\t\t\te_out\t\t\tinclination\t\t\tf_dm\t\t\tP_unstable"+"\n")
        for f_dm in [0, 0.1]:
            triple_list, i_too_large_min_a_in = draw_triples(system, f_dm=f_dm, N_triples=N_samples)
            # Parallel(delayed(wrapper)(triple) for triple in triple_list)
            for i, triple in enumerate(triple_list):
                if i_too_large_min_a_in[i]:
                    stable = -5
                    # print(triple, i_too_large_min_a_in[i], "True")
                else:
                    # min(a_in) was < max_a_out
                    stable = ML_stability_test(triple)
                    # print(triple, stable)
                output_file.writelines("\t\t\t".join(str(x) for x in triple)+"\t\t\t"+str(f_dm)+"\t\t\t"+str(float(stable))+"\n")
            print("done with f_dm =", f_dm)
        print("------------------")
        print("done with", system["name"], "output at", file_name)
        print("------------------")


if __name__ == "__main__":
    """
    This script draws `N_samples` triples randomly that hypothetically may
    be progenitors of the observed SB1 binaries with a fast rotator (but
    the sampling caps a_in to max_a_out=a_observed). It uses Pavan
    Vynatheya et al. 2022, 2023 "ghost orbits" classifier (latest version
    has of Sep 9th, 2023) to assign a probability of being dynamically
    unstable.
    """
    N_samples = int(1e6) # number of triples to draw
    HD46485 = binary(24., 1., 7.0, 0.033, "HD46485")
    HD191595 = binary(15., 1.2, 3.6, 0.0, "HD191595")
    HD25631 = binary(7., 1, 5.2, 0.0, "HD25631")
    systems = [HD46485, HD191595, HD25631]
    Parallel(n_jobs=3)(delayed(wrapper)(system) for system in systems)
