import numpy as np
from math import pi
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from astropy.constants import G
import paths
from ML_stability import binary, get_ZAMS_periastron_RLOF_min_a, get_Rzams_from_mass_radius_rel
import matplotlib as mpl



def plot_min_a_in_vs_q_in_e_in(M_post_merger=24.0 *u.Msun,
                               m3=1.0*u.Msun, f_dm = 0, N=11, ax=None,
                               a_observed=None, cmap=plt.cm.viridis,
                               norm="linear"):
    """
    N.B.: widest separation is for q_in=1 e_in~0.99 to avoid RLOF at periastron
    """
    # correct for mass lost at merger
    mtot = M_post_merger/(1-f_dm)

    q_ins = np.linspace(0.1, 1.0, N)
    e_ins = np.linspace(0, 0.9999, N)
    min_a_in = np.zeros((N,N))
    for q_in in q_ins:
        i = np.argmin(np.absolute(q_ins-q_in))
        for e_in in e_ins:
            j = np.argmin(np.absolute(e_ins-e_in))
            min_a_in[j, i] = get_ZAMS_periastron_RLOF_min_a(q_in, mtot, e=e_in).value
    # reassign units
    min_a_in *= u.Rsun
    if not ax:
        fig = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(100, 100)
        ax = fig.add_subplot(gs[:, :90])
        cax = fig.add_subplot(gs[:, :95])
    # make plot
    p = ax.pcolormesh(q_ins, e_ins, min_a_in.value, norm=norm, cmap=cmap, shading='auto')
    ax.set_xlim(0.1,1)
    ax.set_ylim(0,1)
    if a_observed:
        max_Rzams_sum = 0
        for q_in in q_ins:
            m1zams = mtot/(1+q_in)
            m2zams = m1zams*q_in
            rzams1 = get_Rzams_from_mass_radius_rel(m1zams.value*u.Msun)
            rzams2 = get_Rzams_from_mass_radius_rel(m2zams.value*u.Msun)
            if rzams2+rzams1 >= max_Rzams_sum: max_Rzams_sum = rzams2+rzams1
        contour = ax.contour(q_ins, e_ins, min_a_in.value,
                             levels=[(max_Rzams_sum).value,
                                     a_observed.value], colors=["#808080", "r"])
        fmt = {}
        strings = [r"$contact$",  r"$a_\mathrm{obs}$"]
        for l, s in zip(contour.levels, strings):
            fmt[l] = s
        ax.clabel(contour, inline=1, fontsize=20, fmt=fmt, inline_spacing=1, manual=[(0.3, 0.2), (0.75,0.4)])
        contourf = ax.contourf(q_ins, e_ins, min_a_in.value,
                             levels=[-1, (max_Rzams_sum).value,], colors=["k"])
    return p



def plot_min_a_in(observed_binary, fig_name=None, title=None):
    fig = plt.figure(figsize=(10,5))
    gs = gridspec.GridSpec(100, 110)
    ax1 = fig.add_subplot(gs[:, :40])
    ax2 = fig.add_subplot(gs[:, 60:100])
    # ax3 = fig.add_subplot(gs[:, 100:140])
    cax = fig.add_subplot(gs[:, 105:]) #145:])
    axes = [ax1,ax2] # ,ax3]
    fracs = [0, 0.1] # , 0.5]
    norm = mpl.colors.Normalize(5, 50)
    for f_dm in fracs: # min a_out
        ax = axes[fracs.index(f_dm)]
        # ax.set_title(r"$f_\Delta="+f"{f_dm}$", size=30)
        ax.text(0.95,0.05, r"$f_\Delta="+f"{f_dm}$", va="bottom", ha="right", fontsize=20, transform=ax.transAxes, bbox=dict(facecolor='w', alpha=0.5, edgecolor='black', boxstyle='round,pad=0.1'))
        p = plot_min_a_in_vs_q_in_e_in(M_post_merger =  observed_binary["M1"],
                                       m3 = observed_binary["M2"],
                                      f_dm=f_dm, ax=ax, N=50,
                                      a_observed = observed_binary["a"],
                                      cmap = plt.cm.viridis,
                                      norm = norm)
        ax.set_ylabel(r"$e_\mathrm{in}$",fontsize=25)
        ax.set_xlabel(r"$q_\mathrm{in} = m_2/m_1$",fontsize=25)
    fig.colorbar(p, cax=cax, extend="both")
    cax.set_ylabel(r"$\min(a_{in}) \ [R_\odot]$",fontsize=25)
    cax.axhline(observed_binary["a"].value, 0,1, color="r")
    if title: cax.set_title(title, size=30)
    if fig_name:
        plt.savefig(paths.figures/fig_name, dpi=300)
        plt.close()
    else:
        plt.show()


def plot_min_a_out(observed_binary, fig_name=None, title=None):
    fig = plt.figure(figsize=(20,5))
    gs = gridspec.GridSpec(100, 100)
    ax1 = fig.add_subplot(gs[:, :40])
    ax2 = fig.add_subplot(gs[:, 50:90])
    # ax3 = fig.add_subplot(gs[:, 100:140])
    cax = fig.add_subplot(gs[:, 95:]) # 145:])
    axes = [ax1,ax2] # ,ax3]
    fracs = [0, 0.1]
    norm = mpl.colors.Normalize(5, 50)
    for f_dm in fracs: # min a_out
        ax = axes[fracs.index(f_dm)]
        ax.set_title(r"$f_\Delta="+f"{f_dm}$", size=30)
        p = plot_min_a_out_vs_q_in_e_in(M_post_merger =  observed_binary["M1"],
                                        m3 = observed_binary["M2"],
                                        f_dm=f_dm, ax=ax, N=50,
                                        a_observed = observed_binary["a"],
                                        cmap = plt.cm.viridis,
                                        norm = norm)
        # ax.set_title(r"$f_\Delta=$"+f"{f_dm}", size=30)
    ax1.set_ylabel(r"$e_\mathrm{in}$",fontsize=25)
    ax2.set_xlabel(r"$q_\mathrm{in}$",fontsize=25)
    fig.colorbar(p, cax=cax, extend="both")
    cax.set_ylabel(r"$\min(a_{out}) \ [R_\odot]$",fontsize=25)
    cax.axhline(observed_binary["a"].value, 0,1, color="r")
    if title: cax.set_title(title, size=30)
    if fig_name:
        plt.savefig(paths.figures/fig_name, dpi=300, layout="tight")
        plt.close()
    else:
        plt.show()




# def get_lower_limit_a_out(e_out, rel_incl, q_in, e_in, binary, f_dm=0):
#     """ get the lower limit on the initial ZAMS a_out so that:

#     - no ZAMS RLOF at periastron of outer binary
#     - a_out >= C (RL1+LR2) with C Maardling & Aarseth 1999 criterion for stability

#     """
#     # calculate individual masses
#     m_post_merger = binary["M1"].value
#     m1zams_plus_m2zams =  m_post_merger/(1-f_dm)
#     # define outer binary mass ratio
#     q_out = binary["M2"].value/m1zams_plus_m2zams
#     # get critical separation for non-chaos from Maardling & Aarseth 1999
#     a_crit = critical_a_out_div_a_in_ratio(e_out, q_out, rel_incl)
#     # get minimum separation to avoid inner binary RLOF
#     a_in_min = get_ZAMS_periastron_RLOF_min_a(q_in, m1zams_plus_m2zams, e = e_in)
#     # get minimum outer separation to avoid outer binary RLOF
#     a_out_min_chaos = a_crit*a_in_min
#     # check for outer binary RLOF at periastron
#     a_out_min_outer_RLOF = get_ZAMS_periastron_RLOF_min_a(q_out, m1zams_plus_m2zams+binary["M2"].value, e=e_out)
#     a_out_min = np.minimum(a_out_min_outer_RLOF, a_out_min_chaos)
#     return a_out_min



def get_max_vesc(q, mtot):
    G_cgs = G.to(u.cm**3/(u.g*u.s**2))
    min_a_in = get_ZAMS_periastron_RLOF_min_a(q, mtot)
    v_esc = np.sqrt(2*G_cgs*(mtot*u.Msun).to(u.g)/min_a_in.to(u.cm))
    return v_esc


def plot_v_esc_in_vorb3_vs_q_in(mtot=24.0, m3=1, rel_incl=0, e_out=0,
                                lw=3, color="C{}".format(0), ax=None):
    if not ax:
        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(100, 100)
        ax = fig.add_subplot(gs[:, :])
    G_cgs = G.to(u.cm**3/(u.g*u.s**2))
    q_in = np.linspace(0.1, 1, 100)
    # ax.set_title(f"{mtot}"+r"$\,M_\odot$,"+r"$m_3=$"+f"{m3}"+r"$a_\mathrm{out} = a_\mathrm{crit}(R_\mathrm{1,zams}+R_\mathrm{2,zams})$", size=30)
    ax.plot(q_in, (get_max_vesc(q_in, mtot=mtot).to(u.km/u.s)).value,
            color=color, lw=lw)
    # for legend
    # ax.plot(np.nan, np.nan, ls='--', c='k', label="$v_\mathrm{orb,3}$")
    f_dms = [0.01, 0.1, 0.5]
    line_styles = [':', '--', '-.']
    for f_dm in f_dms:
        # use a_crit for Mardling & Aarset to pick a_out_min
        min_a_in = get_ZAMS_periastron_RLOF_min_a(q_in, mtot)
        critical_ratio = critical_a_out_div_a_in_ratio(e_out, m3/(mtot/(1-f_dm)), rel_incl)
        a_out = critical_ratio*min_a_in.to(u.cm)
        v_orb_3 = np.sqrt(G_cgs*((mtot/(1-f_dm)+m3)*u.Msun).to(u.g)/a_out)
        p = ax.plot(q_in, (v_orb_3.to(u.km/u.s)).value,
                ls=line_styles[f_dms.index(f_dm)],
                color=color, lw=lw)
    return p

def plot_ratio_e_out_qin(mtot=24.0*u.Msun, m3=1.0*u.Msun, rel_incl=0, f_dm=0,ax=None, fig_name=None, N=2):
    if not ax:
        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(100, 100)
        ax = fig.add_subplot(gs[:, :])
    G_cgs = G.to(u.cm**3/(u.g*u.s**2))
    outer_ecc = np.linspace(0, 0.99999, N)
    q_inner = np.linspace(0.1,1, N)
    # initialize ratios to 0
    v_esc_in_div_v_orb3 = np.zeros((N, N))
    for e_out in outer_ecc:
        i = np.argmin(np.absolute(outer_ecc - e_out))
        for q_in in q_inner:
            j = np.argmin(np.absolute(q_inner-q_in))
            # fill ratio
            critical_ratio = critical_a_out_div_a_in_ratio(e_out, m3/(mtot/(1-f_dm)), rel_incl)
            min_a_in = get_ZAMS_periastron_RLOF_min_a(q_in, mtot)
            a_out = critical_ratio*min_a_in.to(u.cm)
            v_orb_3 = np.sqrt(G_cgs*((mtot/(1-f_dm)+m3)).to(u.g)/a_out)
            v_esc_in = get_max_vesc(q_in, mtot=mtot.value)
            v_esc_in_div_v_orb3[j,i] = v_esc_in.to(u.km/u.s)/v_orb_3.to(u.km/u.s)
    norm = mpl.colors.Normalize(vmin=0.1, vmax=30)
    p = ax.pcolormesh(outer_ecc, q_inner, v_esc_in_div_v_orb3, norm=norm, shading='auto')
    ax.text(0.02, 0.02, r"$i_\mathrm{rel}="+f"{rel_incl}"+r"$"+"\n"+r"$f_\Delta ="+f"{f_dm:.1f}"+r"$",
            fontsize=25, va="bottom", ha="left", zorder=10, transform=ax.transAxes,
            bbox=dict(color='white', edgecolor='w', alpha=0.5,  boxstyle='round,pad=0.05'))
    contour = ax.contour(outer_ecc, q_inner, v_esc_in_div_v_orb3, levels=[5, 10], colors=["orange", "r"])
    ax.clabel(contour, inline=1, fontsize=25, fmt="%.0f", inline_spacing=1,)
    return p

def plot_ratio_e_out_rel_incl(mtot=24.0*u.Msun, m3=1.0*u.Msun, f_dm=0,
               q_in=0.1, ax=None, fig_name=None, N=2):
    if not ax:
        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(100, 100)
        ax = fig.add_subplot(gs[:, :])
    G_cgs = G.to(u.cm**3/(u.g*u.s**2))
    outer_ecc = np.linspace(0, 0.9999, N)
    outer_inc = np.linspace(0, pi, N)
    v_esc_in_div_v_orb3 = np.zeros((N, N))
    for e_out in outer_ecc:
        i = np.argmin(np.absolute(outer_ecc - e_out))
        for rel_incl in outer_inc:
            j = np.argmin(np.absolute(outer_inc-rel_incl))
            # fill ratio
            critical_ratio = critical_a_out_div_a_in_ratio(e_out, m3/(mtot/(1-f_dm)), rel_incl)
            min_a_in = get_ZAMS_periastron_RLOF_min_a(q_in, mtot)
            a_out = critical_ratio*min_a_in.to(u.cm)
            v_orb_3 = np.sqrt(G_cgs*((mtot/(1-f_dm)+m3)).to(u.g)/a_out)
            v_esc_in = get_max_vesc(q_in, mtot=mtot.value)
            v_esc_in_div_v_orb3[j,i] = v_esc_in.to(u.km/u.s)/v_orb_3.to(u.km/u.s)
    # define norm
    norm = mpl.colors.Normalize(vmin=0.1, vmax=30)
    p = ax.pcolormesh(outer_ecc, outer_inc, v_esc_in_div_v_orb3, norm=norm, shading='auto')
    ax.text(0.02, 0.02, r"$q_\mathrm{in}="+f"{q_in:.1f}"+r"$"+"\n"+r"$f_\Delta ="+f"{f_dm:.1f}"+r"$",
            fontsize=25, va="bottom", ha="left", zorder=10, transform=ax.transAxes,
            bbox=dict(color='white', edgecolor='w', alpha=0.5,  boxstyle='round,pad=0.05'))
    contour = ax.contour(outer_ecc, outer_inc, v_esc_in_div_v_orb3, levels=[5, 10], colors=["orange", "r"])
    ax.clabel(contour, inline=1, fontsize=25, fmt="%.0f", inline_spacing=1,
              manual=[(0.6, 3*pi/4),(0.9, pi/2)])
    return p



if __name__== "__main__":
    # initialize binaries
    M1 = 24
    M2 = 1
    P = 7
    HD46485 = binary(M1, M2, P)
    # -------------------------------
    M1 = 15
    M2 = 1.5
    P = 3.5
    HD191495 = binary(M1, M2, P)
    # -------------------------------
    M1 = 7
    M2 = 1
    P = 5
    HD25631 = binary(M1, M2, P)
    # -------------------------------
    # -------------------------------
    # Plot for minimum inner separation
    plot_min_a_in(HD46485, fig_name="a_min_RL_24_1_7.pdf", title=r"HD46485") # GOOD
    # plot_min_a_out(HD45485,  title=r"HD45485")
    # # print("24+1, 7d", HD45485["a"])
    # # -------------------------------
    plot_min_a_in(HD191495, fig_name="a_min_RL_15_1p5_3p5.pdf", title=r"HD191495") #GOOD
    # print("15+1.5, 3.5d", HD191495["a"])
    # # -------------------------------
    plot_min_a_in(HD25631, fig_name="a_min_RL_7_1_5.pdf", title=r"HD25631") # GOOD
    # print("7+1, 5d", HD25631["a"])
    # # -------------------------------
    # # -------------------------------
