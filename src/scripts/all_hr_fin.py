import os

import matplotlib.pyplot as plt
import mesaPlot as mp
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from showyourwork.paths import user as Paths

paths = Paths()
plt.style.use(paths.scripts / "matplotlibrc")

if os.path.exists(os.path.join(paths.data, 'HR_HD25631/LOGS1_0/history.data')):
    pass
else:
    os.system(f'python {os.path.join(paths.scripts / "unzip_MESA_output.py")}')


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


m = mp.MESA()
m3 = mp.MESA()

m_50 = mp.MESA()
m3_50 = mp.MESA()

m_100 = mp.MESA()
m3_100 = mp.MESA()

m_150 = mp.MESA()
m3_150 = mp.MESA()

m_200 = mp.MESA()
m3_200 = mp.MESA()

m_230 = mp.MESA()
m3_230 = mp.MESA()

m3_no_wind_no_tides = mp.MESA()
m3_no_tides = mp.MESA()
m_no_tides = mp.MESA()

m3_sync = mp.MESA()
m_sync = mp.MESA()

m.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_0')
m.loadHistory()

m3.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_0')
m3.loadHistory()

m3_sync.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_sync')
m3_sync.loadHistory()
m_sync.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_sync')
m_sync.loadHistory()

m3_no_tides.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_no_tides')
m3_no_tides.loadHistory()
m_no_tides.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_no_tides')
m_no_tides.loadHistory()

# m3_no_wind_no_tides.log_fold='LOGS3_no_wind_no_tides'
# m3_no_wind_no_tides.loadHistory()


m_50.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_50')
m_50.loadHistory()

m3_50.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_50')
m3_50.loadHistory()

m_100.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_100')
m_100.loadHistory()

m3_100.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_100')
m3_100.loadHistory()

m_150.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_150')
m_150.loadHistory()

m3_150.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_150')
m3_150.loadHistory()

m_200.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_200')
m_200.loadHistory()

m3_200.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_200')
m3_200.loadHistory()

m_230.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS1_230')
m_230.loadHistory()

m3_230.log_fold = os.path.join(paths.data, 'HR_HD25631/LOGS3_230')
m3_230.loadHistory()

# age_1_nw_nt=m3_no_wind_no_tides.hist.age
# period_1_nw_nt=m3_no_wind_no_tides.hist.period_days


age_1_nt = m3_no_tides.hist.age
period_1_nt = m3_no_tides.hist.period_days
star_age_1_nt = m_no_tides.hist.star_age

LOGL_1_1_nt = m_no_tides.hist.log_L
log_Teff_1_1_nt = m_no_tides.hist.log_Teff
rl_relative_gap_1_1_nt = m3_no_tides.hist.rl_relative_overflow_1
iRLOF_1_1_nt = rl_relative_gap_1_1_nt > 0
star_age_rlof_ind_1_nt = find_nearest(star_age_1_nt, min(age_1_nt[iRLOF_1_1_nt]))

age_1_s = m3_sync.hist.age
period_1_s = m3_sync.hist.period_days
star_age_1_s = m_sync.hist.star_age

LOGL_1_1_s = m_sync.hist.log_L
log_Teff_1_1_s = m_sync.hist.log_Teff
rl_relative_gap_1_1_s = m3_sync.hist.rl_relative_overflow_1
iRLOF_1_1_s = rl_relative_gap_1_1_s > 0
star_age_rlof_ind_1_s = find_nearest(star_age_1_s, min(age_1_s[iRLOF_1_1_s]))

star_age_1 = m.hist.star_age

surf_avg_vtor_1 = m.hist.surf_avg_v_rot
surf_avg_omega1 = m.hist.surf_avg_omega
star_1_radius1 = m3.hist.star_1_radius
star_1_J_orb_1 = m3.hist.J_orb
star_1_J_spin_1 = m3.hist.J_spin_1
star_2_J_spin_1 = m3.hist.J_spin_2
age_1 = m3.hist.age
jtotal_1 = m3.hist.J_total
log_total_angular_momentum_1 = m.hist.log_total_angular_momentum
surf_avg_j_rot_1 = m.hist.surf_avg_j_rot
center_h1_1 = m.hist.center_h1
indMS_1 = center_h1_1 > 1e-2
LOGL_1_1 = m.hist.log_L
log_Teff_1_1 = m.hist.log_Teff
rl_relative_gap_1_1 = m3.hist.rl_relative_overflow_1
iRLOF_1_1 = rl_relative_gap_1_1 > 0
star_age_rlof_ind_1 = find_nearest(star_age_1, min(age_1[iRLOF_1_1]))
star_age_ms_ind_1 = find_nearest(age_1, star_age_1[indMS_1][-1])
period_1 = m3.hist.period_days

star_age_50 = m_50.hist.star_age

surf_avg_vtor_50 = m_50.hist.surf_avg_v_rot
surf_avg_omega50 = m_50.hist.surf_avg_omega
star_1_radius50 = m3_50.hist.star_1_radius
star_1_J_orb_50 = m3_50.hist.J_orb
star_1_J_spin_50 = m3_50.hist.J_spin_1
star_2_J_spin_50 = m3_50.hist.J_spin_2
age_50 = m3_50.hist.age
jtotal_50 = m3_50.hist.J_total
log_total_angular_momentum_50 = m_50.hist.log_total_angular_momentum
surf_avg_j_rot_50 = m_50.hist.surf_avg_j_rot
center_h1_50 = m_50.hist.center_h1
indMS_50 = center_h1_50 > 1e-2
period_50 = m3_50.hist.period_days

LOGL_1_50 = m_50.hist.log_L
log_Teff_1_50 = m_50.hist.log_Teff
rl_relative_gap_1_50 = m3_50.hist.rl_relative_overflow_1
iRLOF_1_50 = rl_relative_gap_1_50 > 0
star_age_rlof_ind_50 = find_nearest(star_age_50, min(age_50[iRLOF_1_50]))
star_age_ms_ind_50 = find_nearest(age_50, star_age_50[indMS_50][-1])
star_age_100 = m_100.hist.star_age

surf_avg_vtor_100 = m_100.hist.surf_avg_v_rot
surf_avg_omega100 = m_100.hist.surf_avg_omega
star_1_radius100 = m3_100.hist.star_1_radius
star_1_J_orb_100 = m3_100.hist.J_orb
star_1_J_spin_100 = m3_100.hist.J_spin_1
star_2_J_spin_100 = m3_100.hist.J_spin_2
age_100 = m3_100.hist.age
jtotal_100 = m3_100.hist.J_total
log_total_angular_momentum_100 = m_100.hist.log_total_angular_momentum
surf_avg_j_rot_100 = m_100.hist.surf_avg_j_rot
center_h1_100 = m_100.hist.center_h1
indMS_100 = center_h1_100 > 1e-2

period_100 = m3_100.hist.period_days

LOGL_1_100 = m_100.hist.log_L
log_Teff_1_100 = m_100.hist.log_Teff
rl_relative_gap_1_100 = m3_100.hist.rl_relative_overflow_1
iRLOF_1_100 = rl_relative_gap_1_100 > 0
star_age_rlof_ind_100 = find_nearest(star_age_100, min(age_100[iRLOF_1_100]))
star_age_ms_ind_100 = find_nearest(age_100, star_age_100[indMS_100][-1])

star_age_150 = m_150.hist.star_age

surf_avg_vtor_150 = m_150.hist.surf_avg_v_rot
surf_avg_omega150 = m_150.hist.surf_avg_omega
star_1_radius150 = m3_150.hist.star_1_radius
star_1_J_orb_150 = m3_150.hist.J_orb
star_1_J_spin_150 = m3_150.hist.J_spin_1
star_2_J_spin_150 = m3_150.hist.J_spin_2
age_150 = m3_150.hist.age
jtotal_150 = m3_150.hist.J_total
log_total_angular_momentum_150 = m_150.hist.log_total_angular_momentum
surf_avg_j_rot_150 = m_150.hist.surf_avg_j_rot
center_h1_150 = m_150.hist.center_h1
indMS_150 = center_h1_150 > 1e-2
star_age_ms_ind_150 = find_nearest(age_150, star_age_150[indMS_150][-1])
period_150 = m3_150.hist.period_days

LOGL_1_150 = m_150.hist.log_L
log_Teff_1_150 = m_150.hist.log_Teff
rl_relative_gap_1_150 = m3_150.hist.rl_relative_overflow_1
iRLOF_1_150 = rl_relative_gap_1_150 > 0
star_age_rlof_ind_150 = find_nearest(star_age_150, min(age_150[iRLOF_1_150]))

star_age_200 = m_200.hist.star_age

surf_avg_vtor_200 = m_200.hist.surf_avg_v_rot
surf_avg_omega200 = m_200.hist.surf_avg_omega
star_1_radius200 = m3_200.hist.star_1_radius
star_1_J_orb_200 = m3_200.hist.J_orb
star_1_J_spin_200 = m3_200.hist.J_spin_1
star_2_J_spin_200 = m3_200.hist.J_spin_2
age_200 = m3_200.hist.age
jtotal_200 = m3_200.hist.J_total
log_total_angular_momentum_200 = m_200.hist.log_total_angular_momentum
surf_avg_j_rot_200 = m_200.hist.surf_avg_j_rot
center_h1_200 = m_200.hist.center_h1
period_200 = m3_200.hist.period_days
indMS_200 = center_h1_200 > 1e-2

LOGL_1_200 = m_200.hist.log_L
log_Teff_1_200 = m_200.hist.log_Teff
rl_relative_gap_1_200 = m3_200.hist.rl_relative_overflow_1
iRLOF_1_200 = rl_relative_gap_1_200 > 0
star_age_rlof_ind_200 = find_nearest(star_age_200, min(age_200[iRLOF_1_200]))
star_age_ms_ind_200 = find_nearest(age_200, star_age_200[indMS_200][-1])

star_age_230 = m_230.hist.star_age
surf_avg_vtor_230 = m_230.hist.surf_avg_v_rot
age_230 = m3_230.hist.age

center_h1_230 = m_230.hist.center_h1
indMS_230 = center_h1_230 > 1e-2
star_age_ms_ind_230 = find_nearest(age_230, star_age_230[indMS_230][-1])

LOGL_1_230 = m_230.hist.log_L
log_Teff_1_230 = m_230.hist.log_Teff

rl_relative_gap_1_230 = m3_230.hist.rl_relative_overflow_1
iRLOF_1_230 = rl_relative_gap_1_230 > 0
star_age_rlof_ind_230 = find_nearest(star_age_230, min(age_230[iRLOF_1_230]))
period_230 = m3_230.hist.period_days

colors_list = plt.cm.viridis(np.linspace(0, 1, 7))

pp1_all_t = PdfPages(paths.figures / 'age_vsurf_hr_all.pdf')

i_pre_RLOF_1 = age_1 < min(age_1[iRLOF_1_1])
i_pre_RLOF_1_nt = age_1_nt < min(age_1_nt[iRLOF_1_1_nt])
i_pre_RLOF_1_s = age_1_s < min(age_1_s[iRLOF_1_1_s])
i_pre_RLOF_50 = age_50 < min(age_50[iRLOF_1_50])
i_pre_RLOF_100 = age_100 < min(age_100[iRLOF_1_100])
i_pre_RLOF_150 = age_150 < min(age_150[iRLOF_1_150])

i_pre_RLOF_200 = age_200 < min(age_200[iRLOF_1_200])
i_pre_RLOF_230 = age_230 < min(age_230[iRLOF_1_230])

fig, ax1 = plt.subplots(figsize=(10, 10))

# ax1.set_title('M1$_{ini}$=7 M$_{\odot}$, M2$_{ini}$=1 M$_{\odot}$, p$_{ini}$=5.2 d',fontsize=30)
ax1.set_title('$\it{M}_\mathrm{1,ini}$=7, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=5.2 d', fontsize=30)

ax1.plot(age_230[i_pre_RLOF_230] / 1e6, period_230[i_pre_RLOF_230], linestyle='-',
         label='$\it{v}^\mathrm{1,ini}_\mathrm{rot}$ = 230 $km~s^{-1}$', lw=2, c=colors_list[0])
ax1.scatter(age_230[iRLOF_1_230][0] / 1e6, period_230[iRLOF_1_230][0], marker='x', c='k', s=100)

ax1.plot(age_200[i_pre_RLOF_200] / 1e6, period_200[i_pre_RLOF_200], linestyle='-',
         label='$\it{v}^\mathrm{1,ini}_\mathrm{rot}$ = 200 $km~s^{-1}$', lw=2, c=colors_list[1])
ax1.scatter(age_200[iRLOF_1_200][0] / 1e6, period_200[iRLOF_1_200][0], marker='x', c='k', s=100)
# ax1.scatter(star_age_200[indMS_200][-1]/1e6, period_200[star_age_ms_ind_200], marker='^',  c='k',s=100)


ax1.plot(age_150[i_pre_RLOF_150] / 1e6, period_150[i_pre_RLOF_150], linestyle='-',
         label='$\it{v}^\mathrm{1,ini}_\mathrm{rot}$= 150 $km~s^{-1}$', lw=2, c=colors_list[2])
ax1.scatter(age_150[iRLOF_1_150][0] / 1e6, period_150[iRLOF_1_150][0], marker='x', c='k', s=100)
# ax1.scatter(star_age_150[indMS_150][-1]/1e6, period_150[star_age_ms_ind_150], marker='^',  c='k',s=100)


ax1.plot(age_100[i_pre_RLOF_100] / 1e6, period_100[i_pre_RLOF_100], linestyle='-',
         label='$\it{v}^\mathrm{1,ini}_\mathrm{rot}$ = 100 $km~s^{-1}$', lw=2, c=colors_list[3])
ax1.scatter(age_100[iRLOF_1_100][0] / 1e6, period_100[iRLOF_1_100][0], marker='x', c='k', s=100)
# ax1.scatter(star_age_100[indMS_100][-1]/1e6, period_100[star_age_ms_ind_100], marker='^',  c='k',s=100)

ax1.plot(age_50[i_pre_RLOF_50] / 1e6, period_50[i_pre_RLOF_50], linestyle='-',
         label='$\it{v}^\mathrm{1,ini}_\mathrm{rot}$ = 50 $km~s^{-1}$', lw=2, c=colors_list[4])
ax1.scatter(age_50[iRLOF_1_50][0] / 1e6, period_50[iRLOF_1_50][0], marker='o', c='k', s=100)
# ax1.scatter(star_age_50[indMS_50][-1]/1e6, period_50[star_age_ms_ind_50], marker='^',  c='k',s=100)


ax1.plot(age_1[i_pre_RLOF_1] / 1e6, period_1[i_pre_RLOF_1], linestyle='-',
         label='$\it{v}^\mathrm{1,ini}_\mathrm{rot}$ = 1 $km~s^{-1}$', lw=2, c=colors_list[5])
ax1.scatter(age_1[iRLOF_1_1][0] / 1e6, period_1[iRLOF_1_1][0], marker='o', c='k', s=100)

ax1.plot(age_1_nt[i_pre_RLOF_1_nt] / 1e6, period_1_nt[i_pre_RLOF_1_nt], linestyle='--',
         label='$\it{v}^\mathrm{1,ini}_\mathrm{rot}$ = 1 $km~s^{-1}$, no tides', lw=2, color='orange')
ax1.scatter(age_1_nt[iRLOF_1_1_nt][0] / 1e6, period_1_nt[iRLOF_1_1_nt][0], marker='x', c='orange', s=100)

ax1.plot(age_1_s[i_pre_RLOF_1_s] / 1e6, period_1_s[i_pre_RLOF_1_s], linestyle='--', label='initially synchronized',
         lw=2, color='red')
ax1.scatter(age_1_s[iRLOF_1_1_s][0] / 1e6, period_1_s[iRLOF_1_1_s][0], marker='o', c='red', s=100)

# ax1.scatter(star_age_1[indMS_1][-1]/1e6, period_1[star_age_ms_ind_1], marker='^',  c='k',s=100)


plt.plot([20, 80], [5.2, 5.2], color='red', linestyle='--', lw=0.5)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

ax1.text(0.88, 0.90, 'case B', transform=ax1.transAxes, fontsize=15,
         verticalalignment='top', bbox=props)

ax1.text(0.88, 0.80, 'case B', transform=ax1.transAxes, fontsize=15,
         verticalalignment='top', bbox=props)

ax1.text(0.88, 0.60, 'case B', transform=ax1.transAxes, fontsize=15,
         verticalalignment='top', bbox=props)

ax1.text(0.88, 0.43, 'case B', transform=ax1.transAxes, fontsize=15,
         verticalalignment='top', bbox=props)

ax1.text(0.88, 0.21, 'case A', transform=ax1.transAxes, fontsize=15,
         verticalalignment='top', bbox=props)

ax1.text(0.88, 0.06, 'case A', transform=ax1.transAxes, fontsize=15,
         verticalalignment='top', bbox=props)

# ax1.scatter(star_age_230[indMS_230][-1]/1e6, period_230[star_age_ms_ind_230], marker='^',  c='k',s=100)


ax1.plot([star_age_230[indMS_230][-1] / 1e6, star_age_200[indMS_200][-1] / 1e6],
         [period_230[star_age_ms_ind_230], period_200[star_age_ms_ind_200]], linestyle='--', color='gray')

ax1.plot([star_age_200[indMS_200][-1] / 1e6, star_age_150[indMS_150][-1] / 1e6],
         [period_200[star_age_ms_ind_200], period_150[star_age_ms_ind_150]], linestyle='--', color='gray')

ax1.plot([star_age_150[indMS_150][-1] / 1e6, star_age_100[indMS_100][-1] / 1e6],
         [period_150[star_age_ms_ind_150], period_100[star_age_ms_ind_100]], linestyle='--', color='gray')

ax1.plot([star_age_100[indMS_100][-1] / 1e6, star_age_50[indMS_50][-1] / 1e6],
         [period_100[star_age_ms_ind_100], period_50[star_age_ms_ind_50]], linestyle='--', color='gray')

plt.xlim([2e7 / 1e6, 6.3e7 / 1e6])
plt.ylim([4.61, 5.75])

# plt.plot([min(star_age_2),max(star_age_2)],[201,201],color='r',linestyle='--',label='HD 191495 m1=15, m2=1.5 p=3.64d')
# plt.plot([min(star_age_3),max(star_age_3)],[334,334],color='m',linestyle='--',label='HD 46485 m1=24  m2=1 p=6.9d')

# plt.gca().invert_xaxis()
plt.ylabel('Period [days]', fontsize=30)
plt.xlabel('Star age [Myrs]', fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(loc=2, fontsize=16)

left, bottom, width, height = [0.25, 0.2, 0.3, 0.3]
ax3 = fig.add_axes([left, bottom, width, height])

ax3.plot(log_Teff_1_1, LOGL_1_1, linestyle='-', label='vini = 1 [km/s]', lw=2, c=colors_list[5])
ax3.plot(log_Teff_1_1[star_age_rlof_ind_1:-1], LOGL_1_1[star_age_rlof_ind_1:-1], lw=10, c='k', zorder=0)

ax3.plot(log_Teff_1_1_nt, LOGL_1_1_nt, label='vini = 1 [km/s]', lw=2, c='orange', linestyle='--')
ax3.plot(log_Teff_1_1_nt[star_age_rlof_ind_1_nt:-1], LOGL_1_1_nt[star_age_rlof_ind_1_nt:-1], lw=10, c='k', zorder=0)

ax3.plot(log_Teff_1_1_s, LOGL_1_1_s, label='vini = 1 [km/s]', lw=2, c='red', linestyle='--')
ax3.plot(log_Teff_1_1_s[star_age_rlof_ind_1_s:-1], LOGL_1_1_s[star_age_rlof_ind_1_s:-1], lw=10, c='k', zorder=0)

ax3.plot(log_Teff_1_50, LOGL_1_50, linestyle='-', label='vini = 50 [km/s]', lw=2, c=colors_list[4])
ax3.plot(log_Teff_1_50[star_age_rlof_ind_50:-1], LOGL_1_50[star_age_rlof_ind_50:-1], lw=10, c='k', zorder=0)

ax3.plot(log_Teff_1_100, LOGL_1_100, linestyle='-', label='vini = 100 [km/s]', lw=2, c=colors_list[3])
ax3.plot(log_Teff_1_100[star_age_rlof_ind_100:-1], LOGL_1_100[star_age_rlof_ind_100:-1], lw=10, c='k', zorder=0)

ax3.plot(log_Teff_1_150, LOGL_1_150, linestyle='-', label='vini = 150 [km/s]', lw=2, c=colors_list[2])
ax3.plot(log_Teff_1_150[star_age_rlof_ind_150:-1], LOGL_1_150[star_age_rlof_ind_150:-1], lw=10, c='k', zorder=0)

ax3.plot(log_Teff_1_200, LOGL_1_200, linestyle='-', label='vini = 200 [km/s]', lw=2, c=colors_list[1])
ax3.plot(log_Teff_1_200[star_age_rlof_ind_200:-1], LOGL_1_200[star_age_rlof_ind_200:-1], lw=10, c='k', zorder=0)

ax3.plot(log_Teff_1_230, LOGL_1_230, linestyle='-', label='vini = 230 [km/s]', lw=2, c=colors_list[0])
ax3.plot(log_Teff_1_230[star_age_rlof_ind_230:-1], LOGL_1_230[star_age_rlof_ind_230:-1], lw=10, c='k', zorder=0)
# ax3.scatter(log_Teff_1_230[indMS_230][-1], LOGL_1_230[indMS_230][-1], marker='^',  c='k',s=200)
# ax3.scatter(log_Teff_1_200[indMS_200][-1], LOGL_1_200[indMS_200][-1], marker='^',  c='k',s=200)
# ax3.scatter(log_Teff_1_150[indMS_150][-1], LOGL_1_150[indMS_150][-1], marker='^',  c='k',s=200)
# ax3.scatter(log_Teff_1_100[indMS_100][-1], LOGL_1_100[indMS_100][-1], marker='^',  c='k',s=200)
# ax3.scatter(log_Teff_1_50[indMS_50][-1], LOGL_1_50[indMS_50][-1], marker='^',  c='k',s=200)
# ax3.scatter(log_Teff_1_1[indMS_1][-1], LOGL_1_1[indMS_1][-1], marker='^',  c='k',s=200)


ax3.tick_params(labelsize=13)

ax3.set_ylim([3.4, 4])
ax3.set_xlim([4.1, 4.2])
ax3.set_ylabel('log$_{10} (L/L_{\odot}$)', fontsize=18)
ax3.set_xlabel('log$_{10} (T_{eff}/K$)', fontsize=18)
ax3.invert_xaxis()

plt.savefig(pp1_all_t, format='pdf')
pp1_all_t.close()
