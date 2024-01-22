import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import os
import mesaPlot as mp
from showyourwork.paths import user as Paths

paths = Paths()
if os.path.exists(os.path.join(paths.data, 'post_interaction/30_20_10_res/LOGS1/history.data')):
    pass
else:
    os.system(f'python {os.path.join(paths.scripts / "unzip_MESA_output.py")}')

plt.style.use(paths.scripts / "matplotlibrc")


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


m_200 = mp.MESA()
m2_200 = mp.MESA()
m3_200 = mp.MESA()

m_200_newtides = mp.MESA()
m2_200_newtides = mp.MESA()
m3_200_newtides = mp.MESA()

m_200_res = mp.MESA()
m2_200_res = mp.MESA()
m3_200_res = mp.MESA()

m_200_newtides_res = mp.MESA()
m2_200_newtides_res = mp.MESA()
m3_200_newtides_res = mp.MESA()

name = 'post_interaction/30_20_10'

m_200_res.log_fold = os.path.join(paths.data, name + '_res/LOGS1')
m_200_res.loadHistory()

m2_200_res.log_fold = os.path.join(paths.data, name + '_res/LOGS2')
m2_200_res.loadHistory()

m3_200_res.log_fold = os.path.join(paths.data, name + '_res/LOGS3')
m3_200_res.loadHistory()

m_200_newtides_res.log_fold = os.path.join(paths.data, name + '_posydon_res/LOGS1')
m_200_newtides_res.loadHistory()

m2_200_newtides_res.log_fold = os.path.join(paths.data, name + '_posydon_res/LOGS2')
m2_200_newtides_res.loadHistory()

m3_200_newtides_res.log_fold = os.path.join(paths.data, name + '_posydon_res/LOGS3')
m3_200_newtides_res.loadHistory()

m_200.log_fold = os.path.join(paths.data, name + '/LOGS1')
m_200.loadHistory()

m2_200.log_fold = os.path.join(paths.data, name + '/LOGS2')
m2_200.loadHistory()

m3_200.log_fold = os.path.join(paths.data, name + '/LOGS3')
m3_200.loadHistory()

m_200_newtides.log_fold = os.path.join(paths.data, name + '_posydon/LOGS1')
m_200_newtides.loadHistory()

m2_200_newtides.log_fold = os.path.join(paths.data, name + '_posydon/LOGS2')
m2_200_newtides.loadHistory()

m3_200_newtides.log_fold = os.path.join(paths.data, name + '_posydon/LOGS3')
m3_200_newtides.loadHistory()

star_age_200 = m_200.hist.star_age

surf_avg_vtor_1 = m_200.hist.surf_avg_v_rot
surf_avg_vtor_2 = m2_200.hist.surf_avg_v_rot

surf_avg_omega200 = m_200.hist.surf_avg_omega
star_1_radius200 = m3_200.hist.star_1_radius
star_1_J_orb_200 = m3_200.hist.J_orb
star_1_J_spin_200 = m3_200.hist.J_spin_1
star_2_J_spin_200 = m3_200.hist.J_spin_2
star_1_mass = m3_200.hist.star_1_mass
star_2_mass = m3_200.hist.star_2_mass

star_1_mass_res = m3_200_res.hist.star_1_mass
star_2_mass_res = m3_200_res.hist.star_2_mass

age_200 = m3_200.hist.age
jtotal_200 = m3_200.hist.J_total
log_total_angular_momentum_200 = m_200.hist.log_total_angular_momentum
surf_avg_j_rot_200 = m_200.hist.surf_avg_j_rot
center_h1_200 = m_200.hist.center_h1

LOGL_1 = m_200.hist.log_L
LOGL_2 = m2_200.hist.log_L
LOGL_1_newtides = m_200_newtides.hist.log_L
LOGL_2_newtides = m2_200_newtides.hist.log_L

log_Teff_1 = m_200.hist.log_Teff
log_Teff_2 = m2_200.hist.log_Teff

log_Teff_1_newtides = m_200_newtides.hist.log_Teff
log_Teff_2_newtides = m2_200_newtides.hist.log_Teff

star_age_200_newtides = m_200_newtides.hist.star_age

surf_avg_vtor_1_newtides = m_200_newtides.hist.surf_avg_v_rot
surf_avg_vtor_2_newtides = m2_200_newtides.hist.surf_avg_v_rot

surf_avg_omega200_newtides = m_200_newtides.hist.surf_avg_omega
star_1_radius200_newtides = m3_200_newtides.hist.star_1_radius
star_1_J_orb_200_newtides = m3_200_newtides.hist.J_orb
star_1_J_spin_200_newtides = m3_200_newtides.hist.J_spin_1
star_2_J_spin_200_newtides = m3_200_newtides.hist.J_spin_2
age_200_newtides = m3_200_newtides.hist.age
jtotal_200_newtides = m3_200_newtides.hist.J_total
log_total_angular_momentum_200_newtides = m_200_newtides.hist.log_total_angular_momentum
surf_avg_j_rot_200_newtides = m_200_newtides.hist.surf_avg_j_rot
center_h1_200_newtides = m_200_newtides.hist.center_h1

period_class = m3_200.hist.period_days
period_posydon = m3_200_newtides.hist.period_days

period_class_res = m3_200_res.hist.period_days
period_posydon_res = m3_200_newtides_res.hist.period_days

star_1_radius_class = m3_200.hist.star_1_radius
star_1_radius_posydon = m3_200_newtides.hist.star_1_radius

star_2_radius_class = m3_200.hist.star_2_radius
star_2_radius_posydon = m3_200_newtides.hist.star_2_radius

J_orb_class = m3_200.hist.J_orb
J_orb_posydon = m3_200_newtides.hist.J_orb

J_spin2_class = m3_200.hist.J_spin_2
J_spin2_posydon = m3_200_newtides.hist.J_spin_2

J_spin1_class = m3_200.hist.J_spin_1
J_spin1_posydon = m3_200_newtides.hist.J_spin_1
star_1_mass_posydon = m3_200_newtides.hist.star_1_mass
star_2_mass_posydon = m3_200_newtides.hist.star_2_mass

star_1_mass_posydon_res = m3_200_newtides_res.hist.star_1_mass
star_2_mass_posydon_res = m3_200_newtides_res.hist.star_2_mass

surf_avg_omega_1_class = m_200.hist.surf_avg_omega

surf_avg_omega_1_pos = m_200_newtides.hist.surf_avg_omega

surf_avg_omega_2_class = m2_200.hist.surf_avg_omega

surf_avg_omega_2_pos = m2_200_newtides.hist.surf_avg_omega

star_age_pos = m2_200_newtides.hist.star_age
star_age_class = m2_200.hist.star_age

# p1.log_fold='LOGS1'

# p1.loadProfile(num=-1)

# p=mp.plot()


rl_relative_gap_1 = m3_200.hist.rl_relative_overflow_1
rl_relative_gap_1_posydon = m3_200_newtides.hist.rl_relative_overflow_1

rl_relative_gap_1_res = m3_200_res.hist.rl_relative_overflow_1
rl_relative_gap_1_posydon_res = m3_200_newtides_res.hist.rl_relative_overflow_1

age_class = m3_200.hist.age
age_posydon = m3_200_newtides.hist.age

lg_t_sync_2_class = m3_200.hist.lg_t_sync_2
lg_t_sync_2_posydon = m3_200_newtides.hist.lg_t_sync_2

lg_t_sync_1_class = m3_200.hist.lg_t_sync_1
lg_t_sync_1_posydon = m3_200_newtides.hist.lg_t_sync_1

iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_1_posydon = rl_relative_gap_1_posydon > 0

iRLOF_1_res = rl_relative_gap_1_res > 0
iRLOF_1_posydon_res = rl_relative_gap_1_posydon_res > 0

age_class_rlof = age_class[iRLOF_1]
age_posydon_rlof = age_posydon[iRLOF_1_posydon]

J_orb_class_rlof = J_orb_class[iRLOF_1]
J_orb_posydon_rlof = J_orb_posydon[iRLOF_1_posydon]

i_pre_RLOF_class = age_class < min(age_class[iRLOF_1])
i_pre_RLOF_pos = age_posydon < min(age_posydon[iRLOF_1_posydon])

i_post_RLOF_class = age_class > max(age_class[iRLOF_1])
i_post_RLOF_pos = age_posydon > max(age_posydon[iRLOF_1_posydon])

star_age_rlof_ind = find_nearest(star_age_class, min(age_class[iRLOF_1]))
star_age_rlof_ind_pos = find_nearest(star_age_pos, min(age_posydon[iRLOF_1_posydon]))

pp1_all_panel = PdfPages(paths.figures / 'hr_panel_res.pdf')

fig = plt.figure(figsize=(10, 10))

plt.title('$\it{M}_\mathrm{don,ini}$ = 30, $\it{M}_\mathrm{acc,ini}$ = 20, $\it{P}_\mathrm{ini}$ = 10 d', fontsize=30)
plt.plot(star_2_mass / star_1_mass, period_class, color='k', linestyle='-',
         label='MESA default, mesh delta=1, time delta=1', lw=15)
plt.plot(star_2_mass[iRLOF_1] / star_1_mass[iRLOF_1], period_class[iRLOF_1], lw=20, c='k')

plt.plot(star_2_mass_res / star_1_mass_res, period_class_res, color='red', linestyle='--',
         label='MESA default, mesh delta=0.5, time delta=0.75', lw=2)
plt.plot(star_2_mass_res[iRLOF_1_res] / star_1_mass_res[iRLOF_1_res], period_class_res[iRLOF_1_res], lw=4, c='red',
         linestyle='--')

plt.plot(star_2_mass_posydon / star_1_mass_posydon, period_posydon, linestyle='-',
         label='POSYDON, mesh delta=1, time delta=1', color='green', lw=12)
plt.plot(star_2_mass_posydon[iRLOF_1_posydon] / star_1_mass_posydon[iRLOF_1_posydon], period_posydon[iRLOF_1_posydon],
         lw=17, c='green')

plt.plot(star_2_mass_posydon_res / star_1_mass_posydon_res, period_posydon_res, linestyle='--',
         label='POSYDON, mesh delta=0.5, time delta=0.75', color='magenta')
plt.plot(star_2_mass_posydon_res[iRLOF_1_posydon_res] / star_1_mass_posydon_res[iRLOF_1_posydon_res],
         period_posydon_res[iRLOF_1_posydon_res], lw=5, c='magenta', linestyle='--')

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.xlim([5e6,7.5e6])
# plt.ylim([3,17])
plt.xlabel('$\it{Q}$', fontsize=30)
plt.ylabel('Period [days]', fontsize=30)
plt.legend(loc=2, fontsize=18)

plt.xlim([0, 2.5])
###plt.xlim([0.3,1.5])
##plt.xlim([0.64,0.72])
# plt.xlim([0.34,2]) # 5 days
# plt.xlim([0.5,3]) # 3 days

plt.savefig(pp1_all_panel, format='pdf')

pp1_all_panel.close()
