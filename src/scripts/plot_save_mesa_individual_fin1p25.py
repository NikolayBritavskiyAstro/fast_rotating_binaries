import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import mesaPlot as mp
from showyourwork.paths import user as Paths

paths = Paths()

if os.path.exists(os.path.join(paths.data, 'post_interaction/30_20_1p25_g1_new/LOGS3/history.data')):
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

m3_200_g1_new = mp.MESA()
m3_200_g2_new = mp.MESA()
m3_200_g10_new = mp.MESA()

name = 'post_interaction/30_20_1p25'
m3_200_g1_new.log_fold = os.path.join(paths.data, name + '_g1_new/LOGS3')
print(os.path.join(paths.data, name + '_g1_new/LOGS3'))
m3_200_g1_new.loadHistory()

m3_200_g2_new.log_fold = os.path.join(paths.data, name + '_g2_new/LOGS3')
# m3_200_g2_new.log_fold=name+'_g1_pos_new/LOGS3'

m3_200_g2_new.loadHistory()

m3_200_g10_new.log_fold = os.path.join(paths.data, name + '_g10_new/LOGS3')
m3_200_g10_new.loadHistory()

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

rl_relative_gap_1 = m3_200.hist.rl_relative_overflow_1
rl_relative_gap_2 = m3_200.hist.rl_relative_overflow_2
star_1_mass = m3_200.hist.star_1_mass
star_2_mass = m3_200.hist.star_2_mass
iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_2 = rl_relative_gap_2 > 0

period_class = m3_200.hist.period_days

rl_relative_gap_1_g1_new = m3_200_g1_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g1_new = m3_200_g1_new.hist.rl_relative_overflow_2

star_1_mass_g1_new = m3_200_g1_new.hist.star_1_mass
star_2_mass_g1_new = m3_200_g1_new.hist.star_2_mass
iRLOF_1_g1_new = rl_relative_gap_1_g1_new > 0
iRLOF_2_g1_new = rl_relative_gap_2_g1_new > 0
period_class_g1_new = m3_200_g1_new.hist.period_days

rl_relative_gap_1_g2_new = m3_200_g2_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g2_new = m3_200_g2_new.hist.rl_relative_overflow_2
star_1_mass_g2_new = m3_200_g2_new.hist.star_1_mass
star_2_mass_g2_new = m3_200_g2_new.hist.star_2_mass
iRLOF_1_g2_new = rl_relative_gap_1_g2_new > 0
iRLOF_2_g2_new = rl_relative_gap_2_g2_new > 0
period_class_g2_new = m3_200_g2_new.hist.period_days

rl_relative_gap_1_g10_new = m3_200_g10_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g10_new = m3_200_g10_new.hist.rl_relative_overflow_2
star_1_mass_g10_new = m3_200_g10_new.hist.star_1_mass
star_2_mass_g10_new = m3_200_g10_new.hist.star_2_mass
iRLOF_1_g10_new = rl_relative_gap_1_g10_new > 0
iRLOF_2_g10_new = rl_relative_gap_2_g10_new > 0
period_class_g10_new = m3_200_g10_new.hist.period_days

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

surf_avg_omega_1_class = m_200.hist.surf_avg_omega

surf_avg_omega_1_pos = m_200_newtides.hist.surf_avg_omega

surf_avg_omega_2_class = m2_200.hist.surf_avg_omega

surf_avg_omega_2_pos = m2_200_newtides.hist.surf_avg_omega

star_age_pos = m2_200_newtides.hist.star_age
star_age_class = m2_200.hist.star_age

# p1.log_fold='LOGS1'

# p1.loadProfile(num=-1)

# p=mp.plot()


rl_relative_gap_1_posydon = m3_200_newtides.hist.rl_relative_overflow_1
rl_relative_gap_2_posydon = m3_200_newtides.hist.rl_relative_overflow_2

age_class = m3_200.hist.age
age_posydon = m3_200_newtides.hist.age

lg_t_sync_2_class = m3_200.hist.lg_t_sync_2
lg_t_sync_2_posydon = m3_200_newtides.hist.lg_t_sync_2

lg_t_sync_1_class = m3_200.hist.lg_t_sync_1
lg_t_sync_1_posydon = m3_200_newtides.hist.lg_t_sync_1

iRLOF_1_posydon = rl_relative_gap_1_posydon > 0
iRLOF_2_posydon = rl_relative_gap_2_posydon > 0

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

pp1_all_panel = PdfPages(paths.figures / 'p_q_1p25days.pdf')

fig = plt.figure(figsize=(10, 10))

plt.title('$\it{P}_\mathrm{ini}$ = 1.25 [days]', fontsize=30)

plt.plot(star_2_mass / star_1_mass, period_class, color='k', linestyle='-', label='MESA default, $\it{\gamma}$ = 0',
         lw=2)
plt.plot(star_2_mass[iRLOF_1] / star_1_mass[iRLOF_1], period_class[iRLOF_1], lw=7, c='k')

plt.plot(star_2_mass_g1_new / star_1_mass_g1_new, period_class_g1_new, label='MESA default, $\it{\gamma}$ = 1', lw=2,
         linestyle='-', color='orange')
# plt.plot(star_2_mass_g1_new[:-5000]/star_1_mass_g1_new[:-5000], period_class_g1_new[:-5000],label='MESA default, $\gamma$ = 1',lw=2,linestyle='-',color='orange')
plt.plot(star_2_mass_g1_new[iRLOF_1_g1_new] / star_1_mass_g1_new[iRLOF_1_g1_new], period_class_g1_new[iRLOF_1_g1_new],
         lw=7, linestyle='-', color='orange')

plt.plot(star_2_mass_g2_new / star_1_mass_g2_new, period_class_g2_new, label='MESA default, $\it{\gamma}$ = 2', lw=2,
         linestyle='-', color='blue')
plt.plot(star_2_mass_g2_new[iRLOF_1_g2_new] / star_1_mass_g2_new[iRLOF_1_g2_new], period_class_g2_new[iRLOF_1_g2_new],
         lw=7, linestyle='-', color='blue')

plt.plot(star_2_mass_g10_new / star_1_mass_g10_new, period_class_g10_new, label='MESA default, $\it{\gamma}$ = 10',
         lw=2, linestyle='-', color='red')
plt.plot(star_2_mass_g10_new[iRLOF_1_g10_new] / star_1_mass_g10_new[iRLOF_1_g10_new],
         period_class_g10_new[iRLOF_1_g10_new], lw=7, linestyle='-', color='red')

plt.plot(star_2_mass_posydon / star_1_mass_posydon, period_posydon, linestyle='-', label='POSYDON, $\it{\gamma}$ = 0',
         color='green', lw=2)
plt.plot(star_2_mass_posydon[iRLOF_1_posydon] / star_1_mass_posydon[iRLOF_1_posydon], period_posydon[iRLOF_1_posydon],
         lw=7, c='green')

'''


if any(iRLOF_2_g1_new) == True:


   indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
   indx2 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_2_g1_new]))

   print('(len(iRLOF_1_g1_new) > 0):', indx1)
   print('(len(iRLOF_2_g1_new) > 0):', indx2)

   plt.plot(star_2_mass_g1_new[0:indx1]/star_1_mass_g1_new[0:indx1],period_class_g1_new[0:indx1],color='orange',linestyle='-',label='MESA default, $\gamma$ = 1',lw=2)
   plt.plot(star_2_mass_g1_new[indx1:indx2]/star_1_mass_g1_new[indx1:indx2],period_class_g1_new[indx1:indx2], lw=7, c='orange')
   plt.plot(star_2_mass_g1_new[indx2]/star_1_mass_g1_new[indx2],period_class_g1_new[indx2], marker='o', c='orange', mfc = 'orange', ms = 25)


else:
   indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
   print('(len(iRLOF_1_g1_new) > 0):', indx1)

   plt.plot(star_2_mass_g1_new[0:indx1]/star_1_mass_g1_new[0:indx1],period_class_g1_new[0:indx1],color='orange',linestyle='-',label='MESA default, $\gamma$ = 1',lw=2)
   plt.plot(star_2_mass_g1_new[indx1:-1]/star_1_mass_g1_new[indx1:-1],period_class_g1_new[indx1:-1], lw=7, c='orange')
   plt.plot(star_2_mass_g1_new[-2]/star_1_mass_g1_new[-2],period_class_g1_new[-2], marker='o', c='orange',ms=25 ,mfc = 'None')




if any(iRLOF_2_g2_new) == True:


   indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
   indx2 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_2_g2_new]))

   print('(len(iRLOF_1_g2_new) > 0):', indx1)
   print('(len(iRLOF_2_g2_new) > 0):', indx2)

   plt.plot(star_2_mass_g2_new[0:indx1]/star_1_mass_g2_new[0:indx1],period_class_g2_new[0:indx1],color='blue',linestyle='-',label='MESA default, $\gamma$ = 2',lw=2)
   plt.plot(star_2_mass_g2_new[indx1:indx2]/star_1_mass_g2_new[indx1:indx2],period_class_g2_new[indx1:indx2], lw=7, c='blue')
   plt.plot(star_2_mass_g2_new[indx2]/star_1_mass_g2_new[indx2],period_class_g2_new[indx2], marker='o', c='blue', mfc = 'blue', ms = 25)


else:
   indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
   print('(len(iRLOF_1_g2_new) > 0):', indx1)

   plt.plot(star_2_mass_g2_new[0:indx1]/star_1_mass_g2_new[0:indx1],period_class_g2_new[0:indx1],color='blue',linestyle='-',label='MESA default, $\gamma$ = 2',lw=2)
   plt.plot(star_2_mass_g2_new[indx1:-1]/star_1_mass_g2_new[indx1:-1],period_class_g2_new[indx1:-1], lw=7, c='blue')
   plt.plot(star_2_mass_g2_new[-1]/star_1_mass_g2_new[-1],period_class_g2_new[-1], marker='o', c='blue',ms=25 ,mfc = 'None')


if any(iRLOF_2_g10_new) == True:


   indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
   indx2 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_2_g10_new]))

   print('(len(iRLOF_1_g10_new) > 0):', indx1)
   print('(len(iRLOF_2_g10_new) > 0):', indx2)

   plt.plot(star_2_mass_g10_new[0:indx1]/star_1_mass_g10_new[0:indx1],period_class_g10_new[0:indx1],color='red',linestyle='-',label='MESA default, $\gamma$ = 10',lw=2)
   plt.plot(star_2_mass_g10_new[indx1:indx2]/star_1_mass_g10_new[indx1:indx2],period_class_g10_new[indx1:indx2], lw=7, c='red')
   plt.plot(star_2_mass_g10_new[indx2]/star_1_mass[indx2],period_class_g10_new[indx2], marker='o', c='red', mfc = 'red', ms = 25)


else:
   indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
   print('(len(iRLOF_1_g10_new) > 0):', indx1)

   plt.plot(star_2_mass_g10_new[0:indx1]/star_1_mass_g10_new[0:indx1],period_class_g10_new[0:indx1],color='red',linestyle='-',label='MESA default, $\gamma$ = 10',lw=2)
   plt.plot(star_2_mass_g10_new[indx1:-1]/star_1_mass_g10_new[indx1:-1],period_class_g10_new[indx1:-1], lw=7, c='red')
   plt.plot(star_2_mass_g10_new[-1]/star_1_mass_g10_new[-1],period_class_g10_new[-1], marker='o', c='red',ms=25 ,mfc = 'None')



if any(iRLOF_2_posydon) == True:


   indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
   indx2 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_2_posydon]))

   print('(len(iRLOF_1_posydon) > 0):', indx1)
   print('(len(iRLOF_2_posydon) > 0):', indx2)

   plt.plot(star_2_mass_posydon[0:indx1]/star_1_mass_posydon[0:indx1],period_posydon[0:indx1],color='green',linestyle='-',label='POSYDON, $\gamma$ = 0',lw=2)
   plt.plot(star_2_mass_posydon[indx1:indx2]/star_1_mass_posydon[indx1:indx2],period_posydon[indx1:indx2], lw=7, c='green')
   plt.plot(star_2_mass_posydon[indx2]/star_1_mass_posydon[indx2],period_posydon[indx2], marker='o', c='green', mfc = 'green', ms = 25)


else:
   indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
   print('(len(iRLOF_1_posydon) > 0):', indx1)

   plt.plot(star_2_mass_posydon[0:indx1]/star_1_mass_posydon[0:indx1],period_posydon[0:indx1],color='green',linestyle='-',label='POSYDON, $\gamma$ = 0',lw=2)
   plt.plot(star_2_mass_posydon[indx1:-1]/star_1_mass_posydon[indx1:-1],period_posydon[indx1:-1], lw=7, c='green')
   plt.plot(star_2_mass_posydon[-1]/star_1_mass_posydon[-1],period_posydon[-1], marker='o', c='green',ms=25 ,mfc = 'None')




'''

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.xlim([5e6,7.5e6])
plt.ylim([0.9, 1.4])  # 1.25 days

# plt.ylim([1.9,7.5]) # 5 days
# plt.ylim([1,5]) # 3 days
# plt.ylim([0,17]) # 10 days
# plt.ylim([0,125]) # 70 days

# plt.ylim([0,90]) # 50 days


plt.xlabel('$\it{Q}$', fontsize=25)
plt.ylabel('Period [days]', fontsize=25)
plt.legend(loc=2, fontsize=20)

##plt.xlim([0.3,1.8])
###plt.xlim([0.3,1.5])
# plt.xlim([0.64,0.72])
plt.xlim([0.64, 0.72])  # 1.25 days
# plt.xlim([0.25,1.8]) # 5 days
# plt.xlim([0.3,2]) # 3 days

# plt.xlim([0.3,1.8]) # 10 days
# plt.xlim([0.3,1.8]) # 70 days

plt.savefig(pp1_all_panel, format='pdf')

pp1_all_panel.close()

m_200 = mp.MESA()
m2_200 = mp.MESA()
m3_200 = mp.MESA()

m_200_newtides = mp.MESA()
m2_200_newtides = mp.MESA()
m3_200_newtides = mp.MESA()

m3_200_g1_new = mp.MESA()
m3_200_g2_new = mp.MESA()
m3_200_g10_new = mp.MESA()

name = 'post_interaction/30_20_3'

print(os.path.join(paths.data, name + '_g1_new/LOGS3'))
m3_200_g1_new.log_fold = os.path.join(paths.data, name + '_g1_new/LOGS3')
m3_200_g1_new.loadHistory()

m3_200_g2_new.log_fold = os.path.join(paths.data, name + '_g2_new/LOGS3')
# m3_200_g2_new.log_fold=name+'_g1_pos_new/LOGS3'

m3_200_g2_new.loadHistory()

m3_200_g10_new.log_fold = os.path.join(paths.data, name + '_g10_new/LOGS3')
m3_200_g10_new.loadHistory()

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

rl_relative_gap_1 = m3_200.hist.rl_relative_overflow_1
rl_relative_gap_2 = m3_200.hist.rl_relative_overflow_2
star_1_mass = m3_200.hist.star_1_mass
star_2_mass = m3_200.hist.star_2_mass
iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_2 = rl_relative_gap_2 > 0

period_class = m3_200.hist.period_days

rl_relative_gap_1_g1_new = m3_200_g1_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g1_new = m3_200_g1_new.hist.rl_relative_overflow_2

star_1_mass_g1_new = m3_200_g1_new.hist.star_1_mass
star_2_mass_g1_new = m3_200_g1_new.hist.star_2_mass
iRLOF_1_g1_new = rl_relative_gap_1_g1_new > 0
iRLOF_2_g1_new = rl_relative_gap_2_g1_new > 0
period_class_g1_new = m3_200_g1_new.hist.period_days

rl_relative_gap_1_g2_new = m3_200_g2_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g2_new = m3_200_g2_new.hist.rl_relative_overflow_2
star_1_mass_g2_new = m3_200_g2_new.hist.star_1_mass
star_2_mass_g2_new = m3_200_g2_new.hist.star_2_mass
iRLOF_1_g2_new = rl_relative_gap_1_g2_new > 0
iRLOF_2_g2_new = rl_relative_gap_2_g2_new > 0
period_class_g2_new = m3_200_g2_new.hist.period_days

rl_relative_gap_1_g10_new = m3_200_g10_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g10_new = m3_200_g10_new.hist.rl_relative_overflow_2
star_1_mass_g10_new = m3_200_g10_new.hist.star_1_mass
star_2_mass_g10_new = m3_200_g10_new.hist.star_2_mass
iRLOF_1_g10_new = rl_relative_gap_1_g10_new > 0
iRLOF_2_g10_new = rl_relative_gap_2_g10_new > 0
period_class_g10_new = m3_200_g10_new.hist.period_days

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

surf_avg_omega_1_class = m_200.hist.surf_avg_omega

surf_avg_omega_1_pos = m_200_newtides.hist.surf_avg_omega

surf_avg_omega_2_class = m2_200.hist.surf_avg_omega

surf_avg_omega_2_pos = m2_200_newtides.hist.surf_avg_omega

star_age_pos = m2_200_newtides.hist.star_age
star_age_class = m2_200.hist.star_age

# p1.log_fold='LOGS1'

# p1.loadProfile(num=-1)

# p=mp.plot()


rl_relative_gap_1_posydon = m3_200_newtides.hist.rl_relative_overflow_1
rl_relative_gap_2_posydon = m3_200_newtides.hist.rl_relative_overflow_2

age_class = m3_200.hist.age
age_posydon = m3_200_newtides.hist.age

lg_t_sync_2_class = m3_200.hist.lg_t_sync_2
lg_t_sync_2_posydon = m3_200_newtides.hist.lg_t_sync_2

lg_t_sync_1_class = m3_200.hist.lg_t_sync_1
lg_t_sync_1_posydon = m3_200_newtides.hist.lg_t_sync_1

iRLOF_1_posydon = rl_relative_gap_1_posydon > 0
iRLOF_2_posydon = rl_relative_gap_2_posydon > 0

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

pp1_all_panel3 = PdfPages(paths.figures / 'p_q_3days.pdf')

fig = plt.figure(figsize=(10, 10))

plt.title('$\it{P}_\mathrm{ini}$ = 3 [days]', fontsize=30)

# plt.plot(star_2_mass/star_1_mass,period_class,color='k',linestyle='-',label='MESA default, $\gamma$ = 0',lw=2)
# plt.plot(star_2_mass[iRLOF_1]/star_1_mass[iRLOF_1],period_class[iRLOF_1], lw=7, c='k')


if any(iRLOF_2) == True:

    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    indx2 = list(star_2_mass).index(min(star_2_mass[iRLOF_2]))

    print('(len(iRLOF_1) > 0):', indx1)
    print('(len(iRLOF_2) > 0):', indx2)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:indx2] / star_1_mass[indx1:indx2], period_class[indx1:indx2], lw=7, c='k')
    plt.plot(star_2_mass[indx2] / star_1_mass[indx2], period_class[indx2], marker='o', c='k', mfc='k', ms=25)


else:
    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    print('(len(iRLOF_1) > 0):', indx1)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:-1] / star_1_mass[indx1:-1], period_class[indx1:-1], lw=7, c='k')
    plt.plot(star_2_mass[-1] / star_1_mass[-1], period_class[-1], marker='o', c='k', ms=25, mfc='None')

if any(iRLOF_2_g1_new) == True:

    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    indx2 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_2_g1_new]))

    print('(len(iRLOF_1_g1_new) > 0):', indx1)
    print('(len(iRLOF_2_g1_new) > 0):', indx2)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:indx2] / star_1_mass_g1_new[indx1:indx2], period_class_g1_new[indx1:indx2], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[indx2] / star_1_mass_g1_new[indx2], period_class_g1_new[indx2], marker='o', c='orange',
             mfc='orange', ms=25)


else:
    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    print('(len(iRLOF_1_g1_new) > 0):', indx1)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:-1] / star_1_mass_g1_new[indx1:-1], period_class_g1_new[indx1:-1], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[-2] / star_1_mass_g1_new[-2], period_class_g1_new[-2], marker='o', c='orange', ms=25,
             mfc='None')

if any(iRLOF_2_g2_new) == True:

    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    indx2 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_2_g2_new]))

    print('(len(iRLOF_1_g2_new) > 0):', indx1)
    print('(len(iRLOF_2_g2_new) > 0):', indx2)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:indx2] / star_1_mass_g2_new[indx1:indx2], period_class_g2_new[indx1:indx2], lw=7,
             c='blue')
    plt.plot(star_2_mass_g2_new[indx2] / star_1_mass_g2_new[indx2], period_class_g2_new[indx2], marker='o', c='blue',
             mfc='blue', ms=25)


else:
    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    print('(len(iRLOF_1_g2_new) > 0):', indx1)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:-1] / star_1_mass_g2_new[indx1:-1], period_class_g2_new[indx1:-1], lw=7, c='blue')
    plt.plot(star_2_mass_g2_new[-1] / star_1_mass_g2_new[-1], period_class_g2_new[-1], marker='o', c='blue', ms=25,
             mfc='None')

if any(iRLOF_2_g10_new) == True:

    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    indx2 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_2_g10_new]))

    print('(len(iRLOF_1_g10_new) > 0):', indx1)
    print('(len(iRLOF_2_g10_new) > 0):', indx2)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:indx2] / star_1_mass_g10_new[indx1:indx2], period_class_g10_new[indx1:indx2],
             lw=7, c='red')
    plt.plot(star_2_mass_g10_new[indx2] / star_1_mass[indx2], period_class_g10_new[indx2], marker='o', c='red',
             mfc='red', ms=25)


else:
    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    print('(len(iRLOF_1_g10_new) > 0):', indx1)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:-1] / star_1_mass_g10_new[indx1:-1], period_class_g10_new[indx1:-1], lw=7,
             c='red')
    plt.plot(star_2_mass_g10_new[-1] / star_1_mass_g10_new[-1], period_class_g10_new[-1], marker='o', c='red', ms=25,
             mfc='None')

if any(iRLOF_2_posydon) == True:

    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    indx2 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_2_posydon]))

    print('(len(iRLOF_1_posydon) > 0):', indx1)
    print('(len(iRLOF_2_posydon) > 0):', indx2)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:indx2] / star_1_mass_posydon[indx1:indx2], period_posydon[indx1:indx2], lw=7,
             c='green')
    plt.plot(star_2_mass_posydon[indx2] / star_1_mass_posydon[indx2], period_posydon[indx2], marker='o', c='green',
             mfc='green', ms=25)


else:
    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    print('(len(iRLOF_1_posydon) > 0):', indx1)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:-1] / star_1_mass_posydon[indx1:-1], period_posydon[indx1:-1], lw=7, c='green')
    plt.plot(star_2_mass_posydon[-1] / star_1_mass_posydon[-1], period_posydon[-1], marker='o', c='green', ms=25,
             mfc='None')

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.xlim([5e6,7.5e6])
# plt.ylim([0.9,1.4]) # 1.25 days

# plt.ylim([1.9,7.5]) # 5 days
plt.ylim([1, 5])  # 3 days
# plt.ylim([0,17]) # 10 days
# plt.ylim([0,125]) # 70 days

# plt.ylim([0,90]) # 50 days


plt.xlabel('$\it{Q}$', fontsize=25)
plt.ylabel('Period [days]', fontsize=25)
plt.legend(loc=2, fontsize=20)

##plt.xlim([0.3,1.8])
###plt.xlim([0.3,1.5])
# plt.xlim([0.64,0.72])
# plt.xlim([0.64,0.72]) # 1.25 days
# plt.xlim([0.25,1.8]) # 5 days
plt.xlim([0.3, 2])  # 3 days

# plt.xlim([0.3,1.8]) # 10 days
# plt.xlim([0.3,1.8]) # 70 days

plt.savefig(pp1_all_panel3, format='pdf')

pp1_all_panel3.close()

m_200 = mp.MESA()
m2_200 = mp.MESA()
m3_200 = mp.MESA()

m_200_newtides = mp.MESA()
m2_200_newtides = mp.MESA()
m3_200_newtides = mp.MESA()

m3_200_g1_new = mp.MESA()
m3_200_g2_new = mp.MESA()
m3_200_g10_new = mp.MESA()

m_200 = mp.MESA()
m2_200 = mp.MESA()
m3_200 = mp.MESA()

m_200_newtides = mp.MESA()
m2_200_newtides = mp.MESA()
m3_200_newtides = mp.MESA()

m3_200_g1_new = mp.MESA()
m3_200_g2_new = mp.MESA()
m3_200_g10_new = mp.MESA()

name = 'post_interaction/30_20_5'

print(os.path.join(paths.data, name + '_g1_new/LOGS3'))
m3_200_g1_new.log_fold = os.path.join(paths.data, name + '_g1_new/LOGS3')
m3_200_g1_new.loadHistory()

m3_200_g2_new.log_fold = os.path.join(paths.data, name + '_g2_new/LOGS3')
# m3_200_g2_new.log_fold=name+'_g1_pos_new/LOGS3'

m3_200_g2_new.loadHistory()

m3_200_g10_new.log_fold = os.path.join(paths.data, name + '_g10_new/LOGS3')
m3_200_g10_new.loadHistory()

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

rl_relative_gap_1 = m3_200.hist.rl_relative_overflow_1
rl_relative_gap_2 = m3_200.hist.rl_relative_overflow_2
star_1_mass = m3_200.hist.star_1_mass
star_2_mass = m3_200.hist.star_2_mass
iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_2 = rl_relative_gap_2 > 0

period_class = m3_200.hist.period_days

rl_relative_gap_1_g1_new = m3_200_g1_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g1_new = m3_200_g1_new.hist.rl_relative_overflow_2

star_1_mass_g1_new = m3_200_g1_new.hist.star_1_mass
star_2_mass_g1_new = m3_200_g1_new.hist.star_2_mass
iRLOF_1_g1_new = rl_relative_gap_1_g1_new > 0
iRLOF_2_g1_new = rl_relative_gap_2_g1_new > 0
period_class_g1_new = m3_200_g1_new.hist.period_days

rl_relative_gap_1_g2_new = m3_200_g2_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g2_new = m3_200_g2_new.hist.rl_relative_overflow_2
star_1_mass_g2_new = m3_200_g2_new.hist.star_1_mass
star_2_mass_g2_new = m3_200_g2_new.hist.star_2_mass
iRLOF_1_g2_new = rl_relative_gap_1_g2_new > 0
iRLOF_2_g2_new = rl_relative_gap_2_g2_new > 0
period_class_g2_new = m3_200_g2_new.hist.period_days

rl_relative_gap_1_g10_new = m3_200_g10_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g10_new = m3_200_g10_new.hist.rl_relative_overflow_2
star_1_mass_g10_new = m3_200_g10_new.hist.star_1_mass
star_2_mass_g10_new = m3_200_g10_new.hist.star_2_mass
iRLOF_1_g10_new = rl_relative_gap_1_g10_new > 0
iRLOF_2_g10_new = rl_relative_gap_2_g10_new > 0
period_class_g10_new = m3_200_g10_new.hist.period_days

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

surf_avg_omega_1_class = m_200.hist.surf_avg_omega

surf_avg_omega_1_pos = m_200_newtides.hist.surf_avg_omega

surf_avg_omega_2_class = m2_200.hist.surf_avg_omega

surf_avg_omega_2_pos = m2_200_newtides.hist.surf_avg_omega

star_age_pos = m2_200_newtides.hist.star_age
star_age_class = m2_200.hist.star_age

# p1.log_fold='LOGS1'

# p1.loadProfile(num=-1)

# p=mp.plot()


rl_relative_gap_1_posydon = m3_200_newtides.hist.rl_relative_overflow_1
rl_relative_gap_2_posydon = m3_200_newtides.hist.rl_relative_overflow_2

age_class = m3_200.hist.age
age_posydon = m3_200_newtides.hist.age

lg_t_sync_2_class = m3_200.hist.lg_t_sync_2
lg_t_sync_2_posydon = m3_200_newtides.hist.lg_t_sync_2

lg_t_sync_1_class = m3_200.hist.lg_t_sync_1
lg_t_sync_1_posydon = m3_200_newtides.hist.lg_t_sync_1

iRLOF_1_posydon = rl_relative_gap_1_posydon > 0
iRLOF_2_posydon = rl_relative_gap_2_posydon > 0

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

pp1_all_panel5 = PdfPages(paths.figures / 'p_q_5days.pdf')

fig = plt.figure(figsize=(10, 10))

plt.title('$\it{P}_\mathrm{ini}$ = 5 [days]', fontsize=30)

# plt.plot(star_2_mass/star_1_mass,period_class,color='k',linestyle='-',label='MESA default, $\gamma$ = 0',lw=2)
# plt.plot(star_2_mass[iRLOF_1]/star_1_mass[iRLOF_1],period_class[iRLOF_1], lw=7, c='k')


if any(iRLOF_2) == True:

    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    indx2 = list(star_2_mass).index(min(star_2_mass[iRLOF_2]))

    print('(len(iRLOF_1) > 0):', indx1)
    print('(len(iRLOF_2) > 0):', indx2)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:indx2] / star_1_mass[indx1:indx2], period_class[indx1:indx2], lw=7, c='k')
    plt.plot(star_2_mass[indx2] / star_1_mass[indx2], period_class[indx2], marker='o', c='k', mfc='k', ms=25)


else:
    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    print('(len(iRLOF_1) > 0):', indx1)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:-1] / star_1_mass[indx1:-1], period_class[indx1:-1], lw=7, c='k')
    plt.plot(star_2_mass[-1] / star_1_mass[-1], period_class[-1], marker='o', c='k', ms=25, mfc='None')

if any(iRLOF_2_g1_new) == True:

    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    indx2 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_2_g1_new]))

    print('(len(iRLOF_1_g1_new) > 0):', indx1)
    print('(len(iRLOF_2_g1_new) > 0):', indx2)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:indx2] / star_1_mass_g1_new[indx1:indx2], period_class_g1_new[indx1:indx2], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[indx2] / star_1_mass_g1_new[indx2], period_class_g1_new[indx2], marker='o', c='orange',
             mfc='orange', ms=25)


else:
    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    print('(len(iRLOF_1_g1_new) > 0):', indx1)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:-1] / star_1_mass_g1_new[indx1:-1], period_class_g1_new[indx1:-1], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[-2] / star_1_mass_g1_new[-2], period_class_g1_new[-2], marker='o', c='orange', ms=25,
             mfc='None')

if any(iRLOF_2_g2_new) == True:

    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    indx2 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_2_g2_new]))

    print('(len(iRLOF_1_g2_new) > 0):', indx1)
    print('(len(iRLOF_2_g2_new) > 0):', indx2)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:indx2] / star_1_mass_g2_new[indx1:indx2], period_class_g2_new[indx1:indx2], lw=7,
             c='blue')
    plt.plot(star_2_mass_g2_new[indx2] / star_1_mass_g2_new[indx2], period_class_g2_new[indx2], marker='o', c='blue',
             mfc='blue', ms=25)


else:
    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    print('(len(iRLOF_1_g2_new) > 0):', indx1)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:-1] / star_1_mass_g2_new[indx1:-1], period_class_g2_new[indx1:-1], lw=7, c='blue')
    plt.plot(star_2_mass_g2_new[-1] / star_1_mass_g2_new[-1], period_class_g2_new[-1], marker='o', c='blue', ms=25,
             mfc='None')

if any(iRLOF_2_g10_new) == True:

    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    indx2 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_2_g10_new]))

    print('(len(iRLOF_1_g10_new) > 0):', indx1)
    print('(len(iRLOF_2_g10_new) > 0):', indx2)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:indx2] / star_1_mass_g10_new[indx1:indx2], period_class_g10_new[indx1:indx2],
             lw=7, c='red')
    plt.plot(star_2_mass_g10_new[indx2] / star_1_mass[indx2], period_class_g10_new[indx2], marker='o', c='red',
             mfc='red', ms=25)


else:
    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    print('(len(iRLOF_1_g10_new) > 0):', indx1)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:-1] / star_1_mass_g10_new[indx1:-1], period_class_g10_new[indx1:-1], lw=7,
             c='red')
    plt.plot(star_2_mass_g10_new[-1] / star_1_mass_g10_new[-1], period_class_g10_new[-1], marker='o', c='red', ms=25,
             mfc='None')

if any(iRLOF_2_posydon) == True:

    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    indx2 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_2_posydon]))

    print('(len(iRLOF_1_posydon) > 0):', indx1)
    print('(len(iRLOF_2_posydon) > 0):', indx2)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:indx2] / star_1_mass_posydon[indx1:indx2], period_posydon[indx1:indx2], lw=7,
             c='green')
    plt.plot(star_2_mass_posydon[indx2] / star_1_mass_posydon[indx2], period_posydon[indx2], marker='o', c='green',
             mfc='green', ms=25)


else:
    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    print('(len(iRLOF_1_posydon) > 0):', indx1)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:-1] / star_1_mass_posydon[indx1:-1], period_posydon[indx1:-1], lw=7, c='green')
    plt.plot(star_2_mass_posydon[-1] / star_1_mass_posydon[-1], period_posydon[-1], marker='o', c='green', ms=25,
             mfc='None')

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.xlim([5e6,7.5e6])
# plt.ylim([0.9,1.4]) # 1.25 days

plt.ylim([1.9, 7.5])  # 5 days
# plt.ylim([1,5]) # 3 days
# plt.ylim([0,17]) # 10 days
# plt.ylim([0,125]) # 70 days

# plt.ylim([0,90]) # 50 days


plt.xlabel('$\it{Q}$', fontsize=25)
plt.ylabel('Period [days]', fontsize=25)
plt.legend(loc=2, fontsize=20)

##plt.xlim([0.3,1.8])
###plt.xlim([0.3,1.5])
# plt.xlim([0.64,0.72])
# plt.xlim([0.64,0.72]) # 1.25 days
plt.xlim([0.25, 1.8])  # 5 days
# plt.xlim([0.3,2]) # 3 days

# plt.xlim([0.3,1.8]) # 10 days
# plt.xlim([0.3,1.8]) # 70 days

plt.savefig(pp1_all_panel5, format='pdf')

pp1_all_panel5.close()

m_200 = mp.MESA()
m2_200 = mp.MESA()
m3_200 = mp.MESA()

m_200_newtides = mp.MESA()
m2_200_newtides = mp.MESA()
m3_200_newtides = mp.MESA()

m3_200_g1_new = mp.MESA()
m3_200_g2_new = mp.MESA()
m3_200_g10_new = mp.MESA()

name = 'post_interaction/30_20_10'

print(os.path.join(paths.data, name + '_g1_new/LOGS3'))
m3_200_g1_new.log_fold = os.path.join(paths.data, name + '_g1_new/LOGS3')
m3_200_g1_new.loadHistory()

m3_200_g2_new.log_fold = os.path.join(paths.data, name + '_g2_new/LOGS3')
# m3_200_g2_new.log_fold=name+'_g1_pos_new/LOGS3'

m3_200_g2_new.loadHistory()

m3_200_g10_new.log_fold = os.path.join(paths.data, name + '_g10_new/LOGS3')
m3_200_g10_new.loadHistory()

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

rl_relative_gap_1 = m3_200.hist.rl_relative_overflow_1
rl_relative_gap_2 = m3_200.hist.rl_relative_overflow_2
star_1_mass = m3_200.hist.star_1_mass
star_2_mass = m3_200.hist.star_2_mass
iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_2 = rl_relative_gap_2 > 0

period_class = m3_200.hist.period_days

rl_relative_gap_1_g1_new = m3_200_g1_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g1_new = m3_200_g1_new.hist.rl_relative_overflow_2

star_1_mass_g1_new = m3_200_g1_new.hist.star_1_mass
star_2_mass_g1_new = m3_200_g1_new.hist.star_2_mass
iRLOF_1_g1_new = rl_relative_gap_1_g1_new > 0
iRLOF_2_g1_new = rl_relative_gap_2_g1_new > 0
period_class_g1_new = m3_200_g1_new.hist.period_days

rl_relative_gap_1_g2_new = m3_200_g2_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g2_new = m3_200_g2_new.hist.rl_relative_overflow_2
star_1_mass_g2_new = m3_200_g2_new.hist.star_1_mass
star_2_mass_g2_new = m3_200_g2_new.hist.star_2_mass
iRLOF_1_g2_new = rl_relative_gap_1_g2_new > 0
iRLOF_2_g2_new = rl_relative_gap_2_g2_new > 0
period_class_g2_new = m3_200_g2_new.hist.period_days

rl_relative_gap_1_g10_new = m3_200_g10_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g10_new = m3_200_g10_new.hist.rl_relative_overflow_2
star_1_mass_g10_new = m3_200_g10_new.hist.star_1_mass
star_2_mass_g10_new = m3_200_g10_new.hist.star_2_mass
iRLOF_1_g10_new = rl_relative_gap_1_g10_new > 0
iRLOF_2_g10_new = rl_relative_gap_2_g10_new > 0
period_class_g10_new = m3_200_g10_new.hist.period_days

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

surf_avg_omega_1_class = m_200.hist.surf_avg_omega

surf_avg_omega_1_pos = m_200_newtides.hist.surf_avg_omega

surf_avg_omega_2_class = m2_200.hist.surf_avg_omega

surf_avg_omega_2_pos = m2_200_newtides.hist.surf_avg_omega

star_age_pos = m2_200_newtides.hist.star_age
star_age_class = m2_200.hist.star_age

# p1.log_fold='LOGS1'

# p1.loadProfile(num=-1)

# p=mp.plot()


rl_relative_gap_1_posydon = m3_200_newtides.hist.rl_relative_overflow_1
rl_relative_gap_2_posydon = m3_200_newtides.hist.rl_relative_overflow_2

age_class = m3_200.hist.age
age_posydon = m3_200_newtides.hist.age

lg_t_sync_2_class = m3_200.hist.lg_t_sync_2
lg_t_sync_2_posydon = m3_200_newtides.hist.lg_t_sync_2

lg_t_sync_1_class = m3_200.hist.lg_t_sync_1
lg_t_sync_1_posydon = m3_200_newtides.hist.lg_t_sync_1

iRLOF_1_posydon = rl_relative_gap_1_posydon > 0
iRLOF_2_posydon = rl_relative_gap_2_posydon > 0

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

pp1_all_panel10 = PdfPages(paths.figures / 'p_q_10days.pdf')

fig = plt.figure(figsize=(10, 10))

plt.title('$\it{P}_\mathrm{ini}$ = 10 [days]', fontsize=30)

# plt.plot(star_2_mass/star_1_mass,period_class,color='k',linestyle='-',label='MESA default, $\gamma$ = 0',lw=2)
# plt.plot(star_2_mass[iRLOF_1]/star_1_mass[iRLOF_1],period_class[iRLOF_1], lw=7, c='k')


if any(iRLOF_2) == True:

    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    indx2 = list(star_2_mass).index(min(star_2_mass[iRLOF_2]))

    print('(len(iRLOF_1) > 0):', indx1)
    print('(len(iRLOF_2) > 0):', indx2)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:indx2] / star_1_mass[indx1:indx2], period_class[indx1:indx2], lw=7, c='k')
    plt.plot(star_2_mass[indx2] / star_1_mass[indx2], period_class[indx2], marker='o', c='k', mfc='k', ms=25)


else:
    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    print('(len(iRLOF_1) > 0):', indx1)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:-1] / star_1_mass[indx1:-1], period_class[indx1:-1], lw=7, c='k')
    plt.plot(star_2_mass[-1] / star_1_mass[-1], period_class[-1], marker='o', c='k', ms=25, mfc='None')

if any(iRLOF_2_g1_new) == True:

    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    indx2 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_2_g1_new]))

    print('(len(iRLOF_1_g1_new) > 0):', indx1)
    print('(len(iRLOF_2_g1_new) > 0):', indx2)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:indx2] / star_1_mass_g1_new[indx1:indx2], period_class_g1_new[indx1:indx2], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[indx2] / star_1_mass_g1_new[indx2], period_class_g1_new[indx2], marker='o', c='orange',
             mfc='orange', ms=25)


else:
    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    print('(len(iRLOF_1_g1_new) > 0):', indx1)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:-1] / star_1_mass_g1_new[indx1:-1], period_class_g1_new[indx1:-1], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[-2] / star_1_mass_g1_new[-2], period_class_g1_new[-2], marker='o', c='orange', ms=25,
             mfc='None')

if any(iRLOF_2_g2_new) == True:

    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    indx2 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_2_g2_new]))

    print('(len(iRLOF_1_g2_new) > 0):', indx1)
    print('(len(iRLOF_2_g2_new) > 0):', indx2)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:indx2] / star_1_mass_g2_new[indx1:indx2], period_class_g2_new[indx1:indx2], lw=7,
             c='blue')
    plt.plot(star_2_mass_g2_new[indx2] / star_1_mass_g2_new[indx2], period_class_g2_new[indx2], marker='o', c='blue',
             mfc='blue', ms=25)


else:
    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    print('(len(iRLOF_1_g2_new) > 0):', indx1)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:-1] / star_1_mass_g2_new[indx1:-1], period_class_g2_new[indx1:-1], lw=7, c='blue')
    plt.plot(star_2_mass_g2_new[-1] / star_1_mass_g2_new[-1], period_class_g2_new[-1], marker='o', c='blue', ms=25,
             mfc='None')

if any(iRLOF_2_g10_new) == True:

    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    indx2 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_2_g10_new]))

    print('(len(iRLOF_1_g10_new) > 0):', indx1)
    print('(len(iRLOF_2_g10_new) > 0):', indx2)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:indx2] / star_1_mass_g10_new[indx1:indx2], period_class_g10_new[indx1:indx2],
             lw=7, c='red')
    plt.plot(star_2_mass_g10_new[indx2] / star_1_mass[indx2], period_class_g10_new[indx2], marker='o', c='red',
             mfc='red', ms=25)


else:
    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    print('(len(iRLOF_1_g10_new) > 0):', indx1)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:-1] / star_1_mass_g10_new[indx1:-1], period_class_g10_new[indx1:-1], lw=7,
             c='red')
    plt.plot(star_2_mass_g10_new[-1] / star_1_mass_g10_new[-1], period_class_g10_new[-1], marker='o', c='red', ms=25,
             mfc='None')

if any(iRLOF_2_posydon) == True:

    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    indx2 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_2_posydon]))

    print('(len(iRLOF_1_posydon) > 0):', indx1)
    print('(len(iRLOF_2_posydon) > 0):', indx2)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:indx2] / star_1_mass_posydon[indx1:indx2], period_posydon[indx1:indx2], lw=7,
             c='green')
    plt.plot(star_2_mass_posydon[indx2] / star_1_mass_posydon[indx2], period_posydon[indx2], marker='o', c='green',
             mfc='green', ms=25)


else:
    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    print('(len(iRLOF_1_posydon) > 0):', indx1)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:-1] / star_1_mass_posydon[indx1:-1], period_posydon[indx1:-1], lw=7, c='green')
    plt.plot(star_2_mass_posydon[-1] / star_1_mass_posydon[-1], period_posydon[-1], marker='o', c='green', ms=25,
             mfc='None')

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.xlim([5e6,7.5e6])
# plt.ylim([0.9,1.4]) # 1.25 days

# plt.ylim([1.9,7.5]) # 5 days
# plt.ylim([1,5]) # 3 days
plt.ylim([0, 17])  # 10 days
# plt.ylim([0,125]) # 70 days

# plt.ylim([0,90]) # 50 days


plt.xlabel('$\it{Q}$', fontsize=25)
plt.ylabel('Period [days]', fontsize=25)
plt.legend(loc=2, fontsize=20)

##plt.xlim([0.3,1.8])
###plt.xlim([0.3,1.5])
# plt.xlim([0.64,0.72])
# plt.xlim([0.64,0.72]) # 1.25 days
# plt.xlim([0.25,1.8]) # 5 days
# plt.xlim([0.3,2]) # 3 days

plt.xlim([0.3, 1.8])  # 10 days
# plt.xlim([0.3,1.8]) # 70 days

plt.savefig(pp1_all_panel10, format='pdf')

pp1_all_panel10.close()

m_200 = mp.MESA()
m2_200 = mp.MESA()
m3_200 = mp.MESA()

m_200_newtides = mp.MESA()
m2_200_newtides = mp.MESA()
m3_200_newtides = mp.MESA()

m3_200_g1_new = mp.MESA()
m3_200_g2_new = mp.MESA()
m3_200_g10_new = mp.MESA()

name = 'post_interaction/30_20_50'

print(os.path.join(paths.data, name + '_g1_new/LOGS3'))
m3_200_g1_new.log_fold = os.path.join(paths.data, name + '_g1_new/LOGS3')
m3_200_g1_new.loadHistory()

m3_200_g2_new.log_fold = os.path.join(paths.data, name + '_g2_new/LOGS3')
# m3_200_g2_new.log_fold=name+'_g1_pos_new/LOGS3'

m3_200_g2_new.loadHistory()

m3_200_g10_new.log_fold = os.path.join(paths.data, name + '_g10_new/LOGS3')
m3_200_g10_new.loadHistory()

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

rl_relative_gap_1 = m3_200.hist.rl_relative_overflow_1
rl_relative_gap_2 = m3_200.hist.rl_relative_overflow_2
star_1_mass = m3_200.hist.star_1_mass
star_2_mass = m3_200.hist.star_2_mass
iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_2 = rl_relative_gap_2 > 0

period_class = m3_200.hist.period_days

rl_relative_gap_1_g1_new = m3_200_g1_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g1_new = m3_200_g1_new.hist.rl_relative_overflow_2

star_1_mass_g1_new = m3_200_g1_new.hist.star_1_mass
star_2_mass_g1_new = m3_200_g1_new.hist.star_2_mass
iRLOF_1_g1_new = rl_relative_gap_1_g1_new > 0
iRLOF_2_g1_new = rl_relative_gap_2_g1_new > 0
period_class_g1_new = m3_200_g1_new.hist.period_days

rl_relative_gap_1_g2_new = m3_200_g2_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g2_new = m3_200_g2_new.hist.rl_relative_overflow_2
star_1_mass_g2_new = m3_200_g2_new.hist.star_1_mass
star_2_mass_g2_new = m3_200_g2_new.hist.star_2_mass
iRLOF_1_g2_new = rl_relative_gap_1_g2_new > 0
iRLOF_2_g2_new = rl_relative_gap_2_g2_new > 0
period_class_g2_new = m3_200_g2_new.hist.period_days

rl_relative_gap_1_g10_new = m3_200_g10_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g10_new = m3_200_g10_new.hist.rl_relative_overflow_2
star_1_mass_g10_new = m3_200_g10_new.hist.star_1_mass
star_2_mass_g10_new = m3_200_g10_new.hist.star_2_mass
iRLOF_1_g10_new = rl_relative_gap_1_g10_new > 0
iRLOF_2_g10_new = rl_relative_gap_2_g10_new > 0
period_class_g10_new = m3_200_g10_new.hist.period_days

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

surf_avg_omega_1_class = m_200.hist.surf_avg_omega

surf_avg_omega_1_pos = m_200_newtides.hist.surf_avg_omega

surf_avg_omega_2_class = m2_200.hist.surf_avg_omega

surf_avg_omega_2_pos = m2_200_newtides.hist.surf_avg_omega

star_age_pos = m2_200_newtides.hist.star_age
star_age_class = m2_200.hist.star_age

# p1.log_fold='LOGS1'

# p1.loadProfile(num=-1)

# p=mp.plot()


rl_relative_gap_1_posydon = m3_200_newtides.hist.rl_relative_overflow_1
rl_relative_gap_2_posydon = m3_200_newtides.hist.rl_relative_overflow_2

age_class = m3_200.hist.age
age_posydon = m3_200_newtides.hist.age

lg_t_sync_2_class = m3_200.hist.lg_t_sync_2
lg_t_sync_2_posydon = m3_200_newtides.hist.lg_t_sync_2

lg_t_sync_1_class = m3_200.hist.lg_t_sync_1
lg_t_sync_1_posydon = m3_200_newtides.hist.lg_t_sync_1

iRLOF_1_posydon = rl_relative_gap_1_posydon > 0
iRLOF_2_posydon = rl_relative_gap_2_posydon > 0

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

pp1_all_panel50 = PdfPages(paths.figures / 'p_q_50days.pdf')

fig = plt.figure(figsize=(10, 10))

plt.title('$\it{P}_\mathrm{ini}$ = 50 [days]', fontsize=30)

# plt.plot(star_2_mass/star_1_mass,period_class,color='k',linestyle='-',label='MESA default, $\gamma$ = 0',lw=2)
# plt.plot(star_2_mass[iRLOF_1]/star_1_mass[iRLOF_1],period_class[iRLOF_1], lw=7, c='k')


if any(iRLOF_2) == True:

    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    indx2 = list(star_2_mass).index(min(star_2_mass[iRLOF_2]))

    print('(len(iRLOF_1) > 0):', indx1)
    print('(len(iRLOF_2) > 0):', indx2)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:indx2] / star_1_mass[indx1:indx2], period_class[indx1:indx2], lw=7, c='k')
    plt.plot(star_2_mass[indx2] / star_1_mass[indx2], period_class[indx2], marker='o', c='k', mfc='k', ms=25)


else:
    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    print('(len(iRLOF_1) > 0):', indx1)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:-1] / star_1_mass[indx1:-1], period_class[indx1:-1], lw=7, c='k')
    plt.plot(star_2_mass[-1] / star_1_mass[-1], period_class[-1], marker='o', c='k', ms=25, mfc='None')

if any(iRLOF_2_g1_new) == True:

    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    indx2 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_2_g1_new]))

    print('(len(iRLOF_1_g1_new) > 0):', indx1)
    print('(len(iRLOF_2_g1_new) > 0):', indx2)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:indx2] / star_1_mass_g1_new[indx1:indx2], period_class_g1_new[indx1:indx2], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[indx2] / star_1_mass_g1_new[indx2], period_class_g1_new[indx2], marker='o', c='orange',
             mfc='orange', ms=25)


else:
    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    print('(len(iRLOF_1_g1_new) > 0):', indx1)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:-1] / star_1_mass_g1_new[indx1:-1], period_class_g1_new[indx1:-1], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[-2] / star_1_mass_g1_new[-2], period_class_g1_new[-2], marker='o', c='orange', ms=25,
             mfc='None')

if any(iRLOF_2_g2_new) == True:

    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    indx2 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_2_g2_new]))

    print('(len(iRLOF_1_g2_new) > 0):', indx1)
    print('(len(iRLOF_2_g2_new) > 0):', indx2)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:indx2] / star_1_mass_g2_new[indx1:indx2], period_class_g2_new[indx1:indx2], lw=7,
             c='blue')
    plt.plot(star_2_mass_g2_new[indx2] / star_1_mass_g2_new[indx2], period_class_g2_new[indx2], marker='o', c='blue',
             mfc='blue', ms=25)


else:
    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    print('(len(iRLOF_1_g2_new) > 0):', indx1)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$= 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:-1] / star_1_mass_g2_new[indx1:-1], period_class_g2_new[indx1:-1], lw=7, c='blue')
    plt.plot(star_2_mass_g2_new[-1] / star_1_mass_g2_new[-1], period_class_g2_new[-1], marker='o', c='blue', ms=25,
             mfc='None')

if any(iRLOF_2_g10_new) == True:

    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    indx2 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_2_g10_new]))

    print('(len(iRLOF_1_g10_new) > 0):', indx1)
    print('(len(iRLOF_2_g10_new) > 0):', indx2)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:indx2] / star_1_mass_g10_new[indx1:indx2], period_class_g10_new[indx1:indx2],
             lw=7, c='red')
    plt.plot(star_2_mass_g10_new[indx2] / star_1_mass[indx2], period_class_g10_new[indx2], marker='o', c='red',
             mfc='red', ms=25)


else:
    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    print('(len(iRLOF_1_g10_new) > 0):', indx1)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:-1] / star_1_mass_g10_new[indx1:-1], period_class_g10_new[indx1:-1], lw=7,
             c='red')
    plt.plot(star_2_mass_g10_new[-1] / star_1_mass_g10_new[-1], period_class_g10_new[-1], marker='o', c='red', ms=25,
             mfc='None')

if any(iRLOF_2_posydon) == True:

    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    indx2 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_2_posydon]))

    print('(len(iRLOF_1_posydon) > 0):', indx1)
    print('(len(iRLOF_2_posydon) > 0):', indx2)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:indx2] / star_1_mass_posydon[indx1:indx2], period_posydon[indx1:indx2], lw=7,
             c='green')
    plt.plot(star_2_mass_posydon[indx2] / star_1_mass_posydon[indx2], period_posydon[indx2], marker='o', c='green',
             mfc='green', ms=25)


else:
    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    print('(len(iRLOF_1_posydon) > 0):', indx1)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:-1] / star_1_mass_posydon[indx1:-1], period_posydon[indx1:-1], lw=7, c='green')
    plt.plot(star_2_mass_posydon[-1] / star_1_mass_posydon[-1], period_posydon[-1], marker='o', c='green', ms=25,
             mfc='None')

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.xlim([5e6,7.5e6])
# plt.ylim([0.9,1.4]) # 1.25 days

# plt.ylim([1.9,7.5]) # 5 days
# plt.ylim([1,5]) # 3 days
# plt.ylim([0,17]) # 10 days
# plt.ylim([0,125]) # 70 days

plt.ylim([0, 90])  # 50 days

plt.xlabel('$\it{Q}$', fontsize=25)
plt.ylabel('Period [days]', fontsize=25)
plt.legend(loc=2, fontsize=20)

##plt.xlim([0.3,1.8])
###plt.xlim([0.3,1.5])
# plt.xlim([0.64,0.72])
# plt.xlim([0.64,0.72]) # 1.25 days
# plt.xlim([0.25,1.8]) # 5 days
# plt.xlim([0.3,2]) # 3 days

plt.xlim([0.3, 1.8])  # 10 days
# plt.xlim([0.3,1.8]) # 70 days

plt.savefig(pp1_all_panel50, format='pdf')

pp1_all_panel50.close()

m_200 = mp.MESA()
m2_200 = mp.MESA()
m3_200 = mp.MESA()

m_200_newtides = mp.MESA()
m2_200_newtides = mp.MESA()
m3_200_newtides = mp.MESA()

m3_200_g1_new = mp.MESA()
m3_200_g2_new = mp.MESA()
m3_200_g10_new = mp.MESA()

name = 'post_interaction/30_20_70'

print(os.path.join(paths.data, name + '_g1_new/LOGS3'))
m3_200_g1_new.log_fold = os.path.join(paths.data, name + '_g1_new/LOGS3')
m3_200_g1_new.loadHistory()

m3_200_g2_new.log_fold = os.path.join(paths.data, name + '_g2_new/LOGS3')
# m3_200_g2_new.log_fold=name+'_g1_pos_new/LOGS3'

m3_200_g2_new.loadHistory()

m3_200_g10_new.log_fold = os.path.join(paths.data, name + '_g10_new/LOGS3')
m3_200_g10_new.loadHistory()

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

rl_relative_gap_1 = m3_200.hist.rl_relative_overflow_1
rl_relative_gap_2 = m3_200.hist.rl_relative_overflow_2
star_1_mass = m3_200.hist.star_1_mass
star_2_mass = m3_200.hist.star_2_mass
iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_2 = rl_relative_gap_2 > 0

period_class = m3_200.hist.period_days

rl_relative_gap_1_g1_new = m3_200_g1_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g1_new = m3_200_g1_new.hist.rl_relative_overflow_2

star_1_mass_g1_new = m3_200_g1_new.hist.star_1_mass
star_2_mass_g1_new = m3_200_g1_new.hist.star_2_mass
iRLOF_1_g1_new = rl_relative_gap_1_g1_new > 0
iRLOF_2_g1_new = rl_relative_gap_2_g1_new > 0
period_class_g1_new = m3_200_g1_new.hist.period_days

rl_relative_gap_1_g2_new = m3_200_g2_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g2_new = m3_200_g2_new.hist.rl_relative_overflow_2
star_1_mass_g2_new = m3_200_g2_new.hist.star_1_mass
star_2_mass_g2_new = m3_200_g2_new.hist.star_2_mass
iRLOF_1_g2_new = rl_relative_gap_1_g2_new > 0
iRLOF_2_g2_new = rl_relative_gap_2_g2_new > 0
period_class_g2_new = m3_200_g2_new.hist.period_days

rl_relative_gap_1_g10_new = m3_200_g10_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g10_new = m3_200_g10_new.hist.rl_relative_overflow_2
star_1_mass_g10_new = m3_200_g10_new.hist.star_1_mass
star_2_mass_g10_new = m3_200_g10_new.hist.star_2_mass
iRLOF_1_g10_new = rl_relative_gap_1_g10_new > 0
iRLOF_2_g10_new = rl_relative_gap_2_g10_new > 0
period_class_g10_new = m3_200_g10_new.hist.period_days

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

surf_avg_omega_1_class = m_200.hist.surf_avg_omega

surf_avg_omega_1_pos = m_200_newtides.hist.surf_avg_omega

surf_avg_omega_2_class = m2_200.hist.surf_avg_omega

surf_avg_omega_2_pos = m2_200_newtides.hist.surf_avg_omega

star_age_pos = m2_200_newtides.hist.star_age
star_age_class = m2_200.hist.star_age

# p1.log_fold='LOGS1'

# p1.loadProfile(num=-1)

# p=mp.plot()


rl_relative_gap_1_posydon = m3_200_newtides.hist.rl_relative_overflow_1
rl_relative_gap_2_posydon = m3_200_newtides.hist.rl_relative_overflow_2

age_class = m3_200.hist.age
age_posydon = m3_200_newtides.hist.age

lg_t_sync_2_class = m3_200.hist.lg_t_sync_2
lg_t_sync_2_posydon = m3_200_newtides.hist.lg_t_sync_2

lg_t_sync_1_class = m3_200.hist.lg_t_sync_1
lg_t_sync_1_posydon = m3_200_newtides.hist.lg_t_sync_1

iRLOF_1_posydon = rl_relative_gap_1_posydon > 0
iRLOF_2_posydon = rl_relative_gap_2_posydon > 0

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

pp1_all_panel70 = PdfPages(paths.figures / 'p_q_70days.pdf')

fig = plt.figure(figsize=(10, 10))

plt.title('$\it{P}_\mathrm{ini}$ = 70 [days]', fontsize=30)

# plt.plot(star_2_mass/star_1_mass,period_class,color='k',linestyle='-',label='MESA default, $\gamma$ = 0',lw=2)
# plt.plot(star_2_mass[iRLOF_1]/star_1_mass[iRLOF_1],period_class[iRLOF_1], lw=7, c='k')


if any(iRLOF_2) == True:

    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    indx2 = list(star_2_mass).index(min(star_2_mass[iRLOF_2]))

    print('(len(iRLOF_1) > 0):', indx1)
    print('(len(iRLOF_2) > 0):', indx2)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:indx2] / star_1_mass[indx1:indx2], period_class[indx1:indx2], lw=7, c='k')
    plt.plot(star_2_mass[indx2] / star_1_mass[indx2], period_class[indx2], marker='o', c='k', mfc='k', ms=25)


else:
    indx1 = list(star_2_mass).index(min(star_2_mass[iRLOF_1]))
    print('(len(iRLOF_1) > 0):', indx1)

    plt.plot(star_2_mass[0:indx1] / star_1_mass[0:indx1], period_class[0:indx1], color='k', linestyle='-',
             label='MESA default, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass[indx1:-1] / star_1_mass[indx1:-1], period_class[indx1:-1], lw=7, c='k')
    plt.plot(star_2_mass[-1] / star_1_mass[-1], period_class[-1], marker='o', c='k', ms=25, mfc='None')

if any(iRLOF_2_g1_new) == True:

    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    indx2 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_2_g1_new]))

    print('(len(iRLOF_1_g1_new) > 0):', indx1)
    print('(len(iRLOF_2_g1_new) > 0):', indx2)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:indx2] / star_1_mass_g1_new[indx1:indx2], period_class_g1_new[indx1:indx2], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[indx2] / star_1_mass_g1_new[indx2], period_class_g1_new[indx2], marker='o', c='orange',
             mfc='orange', ms=25)


else:
    indx1 = list(star_2_mass_g1_new).index(min(star_2_mass_g1_new[iRLOF_1_g1_new]))
    print('(len(iRLOF_1_g1_new) > 0):', indx1)

    plt.plot(star_2_mass_g1_new[0:indx1] / star_1_mass_g1_new[0:indx1], period_class_g1_new[0:indx1], color='orange',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 1', lw=2)
    plt.plot(star_2_mass_g1_new[indx1:-1] / star_1_mass_g1_new[indx1:-1], period_class_g1_new[indx1:-1], lw=7,
             c='orange')
    plt.plot(star_2_mass_g1_new[-2] / star_1_mass_g1_new[-2], period_class_g1_new[-2], marker='o', c='orange', ms=25,
             mfc='None')

if any(iRLOF_2_g2_new) == True:

    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    indx2 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_2_g2_new]))

    print('(len(iRLOF_1_g2_new) > 0):', indx1)
    print('(len(iRLOF_2_g2_new) > 0):', indx2)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:indx2] / star_1_mass_g2_new[indx1:indx2], period_class_g2_new[indx1:indx2], lw=7,
             c='blue')
    plt.plot(star_2_mass_g2_new[indx2] / star_1_mass_g2_new[indx2], period_class_g2_new[indx2], marker='o', c='blue',
             mfc='blue', ms=25)


else:
    indx1 = list(star_2_mass_g2_new).index(min(star_2_mass_g2_new[iRLOF_1_g2_new]))
    print('(len(iRLOF_1_g2_new) > 0):', indx1)

    plt.plot(star_2_mass_g2_new[0:indx1] / star_1_mass_g2_new[0:indx1], period_class_g2_new[0:indx1], color='blue',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 2', lw=2)
    plt.plot(star_2_mass_g2_new[indx1:-1] / star_1_mass_g2_new[indx1:-1], period_class_g2_new[indx1:-1], lw=7, c='blue')
    plt.plot(star_2_mass_g2_new[-1] / star_1_mass_g2_new[-1], period_class_g2_new[-1], marker='o', c='blue', ms=25,
             mfc='None')

if any(iRLOF_2_g10_new) == True:

    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    indx2 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_2_g10_new]))

    print('(len(iRLOF_1_g10_new) > 0):', indx1)
    print('(len(iRLOF_2_g10_new) > 0):', indx2)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:indx2] / star_1_mass_g10_new[indx1:indx2], period_class_g10_new[indx1:indx2],
             lw=7, c='red')
    plt.plot(star_2_mass_g10_new[indx2] / star_1_mass[indx2], period_class_g10_new[indx2], marker='o', c='red',
             mfc='red', ms=25)


else:
    indx1 = list(star_2_mass_g10_new).index(min(star_2_mass_g10_new[iRLOF_1_g10_new]))
    print('(len(iRLOF_1_g10_new) > 0):', indx1)

    plt.plot(star_2_mass_g10_new[0:indx1] / star_1_mass_g10_new[0:indx1], period_class_g10_new[0:indx1], color='red',
             linestyle='-', label='MESA default, $\it{\gamma}$ = 10', lw=2)
    plt.plot(star_2_mass_g10_new[indx1:-1] / star_1_mass_g10_new[indx1:-1], period_class_g10_new[indx1:-1], lw=7,
             c='red')
    plt.plot(star_2_mass_g10_new[-1] / star_1_mass_g10_new[-1], period_class_g10_new[-1], marker='o', c='red', ms=25,
             mfc='None')

if any(iRLOF_2_posydon) == True:

    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    indx2 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_2_posydon]))

    print('(len(iRLOF_1_posydon) > 0):', indx1)
    print('(len(iRLOF_2_posydon) > 0):', indx2)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:indx2] / star_1_mass_posydon[indx1:indx2], period_posydon[indx1:indx2], lw=7,
             c='green')
    plt.plot(star_2_mass_posydon[indx2] / star_1_mass_posydon[indx2], period_posydon[indx2], marker='o', c='green',
             mfc='green', ms=25)


else:
    indx1 = list(star_2_mass_posydon).index(min(star_2_mass_posydon[iRLOF_1_posydon]))
    print('(len(iRLOF_1_posydon) > 0):', indx1)

    plt.plot(star_2_mass_posydon[0:indx1] / star_1_mass_posydon[0:indx1], period_posydon[0:indx1], color='green',
             linestyle='-', label='POSYDON, $\it{\gamma}$ = 0', lw=2)
    plt.plot(star_2_mass_posydon[indx1:-1] / star_1_mass_posydon[indx1:-1], period_posydon[indx1:-1], lw=7, c='green')
    plt.plot(star_2_mass_posydon[-1] / star_1_mass_posydon[-1], period_posydon[-1], marker='o', c='green', ms=25,
             mfc='None')

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.xlim([5e6,7.5e6])
# plt.ylim([0.9,1.4]) # 1.25 days

# plt.ylim([1.9,7.5]) # 5 days
# plt.ylim([1,5]) # 3 days
# plt.ylim([0,17]) # 10 days
plt.ylim([0, 125])  # 70 days

# plt.ylim([0,90]) # 50 days


plt.xlabel('$\it{Q}$', fontsize=25)
plt.ylabel('Period [days]', fontsize=25)
plt.legend(loc=2, fontsize=20)

##plt.xlim([0.3,1.8])
###plt.xlim([0.3,1.5])
# plt.xlim([0.64,0.72])
# plt.xlim([0.64,0.72]) # 1.25 days
# plt.xlim([0.25,1.8]) # 5 days
# plt.xlim([0.3,2]) # 3 days

plt.xlim([0.3, 1.8])  # 10 days
# plt.xlim([0.3,1.8]) # 70 days

plt.savefig(pp1_all_panel70, format='pdf')

pp1_all_panel70.close()
