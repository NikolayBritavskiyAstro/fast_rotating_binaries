import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np 
from showyourwork.paths import user as Paths
paths = Paths()
import os
import mesaPlot as mp

plt.style.use(paths.scripts / "matplotlibrc")


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



m3_p3_g1_new=mp.MESA()
m3_p5_g1_new=mp.MESA()
m3_p10_g1_new=mp.MESA()
m3_p50_g1_new=mp.MESA()
m3_p70_g1_new=mp.MESA()



m2_p3_g1_new=mp.MESA()
m2_p5_g1_new=mp.MESA()
m2_p10_g1_new=mp.MESA()
m2_p50_g1_new=mp.MESA()
m2_p70_g1_new=mp.MESA()





m3_p3_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p3/LOGS3')
m3_p3_g1_new.loadHistory()

m3_p5_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p5/LOGS3')
m3_p5_g1_new.loadHistory()

m3_p10_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p10/LOGS3')
m3_p10_g1_new.loadHistory()


m3_p50_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p50/LOGS3')
m3_p50_g1_new.loadHistory()


m3_p70_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p70/LOGS3')
m3_p70_g1_new.loadHistory()





m2_p3_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p3/LOGS2')
m2_p3_g1_new.loadHistory()

m2_p5_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p5/LOGS2')
m2_p5_g1_new.loadHistory()

m2_p10_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p10/LOGS2')
m2_p10_g1_new.loadHistory()


m2_p50_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p50/LOGS2')
m2_p50_g1_new.loadHistory()


m2_p70_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p70/LOGS2')
m2_p70_g1_new.loadHistory()






'''
rl_relative_gap_1_g1_new= m3_200_g1_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g1_new= m3_200_g1_new.hist.rl_relative_overflow_2
'''
indx_p3=[]
indx_p5=[]
indx_p10=[]
indx_p50=[]
indx_p70=[]

age_p3_g1_new=m3_p3_g1_new.hist.age
delta_p3_g1_new=m3_p3_g1_new.hist.delta_updated
spin2_p3_g1_new=m3_p3_g1_new.hist.J_spin_2


period_days_p3=m3_p3_g1_new.hist.period_days
omega_sync_p3=2*3.1415926/(period_days_p3) #rad/days



star_age_p3_g1_new=m2_p3_g1_new.hist.star_age
surf_avg_omega_div_omega_crit_p3_g1_new=m2_p3_g1_new.hist.surf_avg_omega_div_omega_crit
surf_avg_omega_crit_p3_g1_new=m2_p3_g1_new.hist.surf_avg_omega_crit

age_p5_g1_new=m3_p5_g1_new.hist.age
delta_p5_g1_new=m3_p5_g1_new.hist.delta_updated
spin2_p5_g1_new=m3_p5_g1_new.hist.J_spin_2

period_days_p5=m3_p5_g1_new.hist.period_days
omega_sync_p5=2*3.1415926/(period_days_p5) #rad/days


star_age_p5_g1_new=m2_p5_g1_new.hist.star_age
surf_avg_omega_div_omega_crit_p5_g1_new=m2_p5_g1_new.hist.surf_avg_omega_div_omega_crit
surf_avg_omega_crit_p5_g1_new=m2_p5_g1_new.hist.surf_avg_omega_crit



age_p10_g1_new=m3_p10_g1_new.hist.age
delta_p10_g1_new=m3_p10_g1_new.hist.delta_updated
spin2_p10_g1_new=m3_p10_g1_new.hist.J_spin_2


star_age_p10_g1_new=m2_p10_g1_new.hist.star_age
surf_avg_omega_div_omega_crit_p10_g1_new=m2_p10_g1_new.hist.surf_avg_omega_div_omega_crit


age_p50_g1_new=m3_p50_g1_new.hist.age
delta_p50_g1_new=m3_p50_g1_new.hist.delta_updated
spin2_p50_g1_new=m3_p50_g1_new.hist.J_spin_2

star_age_p50_g1_new=m2_p50_g1_new.hist.star_age
surf_avg_omega_div_omega_crit_p50_g1_new=m2_p50_g1_new.hist.surf_avg_omega_div_omega_crit


age_p70_g1_new=m3_p70_g1_new.hist.age
delta_p70_g1_new=m3_p70_g1_new.hist.delta_updated
spin2_p70_g1_new=m3_p70_g1_new.hist.J_spin_2
period_days_p70=m3_p70_g1_new.hist.period_days
surf_avg_omega_crit_p70_g1_new=m2_p70_g1_new.hist.surf_avg_omega_crit
omega_sync_p70=2*3.1415926/(period_days_p70) #rad/days


star_age_p70_g1_new=m2_p70_g1_new.hist.star_age
surf_avg_omega_div_omega_crit_p70_g1_new=m2_p70_g1_new.hist.surf_avg_omega_div_omega_crit



for i in range(len(age_p3_g1_new)):
    indx_p3.append(find_nearest(star_age_p3_g1_new, age_p3_g1_new[i]))


for i in range(len(age_p5_g1_new)):
    indx_p5.append(find_nearest(star_age_p5_g1_new, age_p5_g1_new[i]))


for i in range(len(age_p10_g1_new)):
    indx_p10.append(find_nearest(star_age_p10_g1_new, age_p10_g1_new[i]))

for i in range(len(age_p50_g1_new)):
    indx_p50.append(find_nearest(star_age_p50_g1_new, age_p50_g1_new[i]))


for i in range(len(age_p70_g1_new)):
    indx_p70.append(find_nearest(star_age_p70_g1_new, age_p70_g1_new[i]))







'''
star_1_mass_g1_new= m3_200_g1_new.hist.star_1_mass
star_2_mass_g1_new= m3_200_g1_new.hist.star_2_mass
iRLOF_1_g1_new = rl_relative_gap_1_g1_new > 0
iRLOF_2_g1_new = rl_relative_gap_2_g1_new > 0
period_class_g1_new = m3_200_g1_new.hist.period_days

age_200_g2_new=m3_200_g2_new.hist.age
rl_relative_gap_1_g2_new= m3_200_g2_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g2_new= m3_200_g2_new.hist.rl_relative_overflow_2
star_1_mass_g2_new= m3_200_g2_new.hist.star_1_mass
star_2_mass_g2_new= m3_200_g2_new.hist.star_2_mass
iRLOF_1_g2_new = rl_relative_gap_1_g2_new > 0
iRLOF_2_g2_new = rl_relative_gap_2_g2_new > 0
period_class_g2_new = m3_200_g2_new.hist.period_days

age_200_g10_new=m3_200_g10_new.hist.age
rl_relative_gap_1_g10_new= m3_200_g10_new.hist.rl_relative_overflow_1
rl_relative_gap_2_g10_new= m3_200_g10_new.hist.rl_relative_overflow_2
star_1_mass_g10_new= m3_200_g10_new.hist.star_1_mass
star_2_mass_g10_new= m3_200_g10_new.hist.star_2_mass
iRLOF_1_g10_new = rl_relative_gap_1_g10_new > 0
iRLOF_2_g10_new = rl_relative_gap_2_g10_new > 0
period_class_g10_new = m3_200_g10_new.hist.period_days


age_200=m3_200.hist.age
jtotal_200=m3_200.hist.J_total
log_total_angular_momentum_200=m_200.hist.log_total_angular_momentum
surf_avg_j_rot_200=m_200.hist.surf_avg_j_rot
center_h1_200=m_200.hist.center_h1
'''





pp1_period = PdfPages(paths.figures / 'eta_omega.pdf') 



plt.figure(figsize=(10, 10))


plt.title('$\it{M}_\mathrm{don,ini}$ = 30, $\it{M}_\mathrm{acc,ini}$ = 20, $\it{\gamma}$ = 3',fontsize=30)#plt.tick_params(labelsize=18)


print(delta_p3_g1_new)
print(age_p3_g1_new)
param=24*3600

'''
plt.plot(age_p3_g1_new/1e6, 1-delta_p3_g1_new,linestyle='-',label='p = 3d',lw=2)
plt.plot(age_p5_g1_new/1e6, 1-delta_p5_g1_new,linestyle='-',label='p = 5d',lw=2)

plt.plot(age_p10_g1_new/1e6, 1-delta_p10_g1_new,linestyle='-',label='p = 10d',lw=2)
plt.plot(age_p50_g1_new/1e6, 1-delta_p50_g1_new,linestyle='-',label='p = 50d',lw=2)
plt.plot(age_p70_g1_new[:-20]/1e6, 1-delta_p70_g1_new[:-20],linestyle='-',label='p = 70d',lw=2)
'''


'''
plt.scatter(spin2_p3_g1_new, 1-delta_p3_g1_new,label='p = 3d',lw=2)
plt.scatter(spin2_p5_g1_new, 1-delta_p5_g1_new,linestyle='-',label='p = 5 days',lw=2)

plt.scatter(spin2_p10_g1_new, 1-delta_p10_g1_new,linestyle='-',label='p = 10d',lw=2)
plt.scatter(spin2_p50_g1_new, 1-delta_p50_g1_new,linestyle='-',label='p = 50d',lw=2)
plt.scatter(spin2_p70_g1_new[:-20], 1-delta_p70_g1_new[:-20],linestyle='-',label='p = 70d',lw=2)
'''

plt.plot(surf_avg_omega_div_omega_crit_p3_g1_new[indx_p3][:-150], 1-delta_p3_g1_new[:-150],label='$\it{P}_\mathrm{ini}$ = 3 d',lw=3)
#plt.plot(omega_sync_p3/surf_avg_omega_crit_p3_g1_new[indx_p3]/param, 1-delta_p3_g1_new,label='p = 3d ($\omega_{sync}$)',lw=2,linestyle='--')

plt.plot(surf_avg_omega_div_omega_crit_p5_g1_new[indx_p5][:-150], 1-delta_p5_g1_new[:-150],linestyle='-',label='$\it{P}_\mathrm{ini}$ = 5 d',lw=3)
#plt.plot(omega_sync_p5/surf_avg_omega_crit_p5_g1_new[indx_p5]/param, 1-delta_p5_g1_new,label='p = 5d ($\omega_{sync}$)',lw=2,linestyle='--')


plt.plot(surf_avg_omega_div_omega_crit_p10_g1_new[indx_p10], 1-delta_p10_g1_new,linestyle='-',label='$\it{P}_\mathrm{ini}$ = 10 d',lw=3)
plt.plot(surf_avg_omega_div_omega_crit_p50_g1_new[indx_p50][:-30], 1-delta_p50_g1_new[:-30],linestyle='-',label='$\it{P}_\mathrm{ini}$ = 50 d',lw=3)
plt.plot(surf_avg_omega_div_omega_crit_p70_g1_new[indx_p70][:-30], 1-delta_p70_g1_new[:-30],linestyle='-',label='$\it{P}_\mathrm{ini}$ = 70 d',lw=3)

#plt.plot(omega_sync_p70/surf_avg_omega_crit_p70_g1_new[indx_p70]/param, 1-delta_p70_g1_new,label='p = 5d ($\omega_{sync}$)',lw=2,linestyle='--')




#print(star_age_p3_g1_new[indx_p3], surf_avg_omega_div_omega_crit_p3_g1_new[indx_p3])


plt.xlabel('$\it{\omega}/\it{\omega}_\mathrm{crit}$',fontsize=30)


plt.ylabel('$\it{\eta} = \Delta \it{M}_\mathrm{acc}/\Delta \it{M}_\mathrm{don}$',fontsize=30)
#plt.xlabel('Star age [Myrs]',fontsize=30)
plt.legend(loc=1,fontsize=23,frameon=True)
#plt.xlim([5e6,6.5e6])

plt.ylim([-0.02,0.25])

plt.savefig(pp1_period, format='pdf')


pp1_period.close()






