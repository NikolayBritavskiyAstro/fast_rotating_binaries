import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np 
from showyourwork.paths import user as Paths
paths = Paths()
import os
import mesaPlot as mp

if os.path.exists(os.path.join(paths.data,'acceleration/HD191495_mesa/LOGS1/history.data')):
	pass
else:
	os.system(f'python {os.path.join(paths.scripts / "unzip_MESA_output.py")}')

plt.style.use(paths.scripts / "matplotlibrc")


m_200_191=mp.MESA()
m2_200_191=mp.MESA()
m3_200_191=mp.MESA()


m_200_191_newtides=mp.MESA()
m2_200_191_newtides=mp.MESA()
m3_200_191_newtides=mp.MESA()


m_200_nowind_191=mp.MESA()
m2_200_nowind_191=mp.MESA()
m3_200_nowind_191=mp.MESA()


m_200_46=mp.MESA()
m2_200_46=mp.MESA()
m3_200_46=mp.MESA()

m_200_46_newtides=mp.MESA()
m2_200_46_newtides=mp.MESA()
m3_200_46_newtides=mp.MESA()


m_200_nowind_46=mp.MESA()
m2_200_nowind_46=mp.MESA()
m3_200_nowind_46=mp.MESA()


m_200_25=mp.MESA()
m2_200_25=mp.MESA()
m3_200_25=mp.MESA()

m_200_25_newtides=mp.MESA()
m2_200_25_newtides=mp.MESA()
m3_200_25_newtides=mp.MESA()


m_200_nowind_25=mp.MESA()
m2_200_nowind_25=mp.MESA()
m3_200_nowind_25=mp.MESA()



m_200_191.log_fold=os.path.join(paths.data,'acceleration/HD191495_mesa/LOGS1')
m_200_191.loadHistory()

m2_200_191.log_fold=os.path.join(paths.data,'acceleration/HD191495_mesa/LOGS2')
m2_200_191.loadHistory()


m3_200_191.log_fold=os.path.join(paths.data,'acceleration/HD191495_mesa/LOGS3')
m3_200_191.loadHistory()


m_200_191_newtides.log_fold=os.path.join(paths.data,'acceleration/HD191495_pos/LOGS1')
m_200_191_newtides.loadHistory()

m2_200_191_newtides.log_fold=os.path.join(paths.data,'acceleration/HD191495_pos/LOGS2')
m2_200_191_newtides.loadHistory()


m3_200_191_newtides.log_fold=os.path.join(paths.data,'acceleration/HD191495_pos/LOGS3')
m3_200_191_newtides.loadHistory()



m_200_nowind_191.log_fold=os.path.join(paths.data,'acceleration/HD191495_notides/LOGS1')
m_200_nowind_191.loadHistory()

m2_200_nowind_191.log_fold=os.path.join(paths.data,'acceleration/HD191495_notides/LOGS2')
m2_200_nowind_191.loadHistory()


m3_200_nowind_191.log_fold=os.path.join(paths.data,'acceleration/HD191495_notides/LOGS3')
m3_200_nowind_191.loadHistory()




m_200_46_newtides.log_fold=os.path.join(paths.data,'acceleration/HD46485_pos/LOGS1')
m_200_46_newtides.loadHistory()

m2_200_46_newtides.log_fold=os.path.join(paths.data,'acceleration/HD46485_pos/LOGS2')
m2_200_46_newtides.loadHistory()


m3_200_46_newtides.log_fold=os.path.join(paths.data,'acceleration/HD46485_pos/LOGS3')
m3_200_46_newtides.loadHistory()



m_200_46.log_fold=os.path.join(paths.data,'acceleration/HD46485_mesa/LOGS1')
m_200_46.loadHistory()

m2_200_46.log_fold=os.path.join(paths.data,'acceleration/HD46485_mesa/LOGS2')
m2_200_46.loadHistory()


m3_200_46.log_fold=os.path.join(paths.data,'acceleration/HD46485_mesa/LOGS3')
m3_200_46.loadHistory()



m_200_nowind_46.log_fold=os.path.join(paths.data,'acceleration/HD46485_notides/LOGS1')
m_200_nowind_46.loadHistory()

m2_200_nowind_46.log_fold=os.path.join(paths.data,'acceleration/HD46485_notides/LOGS2')
m2_200_nowind_46.loadHistory()


m3_200_nowind_46.log_fold=os.path.join(paths.data,'acceleration/HD46485_notides/LOGS3')
m3_200_nowind_46.loadHistory()



m_200_25.log_fold=os.path.join(paths.data,'acceleration/HD25631_mesa/LOGS1')
m_200_25.loadHistory()

m2_200_25.log_fold=os.path.join(paths.data,'acceleration/HD25631_mesa/LOGS2')
m2_200_25.loadHistory()

m3_200_25.log_fold=os.path.join(paths.data,'acceleration/HD25631_mesa/LOGS3')
m3_200_25.loadHistory()

m_200_25_newtides.log_fold=os.path.join(paths.data,'acceleration/HD25631_pos/LOGS1')
m_200_25_newtides.loadHistory()

m2_200_25_newtides.log_fold=os.path.join(paths.data,'acceleration/HD25631_pos/LOGS2')
m2_200_25_newtides.loadHistory()


m3_200_25_newtides.log_fold=os.path.join(paths.data,'acceleration/HD25631_pos/LOGS3')
m3_200_25_newtides.loadHistory()


m_200_nowind_25.log_fold=os.path.join(paths.data,'acceleration/HD25631_notides/LOGS1')
m_200_nowind_25.loadHistory()

m2_200_nowind_25.log_fold=os.path.join(paths.data,'acceleration/HD25631_notides/LOGS2')
m2_200_nowind_25.loadHistory()


m3_200_nowind_25.log_fold=os.path.join(paths.data,'acceleration/HD25631_notides/LOGS3')
m3_200_nowind_25.loadHistory()


rsun=696000 #km


star_age_200_191=m_200_191.hist.star_age
period_days_191=m3_200_191.hist.period_days
surf_avg_vtor_200_191=m_200_191.hist.surf_avg_v_rot
surf_avg_omega200_191=m_200_191.hist.surf_avg_omega
surf_avg_omega_crit200_191=m_200_191.hist.surf_avg_omega_crit
star_1_radius200_191=m3_200_191.hist.star_1_radius
star_1_J_orb_200_191=m3_200_191.hist.J_orb
star_1_J_spin_200_191=m3_200_191.hist.J_spin_1
star_2_J_spin_200_191=m3_200_191.hist.J_spin_2
age_200_191=m3_200_191.hist.age
jtotal_200_191=m3_200_191.hist.J_total
log_total_angular_momentum_200_191=m_200_191.hist.log_total_angular_momentum
surf_avg_j_rot_200_191=m_200_191.hist.surf_avg_j_rot
center_h1_200_191=m_200_191.hist.center_h1
iRLOF_200_191 = np.argmax(center_h1_200_191 < 1e-2)
stable_tides=[3*el for el in np.add(star_1_J_spin_200_191,star_2_J_spin_200_191)]

log_Teff_1_200_191=m_200_191.hist.log_Teff
LOGL_1_200_46_191= m_200_191.hist.log_L





star_age_200_nowind_191=m_200_nowind_191.hist.star_age

surf_avg_vtor_200_nowind_191=m_200_nowind_191.hist.surf_avg_v_rot
surf_avg_omega200_nowind_191=m_200_nowind_191.hist.surf_avg_omega
surf_avg_omega_crit200_191_nowind=m_200_nowind_191.hist.surf_avg_omega_crit
star_1_radius200_nowind_191=m3_200_nowind_191.hist.star_1_radius
star_1_J_orb_200_nowind_191=m3_200_nowind_191.hist.J_orb
star_1_J_spin_200_nowind_191=m3_200_nowind_191.hist.J_spin_1
star_2_J_spin_200_nowind_191=m3_200_nowind_191.hist.J_spin_2
age_200_nowind_191=m3_200_nowind_191.hist.age
jtotal_200_nowind_191=m3_200_nowind_191.hist.J_total
log_total_angular_momentum_200_nowind_191=m_200_nowind_191.hist.log_total_angular_momentum
surf_avg_j_rot_200_nowind_191=m_200_nowind_191.hist.surf_avg_j_rot
center_h1_200_nowind_191=m_200_nowind_191.hist.center_h1









star_age_200_46=m_200_46.hist.star_age
period_days_46=m3_200_46.hist.period_days

surf_avg_vtor_200_46=m_200_46.hist.surf_avg_v_rot
surf_avg_omega200_46=m_200_46.hist.surf_avg_omega
surf_avg_omega_crit200_46=m_200_46.hist.surf_avg_omega_crit
star_1_radius200_46=m3_200_46.hist.star_1_radius
star_1_J_orb_200_46=m3_200_46.hist.J_orb
star_1_J_spin_200_46=m3_200_46.hist.J_spin_1
star_2_J_spin_200_46=m3_200_46.hist.J_spin_2
age_200_46=m3_200_46.hist.age
jtotal_200_46=m3_200_46.hist.J_total
log_total_angular_momentum_200_46=m_200_46.hist.log_total_angular_momentum
surf_avg_j_rot_200_46=m_200_46.hist.surf_avg_j_rot
center_h1_200_46=m_200_46.hist.center_h1
iRLOF_200_46 = np.argmax(center_h1_200_46 < 1e-2)


stable_tides=[3*el for el in np.add(star_1_J_spin_200_46,star_2_J_spin_200_46)]

star_age_200_nowind_46=m_200_nowind_46.hist.star_age

surf_avg_vtor_200_nowind_46=m_200_nowind_46.hist.surf_avg_v_rot
surf_avg_omega200_nowind_46=m_200_nowind_46.hist.surf_avg_omega
surf_avg_omega_crit200_46_nowind=m_200_nowind_46.hist.surf_avg_omega_crit
star_1_radius200_nowind_46=m3_200_nowind_46.hist.star_1_radius
star_1_J_orb_200_nowind_46=m3_200_nowind_46.hist.J_orb
star_1_J_spin_200_nowind_46=m3_200_nowind_46.hist.J_spin_1
star_2_J_spin_200_nowind_46=m3_200_nowind_46.hist.J_spin_2
age_200_nowind_46=m3_200_nowind_46.hist.age
jtotal_200_nowind_46=m3_200_nowind_46.hist.J_total
log_total_angular_momentum_200_nowind_46=m_200_nowind_46.hist.log_total_angular_momentum
surf_avg_j_rot_200_nowind_46=m_200_nowind_46.hist.surf_avg_j_rot
center_h1_200_nowind_46=m_200_nowind_46.hist.center_h1




star_age_200_25=m_200_25.hist.star_age
period_days_25=m3_200_25.hist.period_days

surf_avg_vtor_200_25=m_200_25.hist.surf_avg_v_rot
surf_avg_omega200_25=m_200_25.hist.surf_avg_omega
surf_avg_omega_crit200_25=m_200_25.hist.surf_avg_omega_crit
star_1_radius200_25=m3_200_25.hist.star_1_radius
star_1_J_orb_200_25=m3_200_25.hist.J_orb
star_1_J_spin_200_25=m3_200_25.hist.J_spin_1
star_2_J_spin_200_25=m3_200_25.hist.J_spin_2
age_200_25=m3_200_25.hist.age
jtotal_200_25=m3_200_25.hist.J_total
log_total_angular_momentum_200_25=m_200_25.hist.log_total_angular_momentum
surf_avg_j_rot_200_25=m_200_25.hist.surf_avg_j_rot
center_h1_200_25=m_200_25.hist.center_h1
iRLOF_200_25 = np.argmax(center_h1_200_25 < 1e-2)






stable_tides=[3*el for el in np.add(star_1_J_spin_200_25,star_2_J_spin_200_25)]

star_age_200_nowind_25=m_200_nowind_25.hist.star_age


surf_avg_vtor_200_nowind_25=m_200_nowind_25.hist.surf_avg_v_rot
surf_avg_omega200_nowind_25=m_200_nowind_25.hist.surf_avg_omega
surf_avg_omega_crit200_25_nowind=m_200_nowind_25.hist.surf_avg_omega_crit
star_1_radius200_nowind_25=m3_200_nowind_25.hist.star_1_radius
star_1_J_orb_200_nowind_25=m3_200_nowind_25.hist.J_orb
star_1_J_spin_200_nowind_25=m3_200_nowind_25.hist.J_spin_1
star_2_J_spin_200_nowind_25=m3_200_nowind_25.hist.J_spin_2
age_200_nowind_25=m3_200_nowind_25.hist.age
jtotal_200_nowind_25=m3_200_nowind_25.hist.J_total
log_total_angular_momentum_200_nowind_25=m_200_nowind_25.hist.log_total_angular_momentum
surf_avg_j_rot_200_nowind_25=m_200_nowind_25.hist.surf_avg_j_rot
center_h1_200_nowind_25=m_200_nowind_25.hist.center_h1


log_Teff_1_200_46_notides=m_200_46_newtides.hist.log_Teff
LOGL_1_200_46_notides=m_200_46_newtides.hist.log_L



star_age_200_46_newtides=m_200_46_newtides.hist.star_age
surf_avg_vtor_200_46_newtides=m_200_46_newtides.hist.surf_avg_v_rot
surf_avg_omega200_46_newtides=m_200_46_newtides.hist.surf_avg_omega
surf_avg_omega_crit200_46_newtides=m_200_46_newtides.hist.surf_avg_omega_crit

star_1_radius200_46_newtides=m3_200_46_newtides.hist.star_1_radius
star_1_J_orb_200_46_newtides=m3_200_46_newtides.hist.J_orb
star_1_J_spin_200_46_newtides=m3_200_46_newtides.hist.J_spin_1
star_2_J_spin_200_46_newtides=m3_200_46_newtides.hist.J_spin_2
age_200_46_newtides=m3_200_46_newtides.hist.age
jtotal_200_46_newtides=m3_200_46_newtides.hist.J_total
log_total_angular_momentum_200_46_newtides=m_200_46_newtides.hist.log_total_angular_momentum
surf_avg_j_rot_200_46_newtides=m_200_46_newtides.hist.surf_avg_j_rot
center_h1_200_46_newtides=m_200_46_newtides.hist.center_h1
stable_tides_newtides=[3*el for el in np.add(star_1_J_spin_200_46_newtides,star_2_J_spin_200_46_newtides)]


iRLOF_200_46_newtides = np.argmax(center_h1_200_46_newtides < 1e-2)





star_age_200_25_newtides=m_200_25_newtides.hist.star_age
surf_avg_vtor_200_25_newtides=m_200_25_newtides.hist.surf_avg_v_rot
surf_avg_omega200_25_newtides=m_200_25_newtides.hist.surf_avg_omega
surf_avg_omega_crit200_25_newtides=m_200_25_newtides.hist.surf_avg_omega_crit
star_1_radius200_25_newtides=m3_200_25_newtides.hist.star_1_radius
star_1_J_orb_200_25_newtides=m3_200_25_newtides.hist.J_orb
star_1_J_spin_200_25_newtides=m3_200_25_newtides.hist.J_spin_1
star_2_J_spin_200_25_newtides=m3_200_25_newtides.hist.J_spin_2
age_200_25_newtides=m3_200_25_newtides.hist.age
jtotal_200_25_newtides=m3_200_25_newtides.hist.J_total
log_total_angular_momentum_200_25_newtides=m_200_25_newtides.hist.log_total_angular_momentum
surf_avg_j_rot_200_25_newtides=m_200_25_newtides.hist.surf_avg_j_rot
center_h1_200_25_newtides=m_200_25_newtides.hist.center_h1
stable_tides_newtides=[3*el for el in np.add(star_1_J_spin_200_25_newtides,star_2_J_spin_200_25_newtides)]

iRLOF_200_25_newtides = np.argmax(center_h1_200_25_newtides < 1e-2)




star_age_200_191_newtides=m_200_191_newtides.hist.star_age
surf_avg_vtor_200_191_newtides=m_200_191_newtides.hist.surf_avg_v_rot
surf_avg_omega200_191_newtides=m_200_191_newtides.hist.surf_avg_omega
surf_avg_omega_crit200_191_newtides=m_200_191_newtides.hist.surf_avg_omega_crit
star_1_radius200_191_newtides=m3_200_191_newtides.hist.star_1_radius
star_1_J_orb_200_191_newtides=m3_200_191_newtides.hist.J_orb
star_1_J_spin_200_191_newtides=m3_200_191_newtides.hist.J_spin_1
star_2_J_spin_200_191_newtides=m3_200_191_newtides.hist.J_spin_2
age_200_191_newtides=m3_200_191_newtides.hist.age
jtotal_200_191_newtides=m3_200_191_newtides.hist.J_total
log_total_angular_momentum_200_191_newtides=m_200_191_newtides.hist.log_total_angular_momentum
surf_avg_j_rot_200_191_newtides=m_200_191_newtides.hist.surf_avg_j_rot
center_h1_200_191_newtides=m_200_191_newtides.hist.center_h1

log_Teff_1_200_191_newtides=m_200_191_newtides.hist.log_Teff
LOGL_1_200_46_191_newtides= m_200_191_newtides.hist.log_L



stable_tides=[3*el for el in np.add(star_1_J_spin_200_191_newtides,star_2_J_spin_200_191_newtides)]

iRLOF_200_191_newtides = np.argmax(center_h1_200_191_newtides < 1e-2)




pp1_all = PdfPages(paths.figures / 'age_vsurf_acceleration_panel.pdf') 
fig=plt.figure(figsize=(10, 10))

ax1=fig.add_subplot(311)
#ax1.title.set_text('HD46485, $v^{ini}_{rot}$=1 $km~s^{-1}$, M1$_{ini}$=24, M2$_{ini}$=1, p$_{ini}$=6.9d',fontsize=24)
ax1.set_title('HD46485, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=1 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=24, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=6.9 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})

plt.plot(star_age_200_46/1e6, surf_avg_vtor_200_46,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_46_newtides/1e6, surf_avg_vtor_200_46_newtides, linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_nowind_46/1e6, surf_avg_vtor_200_nowind_46,linestyle='-.',color='r',label='No tides')


v_sync_191=2*3.1415926*star_1_radius200_191*rsun/(period_days_191*24*3600)
v_sync_46=2*3.1415926*star_1_radius200_46*rsun/(period_days_46*24*3600)
v_sync_25=2*3.1415926*star_1_radius200_25*rsun/(period_days_25*24*3600)


plt.plot(age_200_46/1e6, v_sync_46,linestyle='--',color='gray',label='$\it{v}_\mathrm{sync}$')




#plt.plot(star_age_200_46_newtides_res05/1e6, surf_avg_vtor_200_46_newtides_res05, #linestyle='--',color='b',label='POSYDON, res=05',lw=4)
#plt.plot(star_age_200_46_notides_res05/1e6, surf_avg_vtor_200_46_notides_res05,linestyle='-.',color='r',label='No tides, res=05',lw=4)


#plt.plot(star_age_200_nowind_46/1e6, surf_avg_vtor_200_nowind_46,linestyle='-',color='black',label='no wind',lw=1)


#plt.ylabel('Average rotational velocity [$km~s^{-1}$]',fontsize=24)
#plt.xlabel('Star age [Myrs]',fontsize=24)
plt.legend(loc=2,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)



ax2=fig.add_subplot(312)



ax2.set_title('HD191495, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=1 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=15, $\it{M}_\mathrm{2,ini}$=1.5, $\it{P}_\mathrm{ini}$=3.6 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(star_age_200_191/1e6, surf_avg_vtor_200_191,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_191_newtides/1e6, surf_avg_vtor_200_191_newtides,linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_nowind_191/1e6, surf_avg_vtor_200_nowind_191,linestyle='-.',color='r',label='No tides')
plt.plot(age_200_191/1e6, v_sync_191,linestyle='--',color='gray',label='$\it{v}_\mathrm{sync}$')


#plt.plot(star_age_200_nowind_191/1e6, surf_avg_vtor_200_nowind_191,linestyle='-',label='no wind',color='b',lw=)
#plt.plot([min(star_age_200_25),max(star_age_200_25)],[221,221],color='g',linestyle='--',label='HD 25631, vsini=221')
#plt.plot([min(star_age_1),max(star_age_1)],[224,224],linestyle='--',label=' veq*sin(80)/sin(80)=224')
#plt.plot([min(star_age_200_25),max(star_age_200_25)],[201,201],color='r',linestyle='--',label='HD 191495, vsini=201')
#plt.plot([min(star_age_200_25),max(star_age_200_25)],[334,334],color='m',linestyle='--',label='HD 46485, vsini=334')
plt.xlim([-0.5,14])
#plt.ylabel('Average rotational velocity [$km~s^{-1}$]',fontsize=24)
plt.legend(loc=2,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('Surface rotational velocity [$km~s^{-1}$]',fontsize=24)





ax3=fig.add_subplot(313)


ax3.set_title('HD25631, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=1 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=7, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=5.2 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(star_age_200_25/1e6, surf_avg_vtor_200_25,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_25_newtides/1e6, surf_avg_vtor_200_25_newtides,linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_nowind_25/1e6, surf_avg_vtor_200_nowind_25,linestyle='-.',color='r',label='No tides')

plt.plot(age_200_25/1e6, v_sync_25,linestyle='--',color='gray',label='$\it{v}_\mathrm{sync}$')

#plt.plot(star_age_200_nowind_25/1e6, surf_avg_vtor_200_nowind_25,linestyle='-',color='r',label='no wind',lw=1)


#plt.xlabel('Star age [Myrs]',fontsize=24)
plt.legend(loc=2,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Star age [Myrs]',fontsize=24)











plt.savefig(pp1_all, format='pdf')

pp1_all.close()


param=24*3600


pp2_all = PdfPages(paths.figures / 'age_omega_comparison_panel.pdf') 


fig=plt.figure(figsize=(10, 10))

ax1=fig.add_subplot(311)
#ax1.title.set_text('HD46485, $v^{ini}_{rot}$=1 $km~s^{-1}$, M1$_{ini}$=24, M2$_{ini}$=1, p=6.9d',fontsize=24)
ax1.set_title('HD46485, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=1 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=24, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=6.9 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})

plt.plot(star_age_200_46/1e6, surf_avg_omega200_46*param,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_46_newtides/1e6, surf_avg_omega200_46_newtides*param, linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_nowind_46/1e6, surf_avg_omega200_nowind_46*param,linestyle='-.',color='r',label='No tides')


#plt.plot(star_age_200_46/1e6, surf_avg_omega_crit200_46,'--',color='gray',label='$\omega_{crit}$')
plt.plot(age_200_46/1e6, v_sync_46/star_1_radius200_46/rsun*param,'-.',color='gray',label='$\it{\omega}_\mathrm{sync}$')






#plt.plot(star_age_200_46_newtides_res05/1e6, surf_avg_vtor_200_46_newtides_res05, #linestyle='--',color='b',label='POSYDON, res=05',lw=4)
#plt.plot(star_age_200_46_notides_res05/1e6, surf_avg_vtor_200_46_notides_res05,linestyle='-.',color='r',label='No tides, res=05',lw=4)


#plt.plot(star_age_200_nowind_46/1e6, surf_avg_vtor_200_nowind_46,linestyle='-',color='black',label='no wind',lw=1)


#plt.ylabel('Average rotational velocity [$km~s^{-1}$]',fontsize=24)
#plt.xlabel('Star age [Myrs]',fontsize=24)
plt.legend(loc=2,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)



ax2=fig.add_subplot(312)



ax2.set_title('HD191495, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=1 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=15, $\it{M}_\mathrm{2,ini}$=1.5, $\it{P}_\mathrm{ini}$=3.6 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(star_age_200_191/1e6, surf_avg_omega200_191*param,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_191_newtides/1e6, surf_avg_omega200_191_newtides*param,linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_nowind_191/1e6, surf_avg_omega200_nowind_191*param,linestyle='-.',color='r',label='No tides')


#plt.plot(star_age_200_191/1e6, surf_avg_omega_crit200_191,'--',color='gray',label='$\omega_{crit}$')
plt.plot(age_200_191/1e6, v_sync_191/star_1_radius200_191/rsun*param,'-.',color='gray',label='$\it{\omega}_\mathrm{sync}$')






#plt.plot(star_age_200_nowind_191/1e6, surf_avg_vtor_200_nowind_191,linestyle='-',label='no wind',color='b',lw=)
#plt.plot([min(star_age_200_25),max(star_age_200_25)],[221,221],color='g',linestyle='--',label='HD 25631, vsini=221')
#plt.plot([min(star_age_1),max(star_age_1)],[224,224],linestyle='--',label=' veq*sin(80)/sin(80)=224')
#plt.plot([min(star_age_200_25),max(star_age_200_25)],[201,201],color='r',linestyle='--',label='HD 191495, vsini=201')
#plt.plot([min(star_age_200_25),max(star_age_200_25)],[334,334],color='m',linestyle='--',label='HD 46485, vsini=334')
plt.xlim([-0.5,14])
#plt.ylabel('Average rotational velocity [$km~s^{-1}$]',fontsize=24)
plt.legend(loc=2,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('Surface average angular velocity [rad/day]',fontsize=24)





ax3=fig.add_subplot(313)


ax3.set_title('HD25631, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=1 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=7, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=5.2 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(star_age_200_25/1e6, surf_avg_omega200_25*param,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_25_newtides/1e6, surf_avg_omega200_25_newtides*param,linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_nowind_25/1e6, surf_avg_omega200_nowind_25*param,linestyle='-.',color='r',label='No tides')



#plt.plot(star_age_200_25/1e6, surf_avg_omega_crit200_25,'--',color='gray',label='$\omega_{crit}$')
plt.plot(age_200_25/1e6, v_sync_25/star_1_radius200_25/rsun*param,'-.',color='gray',label='$\it{\omega}_\mathrm{sync}$')




#plt.plot(star_age_200_nowind_25/1e6, surf_avg_vtor_200_nowind_25,linestyle='-',color='r',label='no wind',lw=1)


#plt.xlabel('Star age [Myrs]',fontsize=24)
plt.legend(loc=2,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Star age [Myrs]',fontsize=24)


plt.savefig(pp2_all, format='pdf')

pp2_all.close()




