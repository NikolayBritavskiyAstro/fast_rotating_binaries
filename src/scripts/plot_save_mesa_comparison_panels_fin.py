import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np 
from showyourwork.paths import user as Paths
paths = Paths()
import os
import mesaPlot as mp

plt.style.use(paths.scripts / "matplotlibrc")


rsun=696000 #km


m_200_191=mp.MESA()
m3_200_191=mp.MESA()


m_200_191_newtides=mp.MESA()
m3_200_191_newtides=mp.MESA()


m_200_191_newtides_str=mp.MESA()
m3_200_191_newtides_str=mp.MESA()


m_200_46=mp.MESA()
m3_200_46=mp.MESA()

m_200_46_newtides=mp.MESA()
m3_200_46_newtides=mp.MESA()

m_200_46_newtides_str=mp.MESA()
m3_200_46_newtides_str=mp.MESA()



m_200_25=mp.MESA()
m3_200_25=mp.MESA()

m_200_25_newtides=mp.MESA()
m3_200_25_newtides=mp.MESA()


m_200_25_newtides_str=mp.MESA()
m3_200_25_newtides_str=mp.MESA()


m3_200_46_newtides_str_res04=mp.MESA()


m_200_191.log_fold= os.path.join(paths.data,'HD191495/LOGS1_mesa')
print(m_200_191.log_fold)
m_200_191.loadHistory()


m3_200_191.log_fold=os.path.join(paths.data,'HD191495/LOGS3_mesa')
m3_200_191.loadHistory()


m_200_191_newtides.log_fold=os.path.join(paths.data,'HD191495/LOGS1_notides')
m_200_191_newtides.loadHistory()


m3_200_191_newtides.log_fold=os.path.join(paths.data,'HD191495/LOGS3_notides')
m3_200_191_newtides.loadHistory()


m_200_191_newtides_str.log_fold=os.path.join(paths.data,'HD191495/LOGS1_posydon')
m_200_191_newtides_str.loadHistory()


m3_200_191_newtides_str.log_fold=os.path.join(paths.data,'HD191495/LOGS3_posydon')
m3_200_191_newtides_str.loadHistory()




m_200_46_newtides.log_fold=os.path.join(paths.data,'HD46485/LOGS1_notides')
m_200_46_newtides.loadHistory()


m3_200_46_newtides.log_fold=os.path.join(paths.data,'HD46485/LOGS3_notides')
m3_200_46_newtides.loadHistory()




m_200_46_newtides_str.log_fold=os.path.join(paths.data,'HD46485/LOGS1_posydon')
m_200_46_newtides_str.loadHistory()


m3_200_46_newtides_str.log_fold=os.path.join(paths.data,'HD46485/LOGS3_posydon')
m3_200_46_newtides_str.loadHistory()


m3_200_46_newtides_str_res04.log_fold=os.path.join(paths.data,'HD46485/LOGS1_posydon_res')
m3_200_46_newtides_str_res04.loadHistory()



m_200_46.log_fold=os.path.join(paths.data,'HD46485/LOGS1_mesa')
m_200_46.loadHistory()

m3_200_46.log_fold=os.path.join(paths.data,'HD46485/LOGS3_mesa')
m3_200_46.loadHistory()




m_200_25.log_fold=os.path.join(paths.data,'HD25631/LOGS1_mesa')
m_200_25.loadHistory()


m3_200_25.log_fold=os.path.join(paths.data,'HD25631/LOGS3_mesa')
m3_200_25.loadHistory()

m_200_25_newtides.log_fold=os.path.join(paths.data,'HD25631/LOGS1_notides')
m_200_25_newtides.loadHistory()



m3_200_25_newtides.log_fold=os.path.join(paths.data,'HD25631/LOGS3_notides')
m3_200_25_newtides.loadHistory()


m_200_25_newtides_str.log_fold=os.path.join(paths.data,'HD25631/LOGS1_posydon')
m_200_25_newtides_str.loadHistory()


m3_200_25_newtides_str.log_fold=os.path.join(paths.data,'HD25631/LOGS3_posydon')
m3_200_25_newtides_str.loadHistory()






star_age_200_191=m_200_191.hist.star_age
period_days_191=m3_200_191.hist.period_days

surf_avg_vtor_200_191=m_200_191.hist.surf_avg_v_rot
surf_avg_omega200_191=m_200_191.hist.surf_avg_omega
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










star_age_200_46=m_200_46.hist.star_age
period_days_46=m3_200_46.hist.period_days
surf_avg_vtor_200_46=m_200_46.hist.surf_avg_v_rot
surf_avg_omega200_46=m_200_46.hist.surf_avg_omega
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


star_age_200_25=m_200_25.hist.star_age
period_days_25=m3_200_25.hist.period_days
surf_avg_vtor_200_25=m_200_25.hist.surf_avg_v_rot
surf_avg_omega200_25=m_200_25.hist.surf_avg_omega
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




log_Teff_1_200_46_notides=m_200_46_newtides.hist.log_Teff
LOGL_1_200_46_notides=m_200_46_newtides.hist.log_L


log_Teff_1_200_46_newtides=m_200_46_newtides_str.hist.log_Teff
LOGL_1_200_46_newtides=m_200_46_newtides_str.hist.log_L




star_age_200_46_newtides_res04=m3_200_46_newtides_str_res04.hist.star_age
surf_avg_vtor_200_46_newtides_res04=m3_200_46_newtides_str_res04.hist.surf_avg_v_rot
log_Teff_1_200_46_newtides_res04=m3_200_46_newtides_str_res04.hist.log_Teff
LOGL_1_200_46_newtides_res04=m3_200_46_newtides_str_res04.hist.log_L


star_age_200_46_newtides=m_200_46_newtides.hist.star_age
surf_avg_vtor_200_46_newtides=m_200_46_newtides.hist.surf_avg_v_rot
surf_avg_omega200_46_newtides=m_200_46_newtides.hist.surf_avg_omega
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



star_age_200_46_newtides_str=m_200_46_newtides_str.hist.star_age
surf_avg_vtor_200_46_newtides_str=m_200_46_newtides_str.hist.surf_avg_v_rot
surf_avg_omega200_46_newtides_str=m_200_46_newtides_str.hist.surf_avg_omega
star_1_radius200_46_newtides_str=m3_200_46_newtides_str.hist.star_1_radius
star_1_J_orb_200_46_newtides_str=m3_200_46_newtides_str.hist.J_orb
star_1_J_spin_200_46_newtides_str=m3_200_46_newtides_str.hist.J_spin_1
star_2_J_spin_200_46_newtides_str=m3_200_46_newtides_str.hist.J_spin_2
age_200_46_newtides_str=m3_200_46_newtides_str.hist.age
jtotal_200_46_newtides_str=m3_200_46_newtides_str.hist.J_total
log_total_angular_momentum_200_46_newtides_str=m_200_46_newtides_str.hist.log_total_angular_momentum
surf_avg_j_rot_200_46_newtides_str=m_200_46_newtides_str.hist.surf_avg_j_rot
center_h1_200_46_newtides_str=m_200_46_newtides_str.hist.center_h1
stable_tides_200_46_newtides_str=[3*el for el in np.add(star_1_J_spin_200_46_newtides_str,star_2_J_spin_200_46_newtides_str)]

iRLOF_200_46_newtides_str = np.argmax(center_h1_200_46_newtides_str < 1e-2)







star_age_200_25_newtides=m_200_25_newtides.hist.star_age
surf_avg_vtor_200_25_newtides=m_200_25_newtides.hist.surf_avg_v_rot
surf_avg_omega200_25_newtides=m_200_25_newtides.hist.surf_avg_omega
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




star_age_200_25_newtides_str=m_200_25_newtides_str.hist.star_age
surf_avg_vtor_200_25_newtides_str=m_200_25_newtides_str.hist.surf_avg_v_rot
surf_avg_omega200_25_newtides_str=m_200_25_newtides_str.hist.surf_avg_omega
star_1_radius200_25_newtides_str=m3_200_25_newtides_str.hist.star_1_radius
star_1_J_orb_200_25_newtides_str=m3_200_25_newtides_str.hist.J_orb
star_1_J_spin_200_25_newtides_str=m3_200_25_newtides_str.hist.J_spin_1
star_2_J_spin_200_25_newtides_str=m3_200_25_newtides_str.hist.J_spin_2
age_200_25_newtides_str=m3_200_25_newtides_str.hist.age
jtotal_200_25_newtides_str=m3_200_25_newtides_str.hist.J_total
log_total_angular_momentum_200_25_newtides_str=m_200_25_newtides_str.hist.log_total_angular_momentum
surf_avg_j_rot_200_25_newtides_str=m_200_25_newtides_str.hist.surf_avg_j_rot
center_h1_200_25_newtides_str=m_200_25_newtides_str.hist.center_h1

stable_tides_200_20_newtides_str=[3*el for el in np.add(star_1_J_spin_200_25_newtides_str,star_2_J_spin_200_25_newtides_str)]
iRLOF_200_25_newtides_str = np.argmax(center_h1_200_25_newtides_str < 1e-2)




star_age_200_191_newtides=m_200_191_newtides.hist.star_age
surf_avg_vtor_200_191_newtides=m_200_191_newtides.hist.surf_avg_v_rot
surf_avg_omega200_191_newtides=m_200_191_newtides.hist.surf_avg_omega
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






star_age_200_191_newtides_str=m_200_191_newtides_str.hist.star_age
surf_avg_vtor_200_191_newtides_str=m_200_191_newtides_str.hist.surf_avg_v_rot
surf_avg_omega200_191_newtides_str=m_200_191_newtides_str.hist.surf_avg_omega
star_1_radius200_191_newtides_str=m3_200_191_newtides_str.hist.star_1_radius
star_1_J_orb_200_191_newtides_str=m3_200_191_newtides_str.hist.J_orb
star_1_J_spin_200_191_newtides_str=m3_200_191_newtides_str.hist.J_spin_1
star_2_J_spin_200_191_newtides_str=m3_200_191_newtides_str.hist.J_spin_2
age_200_191_newtides_str=m3_200_191_newtides_str.hist.age
jtotal_200_191_newtides_str=m3_200_191_newtides_str.hist.J_total
log_total_angular_momentum_200_191_newtides_str=m_200_191_newtides_str.hist.log_total_angular_momentum
surf_avg_j_rot_200_191_newtides_str=m_200_191_newtides_str.hist.surf_avg_j_rot
center_h1_200_191_newtides_str=m_200_191_newtides_str.hist.center_h1

iRLOF_200_191_newtides_str = np.argmax(center_h1_200_191_newtides_str < 1e-2)

log_Teff_1_200_191_newtides_str=m_200_191_newtides_str.hist.log_Teff
LOGL_1_200_46_191_newtides_str= m_200_191_newtides_str.hist.log_L


v_sync_191=2*3.1415926*star_1_radius200_191*rsun/(period_days_191*24*3600)
v_sync_46=2*3.1415926*star_1_radius200_46*rsun/(period_days_46*24*3600)
v_sync_25=2*3.1415926*star_1_radius200_25*rsun/(period_days_25*24*3600)




pp1_all_res = PdfPages(paths.figures / 'age_vsurf_comparison_resolution_test.pdf') 

fig=plt.figure(figsize=(10, 10))



#ax1.title.set_text('HD46485, $v^{ini}_{rot}$=350 $km~s^{-1}$, M1$_{ini}$=24, M2$_{ini}$=1, p=6.9d',fontsize=24)
plt.title('HD46485, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=350 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=24, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=6.9 d', fontsize=22)


plt.plot(star_age_200_46_newtides_str/1e6, surf_avg_vtor_200_46_newtides_str, linestyle='-',color='black',label='POSYDON, mesh delta = 1, time delta = 1')


plt.plot(star_age_200_46_newtides_res04/1e6, surf_avg_vtor_200_46_newtides_res04,linestyle='--',label='POSYDON, mesh delta = 0.5, time delta = 0.75',color='red')



plt.ylabel('Surface rotational velocity [$km~s^{-1}$]',fontsize=30)
plt.xlabel('Star age [Myrs]',fontsize=30)
plt.legend(loc=3,fontsize=20)

plt.savefig(pp1_all_res, format='pdf')
pp1_all_res.close()




pp1_all = PdfPages(paths.figures / 'age_vsurf_comparison_panel.pdf') 


fig=plt.figure(figsize=(10, 10))

ax1=fig.add_subplot(311)
#ax1.title.set_text('HD46485, $v^{ini}_{rot}$=350 $km~s^{-1}$, M1$_{ini}$=24, M2$_{ini}$=1, p=6.9d',fontsize=24)
ax1.set_title('HD46485, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=350 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=24, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=6.9 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})

plt.plot(star_age_200_46/1e6, surf_avg_vtor_200_46,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_46_newtides_str/1e6, surf_avg_vtor_200_46_newtides_str, linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_46_newtides/1e6, surf_avg_vtor_200_46_newtides,linestyle='-.',color='r',label='No tides')
plt.plot(age_200_46/1e6, v_sync_46,linestyle='--',color='gray',label='$\it{v}_\mathrm{sync}$')

plt.legend(loc=3,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)



ax2=fig.add_subplot(312)



ax2.set_title('HD191495, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=200 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=15, $\it{M}_\mathrm{2,ini}$=1.5, $\it{P}_\mathrm{ini}$=3.6 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(star_age_200_191/1e6, surf_avg_vtor_200_191,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_191_newtides_str/1e6, surf_avg_vtor_200_191_newtides_str,linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_191_newtides/1e6, surf_avg_vtor_200_191_newtides,linestyle='-.',color='r',label='No tides')
plt.plot(age_200_191/1e6, v_sync_191,linestyle='--',color='gray',label='$\it{v}_\mathrm{sync}$')

plt.legend(loc=3,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('Surface rotational velocity [$km~s^{-1}$]',fontsize=24)





ax3=fig.add_subplot(313)


ax3.set_title('HD25631, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=220 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=7, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=5.2 d', fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(star_age_200_25/1e6, surf_avg_vtor_200_25,linestyle='-',color='black',label='MESA default')

plt.plot(star_age_200_25_newtides_str/1e6, surf_avg_vtor_200_25_newtides_str,linestyle='--',color='b',label='POSYDON')

plt.plot(star_age_200_25_newtides/1e6, surf_avg_vtor_200_25_newtides,linestyle='-.',color='r',label='No tides')

plt.plot(age_200_25/1e6, v_sync_25,linestyle='--',color='gray',label='$\it{v}_\mathrm{sync}$')


plt.legend(loc=3,fontsize=16,frameon=True)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Star age [Myrs]',fontsize=24)


plt.savefig(pp1_all, format='pdf')

pp1_all.close()








pp1_hr = PdfPages(paths.figures / 'HR_primary.pdf') 

fig=plt.figure(figsize=(10, 10))


plt.title('HD46485, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=350 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=24, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=6.9 d', fontsize=22)


#plt.title('HD46485 POSYDON tides', fontsize=17)
#plt.tick_params(labelsize=18)


#fig=plt.figure(figsize=(10, 10))



#plt.plot(log_Teff_1_200_46_newtides[0:iRLOF_200_46_newtides], #LOGL_1_200_46_newtides[0:iRLOF_200_46_newtides],color='k',linestyle='-',label='POSYDON, mesh delta = 1, time delta = 1',lw=2)



plt.plot(10 ** (log_Teff_1_200_46_newtides) / 1000, LOGL_1_200_46_newtides,color='k',linestyle='-',label='POSYDON, mesh delta = 1, time delta = 1',lw=2)


#plt.plot(log_Teff_1_200_46_newtides_res05, LOGL_1_200_46_newtides_res05,linestyle='--',label='POSYDON, mesh coeff = 0.75, time delta  = 0.75',color='green')


plt.plot(10 ** (log_Teff_1_200_46_newtides_res04) / 1000, LOGL_1_200_46_newtides_res04,linestyle='--',label='POSYDON, mesh delta = 0.5, time delta = 0.75',color='red')



#plt.plot(log_Teff_1_newtides[star_age_rlof_ind_pos:-1], LOGL_1_newtides[star_age_rlof_ind_pos:-1], lw=5, c='green', zorder=0)


plt.gca().invert_xaxis()
plt.legend(loc=4,fontsize=20)
plt.ylabel('log$_{10} (L/L_{\odot}$)',fontsize=30)
plt.xlabel('$T_{eff}$ [kK]',fontsize=30)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.xlim([5e6,6.5e6])


plt.savefig(pp1_hr, format='pdf')

pp1_hr.close()








rl_relative_gap_1=m3_200_191.hist.rl_relative_overflow_1
rl_relative_gap_1_posydon=m3_200_191_newtides_str.hist.rl_relative_overflow_1

lg_t_sync_1_class=m3_200_191.hist.lg_t_sync_1
lg_t_sync_1_posydon=m3_200_191_newtides_str.hist.lg_t_sync_1

iRLOF_1 = rl_relative_gap_1 > 0
iRLOF_1_posydon = rl_relative_gap_1_posydon > 0





pp1_lgt = PdfPages(paths.figures / 'LOGS_lg_tsync_HD191495.pdf') 


plt.figure(figsize=(10, 10))


plt.title('HD191495, $\it{v}^\mathrm{1,ini}_\mathrm{rot}$=200 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=15, $\it{M}_\mathrm{2,ini}$=1.5, $\it{P}_\mathrm{ini}$=3.6 d', fontsize=22)
#plt.tick_params(labelsize=18)


plt.plot(star_1_radius200_191, lg_t_sync_1_class,color='black',linestyle='-',label='MESA default',lw=2)
plt.plot(star_1_radius200_191[iRLOF_1], lg_t_sync_1_class[iRLOF_1], lw=5, c='black', zorder=0)

plt.plot(star_1_radius200_191_newtides_str, lg_t_sync_1_posydon,linestyle='--',label='POSYDON',color='b')
plt.plot(star_1_radius200_191_newtides_str[iRLOF_1_posydon], lg_t_sync_1_posydon[iRLOF_1_posydon], lw=5, c='b', zorder=0)



plt.ylabel('log$_{10}$ ($\it{t}_\mathrm{1,sync}$)',fontsize=30)
plt.xlabel('$\it{R}_\mathrm{1}$ [$\it{R}_\mathrm{\odot}$]',fontsize=30)
plt.legend(loc=3,fontsize=20,frameon=True)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.xlim([5e6,6.5e6])

plt.savefig(pp1_lgt, format='pdf')
pp1_lgt.close()














