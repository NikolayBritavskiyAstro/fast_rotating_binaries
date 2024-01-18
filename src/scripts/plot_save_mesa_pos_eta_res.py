import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np 
from showyourwork.paths import user as Paths
paths = Paths()
import os
import mesaPlot as mp

plt.style.use(paths.scripts / "matplotlibrc")

if os.path.exists(os.path.join(paths.data,'mass_transfer_efficiency/p3_pos/LOGS3/history.data')):
	pass
else:
	os.system(f'python {os.path.join(paths.scripts / "unzip_MESA_output.py")}')



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



m3_p3_g1_res=mp.MESA()
m3_p5_g1_res=mp.MESA()
m3_p10_g1_res=mp.MESA()
m3_p50_g1_res=mp.MESA()
m3_p70_g1_res=mp.MESA()



m2_p3_g1_res=mp.MESA()
m2_p5_g1_res=mp.MESA()
m2_p10_g1_res=mp.MESA()
m2_p50_g1_res=mp.MESA()
m2_p70_g1_res=mp.MESA()

#name=sys.argv[1]
#print(name+'/LOGS1')







m3_p3_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p3_pos/LOGS3')
m3_p3_g1_new.loadHistory()

m3_p5_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p5_pos/LOGS3')
m3_p5_g1_new.loadHistory()

m3_p10_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p10_pos/LOGS3')
m3_p10_g1_new.loadHistory()


m3_p50_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p50_pos/LOGS3')
m3_p50_g1_new.loadHistory()


m3_p70_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p70_pos/LOGS3')
m3_p70_g1_new.loadHistory()





m2_p3_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p3_pos/LOGS2')
m2_p3_g1_new.loadHistory()

m2_p5_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p5_pos/LOGS2')
m2_p5_g1_new.loadHistory()

m2_p10_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p10_pos/LOGS2')
m2_p10_g1_new.loadHistory()


m2_p50_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p50_pos/LOGS2')
m2_p50_g1_new.loadHistory()


m2_p70_g1_new.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p70_pos/LOGS2')
m2_p70_g1_new.loadHistory()



m3_p3_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p3_res/LOGS3')
m3_p3_g1_res.loadHistory()

m3_p5_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p5_res/LOGS3')
m3_p5_g1_res.loadHistory()

m3_p10_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p10_res/LOGS3')
m3_p10_g1_res.loadHistory()


m3_p50_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p50_res/LOGS3')
m3_p50_g1_res.loadHistory()


m3_p70_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p70_res/LOGS3')
m3_p70_g1_res.loadHistory()





m2_p3_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p3_res/LOGS2')
m2_p3_g1_res.loadHistory()

m2_p5_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p5_res/LOGS2')
m2_p5_g1_res.loadHistory()

m2_p10_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p10_res/LOGS2')
m2_p10_g1_res.loadHistory()


m2_p50_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p50_res/LOGS2')
m2_p50_g1_res.loadHistory()


m2_p70_g1_res.log_fold=os.path.join(paths.data,'mass_transfer_efficiency/p70_res/LOGS2')
m2_p70_g1_res.loadHistory()



indx_p3=[]
indx_p5=[]
indx_p10=[]
indx_p50=[]
indx_p70=[]


indx_p3_res=[]
indx_p5_res=[]
indx_p10_res=[]
indx_p50_res=[]
indx_p70_res=[]










age_p3_g1_res=m3_p3_g1_res.hist.age
delta_p3_g1_res=m3_p3_g1_res.hist.delta_updated
spin2_p3_g1_res=m3_p3_g1_res.hist.J_spin_2


period_days_p3_res=m3_p3_g1_res.hist.period_days
omega_sync_p3=2*3.1415926/(period_days_p3_res) #rad/days



star_age_p3_g1_res=m2_p3_g1_res.hist.star_age
surf_avg_omega_div_omega_crit_p3_g1_res=m2_p3_g1_res.hist.surf_avg_omega_div_omega_crit
surf_avg_omega_crit_p3_g1_res=m2_p3_g1_res.hist.surf_avg_omega_crit

age_p5_g1_res=m3_p5_g1_res.hist.age
delta_p5_g1_res=m3_p5_g1_res.hist.delta_updated
spin2_p5_g1_res=m3_p5_g1_res.hist.J_spin_2

period_days_p5_res=m3_p5_g1_res.hist.period_days
omega_sync_p5_res=2*3.1415926/(period_days_p5_res) #rad/days


star_age_p5_g1_res=m2_p5_g1_res.hist.star_age
surf_avg_omega_div_omega_crit_p5_g1_res=m2_p5_g1_res.hist.surf_avg_omega_div_omega_crit
surf_avg_omega_crit_p5_g1_res=m2_p5_g1_res.hist.surf_avg_omega_crit



age_p10_g1_res=m3_p10_g1_res.hist.age
delta_p10_g1_res=m3_p10_g1_res.hist.delta_updated
spin2_p10_g1_res=m3_p10_g1_res.hist.J_spin_2


star_age_p10_g1_res=m2_p10_g1_res.hist.star_age
surf_avg_omega_div_omega_crit_p10_g1_res=m2_p10_g1_res.hist.surf_avg_omega_div_omega_crit


age_p50_g1_res=m3_p50_g1_res.hist.age
delta_p50_g1_res=m3_p50_g1_res.hist.delta_updated
spin2_p50_g1_res=m3_p50_g1_res.hist.J_spin_2

star_age_p50_g1_res=m2_p50_g1_res.hist.star_age
surf_avg_omega_div_omega_crit_p50_g1_res=m2_p50_g1_res.hist.surf_avg_omega_div_omega_crit


age_p70_g1_res=m3_p70_g1_res.hist.age
delta_p70_g1_res=m3_p70_g1_res.hist.delta_updated
spin2_p70_g1_res=m3_p70_g1_res.hist.J_spin_2
period_days_p70_res=m3_p70_g1_res.hist.period_days



star_age_p70_g1_res=m2_p70_g1_res.hist.star_age
surf_avg_omega_div_omega_crit_p70_g1_res=m2_p70_g1_res.hist.surf_avg_omega_div_omega_crit



for i in range(len(age_p3_g1_res)):
    indx_p3_res.append(find_nearest(star_age_p3_g1_res, age_p3_g1_res[i]))


for i in range(len(age_p5_g1_res)):
    indx_p5_res.append(find_nearest(star_age_p5_g1_res, age_p5_g1_res[i]))


for i in range(len(age_p10_g1_res)):
    indx_p10_res.append(find_nearest(star_age_p10_g1_res, age_p10_g1_res[i]))

for i in range(len(age_p50_g1_res)):
    indx_p50_res.append(find_nearest(star_age_p50_g1_res, age_p50_g1_res[i]))


for i in range(len(age_p70_g1_res)):
    indx_p70_res.append(find_nearest(star_age_p70_g1_res, age_p70_g1_res[i]))






















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







pp1_period = PdfPages(paths.figures / 'eta_omega_res.pdf') 



plt.figure(figsize=(10, 10))


plt.title('$\it{M}_\mathrm{don,ini}$ = 30, $\it{M}_\mathrm{acc,ini}$ = 20, $\it{\gamma}$ = 3',fontsize=30)
#plt.tick_params(labelsize=18)


param=24*3600




plt.plot(surf_avg_omega_div_omega_crit_p3_g1_res[indx_p3_res][:-100], 1-delta_p3_g1_res[:-100],label='$\it{P}_\mathrm{ini}$ = 3 d',lw=3)

plt.plot(surf_avg_omega_div_omega_crit_p5_g1_res[indx_p5_res][:-100], 1-delta_p5_g1_res[:-100],linestyle='-',label='$\it{P}_\mathrm{ini}$ = 5 d',lw=3)


plt.plot(surf_avg_omega_div_omega_crit_p10_g1_res[indx_p10_res], 1-delta_p10_g1_res,linestyle='-',label='$\it{P}_\mathrm{ini}$ = 10 d',lw=3)
plt.plot(surf_avg_omega_div_omega_crit_p50_g1_res[indx_p50_res], 1-delta_p50_g1_res,linestyle='-',label='$\it{P}_\mathrm{ini}$ = 50 d',lw=3)
plt.plot(surf_avg_omega_div_omega_crit_p70_g1_res[indx_p70_res], 1-delta_p70_g1_res,linestyle='-',label='$\it{P}_\mathrm{ini}$ = 70 d',lw=3)



plt.text(0.01,0.23,'mesh delta = 0.5, time delta = 0.75',fontsize=18)









plt.plot(surf_avg_omega_div_omega_crit_p3_g1_new[indx_p3][:-150], 1-delta_p3_g1_new[:-150],lw=2,color='gray')

plt.plot(surf_avg_omega_div_omega_crit_p5_g1_new[indx_p5][:-150], 1-delta_p5_g1_new[:-150],linestyle='-',color='gray',lw=2)


plt.plot(surf_avg_omega_div_omega_crit_p10_g1_new[indx_p10], 1-delta_p10_g1_new,linestyle='-',color='gray',lw=2)
plt.plot(surf_avg_omega_div_omega_crit_p50_g1_new[indx_p50][:-30], 1-delta_p50_g1_new[:-30],linestyle='-',color='gray',lw=2)
plt.plot(surf_avg_omega_div_omega_crit_p70_g1_new[indx_p70][:-30], 1-delta_p70_g1_new[:-30],color='gray',lw=2)





plt.xlabel('$\it{\omega}/\it{\omega}_\mathrm{crit}$',fontsize=30)
plt.ylabel('$\it{\eta} = \Delta \it{M}_\mathrm{acc}/\Delta \it{M}_\mathrm{don}$',fontsize=30)

#plt.xlabel('Star age [Myrs]',fontsize=30)
plt.legend(loc=1,fontsize=23,frameon=True)
#plt.xlim([5e6,6.5e6])

plt.ylim([-0.02,0.25])

plt.savefig(pp1_period, format='pdf')


pp1_period.close()






