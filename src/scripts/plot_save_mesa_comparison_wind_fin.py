import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np 
from showyourwork.paths import user as Paths
paths = Paths()
import os
import mesaPlot as mp

plt.style.use(paths.scripts / "matplotlibrc")


colors_list = plt.cm.viridis(np.linspace(0,1, 8))



m_1=mp.MESA()
m_50=mp.MESA()
m_100=mp.MESA()
m_200=mp.MESA()
m_300=mp.MESA()
m_350=mp.MESA()


m_50_b=mp.MESA()
m_100_b=mp.MESA()
m_200_b=mp.MESA()
m_300_b=mp.MESA()
m_340_b=mp.MESA()


m_50_st=mp.MESA()
m_50_f=mp.MESA()
m_350_f=mp.MESA()
m_200_f=mp.MESA()
m_300_f=mp.MESA()
m_100_f=mp.MESA()

m_200_st=mp.MESA()
m_100_st=mp.MESA()
m_350_st=mp.MESA()
m_300_st=mp.MESA()



m_50_st.log_fold=os.path.join(paths.data,'one_star/M24_ST/LOGS50')
m_50_st.loadHistory()
m_200_st.log_fold=os.path.join(paths.data,'one_star/M24_ST/LOGS200')
m_200_st.loadHistory()
m_100_st.log_fold=os.path.join(paths.data,'one_star/M24_ST/LOGS100')
m_100_st.loadHistory()
m_300_st.log_fold=os.path.join(paths.data,'one_star/M24_ST/LOGS300')
m_300_st.loadHistory()
m_350_st.log_fold=os.path.join(paths.data,'one_star/M24_ST/LOGS350')
m_350_st.loadHistory()
m_50_f.log_fold=os.path.join(paths.data,'one_star/M24_F/LOGS_50')
m_50_f.loadHistory()

m_350_f.log_fold=os.path.join(paths.data,'one_star/M24_F/LOGS_350')
m_350_f.loadHistory()

m_200_f.log_fold=os.path.join(paths.data,'one_star/M24_F/LOGS_200')
m_200_f.loadHistory()

m_300_f.log_fold=os.path.join(paths.data,'one_star/M24_F/LOGS_300')
m_300_f.loadHistory()

m_100_f.log_fold=os.path.join(paths.data,'one_star/M24_F/LOGS_100')
m_100_f.loadHistory()


m_50_b.log_fold=os.path.join(paths.data,'one_star/M24_Bjorklund_wind/LOGS50_new_b')
m_50_b.loadHistory()

m_100_b.log_fold=os.path.join(paths.data,'one_star/M24_Bjorklund_wind/LOGS100_new_b')
m_100_b.loadHistory()

m_200_b.log_fold=os.path.join(paths.data,'one_star/M24_Bjorklund_wind/LOGS200_new_b')
m_200_b.loadHistory()

m_300_b.log_fold=os.path.join(paths.data,'one_star/M24_Bjorklund_wind/LOGS300_new_b')
m_300_b.loadHistory()

m_340_b.log_fold=os.path.join(paths.data,'one_star/M24_Bjorklund_wind/LOGS340_new_b')
m_340_b.loadHistory()


m_1.log_fold=os.path.join(paths.data,'one_star/M24_Vink_wind/LOGS200_res05')
m_1.loadHistory()

m_50.log_fold=os.path.join(paths.data,'one_star/M24_Vink_wind/LOGS50')
m_50.loadHistory()

m_100.log_fold=os.path.join(paths.data,'one_star/M24_Vink_wind/LOGS100')
m_100.loadHistory()

m_200.log_fold=os.path.join(paths.data,'one_star/M24_Vink_wind/LOGS200')
m_200.loadHistory()

m_300.log_fold=os.path.join(paths.data,'one_star/M24_Vink_wind/LOGS300')
m_300.loadHistory()

m_350.log_fold=os.path.join(paths.data,'one_star/M24_Vink_wind/LOGS350')
m_350.loadHistory()


star_age_1=m_1.hist.star_age
surf_avg_vtor_1=m_1.hist.surf_avg_v_rot
center_h1_1=m_1.hist.center_h1
surface_n14_1=(1.00797/14.0067)*m_1.hist.surface_n14/m_1.hist.surface_h1
center_n14_1=m_1.hist.center_n14

surface_n14_50=(1.00797/14.0067)*m_50.hist.surface_n14/m_50.hist.surface_h1
center_n14_50=m_50.hist.center_n14

surface_n14_100=(1.00797/14.0067)*m_100.hist.surface_n14/m_100.hist.surface_h1
center_n14_100=m_100.hist.center_n14


surface_n14_200=(1.00797/14.0067)*m_200.hist.surface_n14/m_200.hist.surface_h1
center_n14_200=m_200.hist.center_n14


surface_n14_300=(1.00797/14.0067)*m_300.hist.surface_n14/m_300.hist.surface_h1
center_n14_300=m_300.hist.center_n14


surface_n14_350=(1.00797/14.0067)*m_350.hist.surface_n14/m_350.hist.surface_h1
center_n14_350=m_350.hist.center_n14


surface_n14_50b=(1.00797/14.0067)*m_50_b.hist.surface_n14/m_50_b.hist.surface_h1
surface_n14_100b=(1.00797/14.0067)*m_100_b.hist.surface_n14/m_100_b.hist.surface_h1
surface_n14_200b=(1.00797/14.0067)*m_200_b.hist.surface_n14/m_200_b.hist.surface_h1
surface_n14_300b=(1.00797/14.0067)*m_300_b.hist.surface_n14/m_300_b.hist.surface_h1
surface_n14_340b=(1.00797/14.0067)*m_340_b.hist.surface_n14/m_340_b.hist.surface_h1




surface_n14_50f=(1.00797/14.0067)*m_50_f.hist.surface_n14/m_50_f.hist.surface_h1
surface_n14_100f=(1.00797/14.0067)*m_100_f.hist.surface_n14/m_100_f.hist.surface_h1
surface_n14_200f=(1.00797/14.0067)*m_200_f.hist.surface_n14/m_200_f.hist.surface_h1
surface_n14_300f=(1.00797/14.0067)*m_300_f.hist.surface_n14/m_300_f.hist.surface_h1
surface_n14_350f=(1.00797/14.0067)*m_350_f.hist.surface_n14/m_350_f.hist.surface_h1



surface_c12_1=(1.00797/12.011)*m_1.hist.surface_c12/m_1.hist.surface_h1
surface_c12_50=(1.00797/12.011)*m_50.hist.surface_c12/m_50.hist.surface_h1
surface_c12_100=(1.00797/12.011)*m_100.hist.surface_c12/m_100.hist.surface_h1
surface_c12_200=(1.00797/12.011)*m_200.hist.surface_c12/m_200.hist.surface_h1
surface_c12_300=(1.00797/12.011)*m_300.hist.surface_c12/m_300.hist.surface_h1
surface_c12_350=(1.00797/12.011)*m_350.hist.surface_c12/m_350.hist.surface_h1

surface_o16_1=(1.00797/15.9994)*m_1.hist.surface_o16/m_1.hist.surface_h1
surface_o16_50=(1.00797/15.9994)*m_50.hist.surface_o16/m_50.hist.surface_h1
surface_o16_100=(1.00797/15.9994)*m_100.hist.surface_o16/m_100.hist.surface_h1
surface_o16_200=(1.00797/15.9994)*m_200.hist.surface_o16/m_200.hist.surface_h1
surface_o16_300=(1.00797/15.9994)*m_300.hist.surface_o16/m_300.hist.surface_h1
surface_o16_350=(1.00797/15.9994)*m_350.hist.surface_o16/m_350.hist.surface_h1



LOGL_1=m_1.hist.log_L
log_Teff_1=m_1.hist.log_Teff
iRLOF_1 = np.argmax(center_h1_1 < 1e-2)

star_mdot_300=m_300.hist.star_mdot
radius_300=m_300.hist.surf_avg_omega
star_mdot_300=m_300.hist.star_mdot
jtot_300=m_300.hist.log_total_angular_momentum



star_age_50=m_50.hist.star_age
surf_avg_vtor_50=m_50.hist.surf_avg_v_rot
center_h1_50=m_50.hist.center_h1
iRLOF_50 = np.argmax(center_h1_50 < 1e-2)


star_age_50_b=m_50_b.hist.star_age
surf_avg_vtor_50_b=m_50_b.hist.surf_avg_v_rot
center_h1_50_b=m_50_b.hist.center_h1
iRLOF_50_b = np.argmax(center_h1_50_b < 1e-2)


star_age_50_st=m_50_st.hist.star_age
surf_avg_vtor_50_st=m_50_st.hist.surf_avg_v_rot
center_h1_50_st=m_50_st.hist.center_h1
iRLOF_50_st = np.argmax(center_h1_50_st < 1e-2)


star_age_350_st=m_350_st.hist.star_age
surf_avg_vtor_350_st=m_350_st.hist.surf_avg_v_rot
center_h1_350_st=m_350_st.hist.center_h1
iRLOF_350_st = np.argmax(center_h1_350_st < 1e-2)

star_age_100_st=m_100_st.hist.star_age
surf_avg_vtor_100_st=m_100_st.hist.surf_avg_v_rot
center_h1_100_st=m_100_st.hist.center_h1
iRLOF_100_st = np.argmax(center_h1_100_st < 1e-2)

star_age_300_st=m_300_st.hist.star_age
surf_avg_vtor_300_st=m_300_st.hist.surf_avg_v_rot
center_h1_300_st=m_300_st.hist.center_h1
iRLOF_300_st = np.argmax(center_h1_300_st < 1e-2)
radius_300_st=m_300_st.hist.surf_avg_omega



star_age_50_f=m_50_f.hist.star_age
surf_avg_vtor_50_f=m_50_f.hist.surf_avg_v_rot
center_h1_50_f=m_50_f.hist.center_h1
iRLOF_50_f = np.argmax(center_h1_50_f < 1e-2)


star_age_350_f=m_350_f.hist.star_age
surf_avg_vtor_350_f=m_350_f.hist.surf_avg_v_rot
center_h1_350_f=m_350_f.hist.center_h1
iRLOF_350_f = np.argmax(center_h1_350_f < 1e-2)



star_age_200_f=m_200_f.hist.star_age
surf_avg_vtor_200_f=m_200_f.hist.surf_avg_v_rot
center_h1_200_f=m_200_f.hist.center_h1
iRLOF_200_f = np.argmax(center_h1_200_f < 1e-2)

star_age_100_f=m_100_f.hist.star_age
surf_avg_vtor_100_f=m_100_f.hist.surf_avg_v_rot
center_h1_100_f=m_100_f.hist.center_h1
iRLOF_100_f = np.argmax(center_h1_100_f < 1e-2)

star_age_300_f=m_300_f.hist.star_age
surf_avg_vtor_300_f=m_300_f.hist.surf_avg_v_rot
center_h1_300_f=m_300_f.hist.center_h1
iRLOF_300_f = np.argmax(center_h1_300_f < 1e-2)

star_mdot_300_f=m_300_f.hist.star_mdot
radius_300_f=m_300_f.hist.surf_avg_omega
jtot_300_f=m_300_f.hist.log_total_angular_momentum



star_age_100_b=m_100_b.hist.star_age
surf_avg_vtor_100_b=m_100_b.hist.surf_avg_v_rot
center_h1_100_b=m_100_b.hist.center_h1
iRLOF_100_b = np.argmax(center_h1_100_b < 1e-2)


star_age_200_b=m_200_b.hist.star_age
surf_avg_vtor_200_b=m_200_b.hist.surf_avg_v_rot
center_h1_200_b=m_200_b.hist.center_h1
iRLOF_200_b = np.argmax(center_h1_200_b < 1e-2)


star_age_200_st=m_200_st.hist.star_age
surf_avg_vtor_200_st=m_200_st.hist.surf_avg_v_rot
center_h1_200_st=m_200_st.hist.center_h1
iRLOF_200_st = np.argmax(center_h1_200_st < 1e-2)


star_age_300_b=m_300_b.hist.star_age
surf_avg_vtor_300_b=m_300_b.hist.surf_avg_v_rot
center_h1_300_b=m_300_b.hist.center_h1
iRLOF_300_b = np.argmax(center_h1_300_b < 1e-2)

star_mdot_300_b=m_300_b.hist.star_mdot
radius_300_b=m_300_b.hist.surf_avg_omega
jtot_300_b=m_300_b.hist.log_total_angular_momentum


star_age_340_b=m_340_b.hist.star_age
surf_avg_vtor_340_b=m_340_b.hist.surf_avg_v_rot
center_h1_340_b=m_340_b.hist.center_h1
iRLOF_340_b = np.argmax(center_h1_340_b < 1e-2)



#star_mdot
star_age_100=m_100.hist.star_age
surf_avg_vtor_100=m_100.hist.surf_avg_v_rot
center_h1_100=m_100.hist.center_h1
iRLOF_100 = np.argmax(center_h1_100 < 1e-2)



star_age_200=m_200.hist.star_age
surf_avg_vtor_200=m_200.hist.surf_avg_v_rot
center_h1_200=m_200.hist.center_h1
iRLOF_200 = np.argmax(center_h1_200 < 1e-2)
LOGL_200=m_200.hist.log_L
log_Teff_200=m_200.hist.log_Teff
star_age_200=m_200.hist.star_age


star_age_300=m_300.hist.star_age
surf_avg_vtor_300=m_300.hist.surf_avg_v_rot
center_h1_300=m_300.hist.center_h1
iRLOF_300 = np.argmax(center_h1_300 < 1e-2)


star_age_350=m_350.hist.star_age
surf_avg_vtor_350=m_350.hist.surf_avg_v_rot
center_h1_350=m_350.hist.center_h1
iRLOF_350 = np.argmax(center_h1_350 < 1e-2)



pp1_all = PdfPages(paths.figures / 'age_vsurf_comparison_onestar.pdf') 



plt.figure(figsize=(10, 10))

#plt.title("M$_{ini}$ = 24 M$_{\odot}$", fontsize=30)




#plt.plot(star_age_1[0: iRLOF_1]/1e6, surf_avg_vtor_1[0: iRLOF_1],label='$v^{ini}_{rot}$=1 $km~s^{-1}$',c= colors_list[0])

plt.text(10,100,"$\it{M}_\mathrm{ini}$ = 24 $\it{M}_\mathrm{\odot}$", fontsize=30)


#plt.plot(star_age_50/1e6, surf_avg_vtor_50,linestyle='-',color='magenta',label='$v^{ini}_{rot}$=50 $km~s^{-1}$',lw=3)
#plt.plot(star_age_50[iRLOF_50]/1e6, surf_avg_vtor_50[iRLOF_50],linestyle='-',color='magenta',lw=5)

plt.plot(star_age_50[0:iRLOF_50]/1e6, surf_avg_vtor_50[0:iRLOF_50],label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=50 $km~s^{-1}$',c= colors_list[1],lw=3)


#plt.plot(star_age_100/1e6, surf_avg_vtor_100,linestyle='-',color='gray',label='$v^{ini}_{rot}$=100 $km~s^{-1}$',lw=2)
#plt.plot(star_age_100[iRLOF_100]/1e6, surf_avg_vtor_100[iRLOF_100],linestyle='-',color='gray',lw=5)

plt.plot(star_age_100[0:iRLOF_100]/1e6, surf_avg_vtor_100[0:iRLOF_100],label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=100 $km~s^{-1}$',c= colors_list[2],lw=3)
plt.plot(star_age_200[0:iRLOF_200]/1e6, surf_avg_vtor_200[0:iRLOF_200],label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$',c= colors_list[3],lw=3)

#plt.plot(star_age_150/1e6, surf_avg_vtor_150,linestyle='-',color='orange',label='$v^{ini}_{rot}$=150 $km~s^{-1}$',lw=2)
#plt.plot(star_age_150[iRLOF_150]/1e6, surf_avg_vtor_150[iRLOF_150],linestyle='-',color='orange',lw=5)



#plt.plot(star_age_150[0:iRLOF_150]/1e6, surf_avg_vtor_150[0:iRLOF_150],label='$v^{ini}_{rot}$=150 $km~s^{-1}$',c= colors_list[3])



#plt.plot(star_age_200/1e6, surf_avg_vtor_200,linestyle='-',color='green',label='$v^{ini}_{rot}$=200 $km~s^{-1}$',lw=2)

#plt.plot(star_age_200[iRLOF_200]/1e6, surf_avg_vtor_200[iRLOF_200],linestyle='-',color='green',lw=5)


#plt.plot(star_age_200[0:iRLOF_200]/1e6, surf_avg_vtor_200[0:iRLOF_200],label='$v^{ini}_{rot}$=200 $km~s^{-1} old$',c= colors_list[4])



#plt.plot(star_age_250/1e6, surf_avg_vtor_250,linestyle='-',color='blue',label='$v^{ini}_{rot}$=250 $km~s^{-1}$',lw=2)
#plt.plot(star_age_250[iRLOF_250]/1e6, surf_avg_vtor_250[iRLOF_250],linestyle='-',color='blue',lw=5)


#plt.plot(star_age_250[0:iRLOF_250]/1e6, surf_avg_vtor_250[0:iRLOF_250],label='$v^{ini}_{rot}$=250 $km~s^{-1}$',c= colors_list[5])


#plt.plot(star_age_250[0:iRLOF_250]/1e6, surf_avg_vtor_250[0:iRLOF_250],label='$v^{ini}_{rot}$=100 $km~s^{-1}$ res=1')

#plt.plot(star_age_300/1e6, surf_avg_vtor_300,linestyle='-',color='red',label='$v^{ini}_{rot}$=100 $km~s^{-1}$ res=0.3 time=0.5',lw=2)
plt.plot(star_age_300[0:iRLOF_300]/1e6, surf_avg_vtor_300[0:iRLOF_300],linestyle='-',c= colors_list[4],label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=300 $km~s^{-1}$',lw=3)




plt.plot(star_age_350[0:iRLOF_350]/1e6, surf_avg_vtor_350[0:iRLOF_350],linestyle='-',label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=350 $km~s^{-1}$',lw=3,c= colors_list[5])







plt.plot(star_age_50_st/1e6, surf_avg_vtor_50_st, c= colors_list[1],linestyle='-.',lw=2,label='AM (no S-T dynamo)')
plt.plot(star_age_100_st[0:iRLOF_100_st]/1e6, surf_avg_vtor_100_st[0:iRLOF_100_st],c= colors_list[2],linestyle='-.',lw=2)
plt.plot(star_age_200_st[0:iRLOF_200_st]/1e6, surf_avg_vtor_200_st[0:iRLOF_200_st],c= colors_list[3],linestyle='-.',lw=2)
plt.plot(star_age_300_st[0:iRLOF_300_st]/1e6, surf_avg_vtor_300_st[0:iRLOF_300_st],c= colors_list[4],linestyle='-.',lw=2)
plt.plot(star_age_350_st[0:iRLOF_350_st]/1e6, surf_avg_vtor_350_st[0:iRLOF_350_st],c= colors_list[5],linestyle='-.',lw=2)

#plt.plot(star_age_350_b/1e6, surf_avg_vtor_350_b,c='gray',linestyle='--',lw=1)

plt.plot(star_age_50_f/1e6, surf_avg_vtor_50_f,c= colors_list[1], linestyle=':',lw=2,label='AM (Fuller+19)')
plt.plot(star_age_100_f[0:iRLOF_100_f]/1e6, surf_avg_vtor_100_f[0:iRLOF_100_f],c= colors_list[2], linestyle=':',lw=2)
plt.plot(star_age_200_f[0:iRLOF_200_f]/1e6, surf_avg_vtor_200_f[0:iRLOF_200_f],c= colors_list[3], linestyle=':',lw=2)
plt.plot(star_age_300_f[0:iRLOF_300_f]/1e6, surf_avg_vtor_300_f[0:iRLOF_300_f],c= colors_list[4], linestyle=':',lw=2)
plt.plot(star_age_350_f[0:iRLOF_350_f]/1e6, surf_avg_vtor_350_f[0:iRLOF_350_f],c= colors_list[5], linestyle=':',lw=2)





plt.plot(star_age_50_b[0:iRLOF_50_b]/1e6, surf_avg_vtor_50_b[0:iRLOF_50_b],label='Wind (Bjorklund+23)',c= colors_list[1],linestyle='--',lw=2)
plt.plot(star_age_100_b[0:iRLOF_100_b]/1e6, surf_avg_vtor_100_b[0:iRLOF_100_b],c= colors_list[2],linestyle='--',lw=2)
plt.plot(star_age_200_b[0:iRLOF_200_b]/1e6, surf_avg_vtor_200_b[0:iRLOF_200_b],c= colors_list[3],linestyle='--',lw=2)
plt.plot(star_age_300_b[0:iRLOF_300_b]/1e6, surf_avg_vtor_300_b[0:iRLOF_300_b],c= colors_list[4],linestyle='--',lw=2)
#plt.plot(star_age_350_b[0:iRLOF_350_b]/1e6, surf_avg_vtor_350_b[0:iRLOF_350_b],c='gray',linestyle='--',lw=1)
plt.plot(star_age_340_b[0:iRLOF_340_b]/1e6, surf_avg_vtor_340_b[0:iRLOF_340_b],c= colors_list[5],linestyle='--',lw=2)
#plt.plot(star_age_330_b[0:iRLOF_330_b]/1e6, surf_avg_vtor_330_b[0:iRLOF_330_b],c='gray',linestyle='--',lw=1)
#plt.plot(star_age_320_b[0:iRLOF_320_b]/1e6, surf_avg_vtor_320_b[0:iRLOF_320_b],c='gray',linestyle='--',lw=1)
#plt.plot(star_age_310_b[0:iRLOF_310_b]/1e6, surf_avg_vtor_310_b[0:iRLOF_310_b],c='gray',linestyle='--',lw=1)




plt.errorbar(2,330,20,1.1, mfc='red',ecolor='gray',linestyle='None')


plt.ylabel('Surface rotational velocity [$km~s^{-1}$]',fontsize=30)
plt.xlabel('Star age [Myrs]',fontsize=30)
plt.legend(loc=1,fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim([-1,15.5])
plt.savefig(pp1_all, format='pdf')

pp1_all.close()

pp1_all_res = PdfPages(paths.figures / 'age_vsurf_comparison_onestar_resolution.pdf') 

plt.figure(figsize=(10, 10))
plt.title("$\it{M}_\mathrm{ini}$ = 24 $\it{M}_\mathrm{\odot}$", fontsize=30)
plt.plot(10**log_Teff_200/1000, LOGL_200,linestyle='-',label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=1, time delta=1',lw=1,c= colors_list[3])

plt.plot(10**log_Teff_200[0:iRLOF_200]/1000, LOGL_200[0:iRLOF_200],linestyle='-',c= colors_list[3],lw=5)


plt.plot(10**log_Teff_1/1000, LOGL_1,linestyle='--',color='red',lw=1)

plt.plot(10**log_Teff_1[0:iRLOF_1]/1000, LOGL_1[0:iRLOF_1],linestyle='--',color='red',lw=8,label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=0.5, time delta=0.75')

#plt.plot(10**log_Teff_350_b, LOGL_350_b,linestyle='-',color='green',lw=2,label='$v^{ini}_{rot}$=350 $km~s^{-1}$, B+23 mesh delta=0.5, time delta=0.6')


#mesh_delta_coeff = 0.5d0
#time_delta_coeff = 0.75d0

plt.gca().invert_xaxis()

plt.legend(loc=2,fontsize=16)
plt.ylabel('log$_{10} (L/L_{\odot}$)',fontsize=30)
plt.xlabel('$T_{eff}$ [kK]',fontsize=30)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)

#plt.plot(star_age_200, center_h1_200,linestyle='-',color='green',label='$v^{ini}_{rot}$=200 $km~s^{-1}$',lw=1)
#plt.plot(star_age_250, center_h1_250,linestyle='-',color='blue',label='$v^{ini}_{rot}$=250 $km~s^{-1}$',lw=2)

plt.savefig(pp1_all_res, format='pdf')

pp1_all_res.close()

pp1_all_res_vrot= PdfPages(paths.figures / 'age_vsurf_comparison_onestar_resolution_vrot.pdf') 


plt.figure(figsize=(10, 10))
plt.title("$\it{M}_\mathrm{ini}$ = 24 $\it{M}_\mathrm{\odot}$", fontsize=30)
plt.plot(star_age_1[0:iRLOF_1]/1e6, surf_avg_vtor_1[0:iRLOF_1],color='red',label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=0.5, time delta=0.75',lw=8,linestyle='--')
#plt.plot(star_age_1[iRLOF_1]/1e6, surf_avg_vtor_1[iRLOF_1],color='red',lw=5,linestyle='-')

plt.plot(star_age_200[0:iRLOF_200]/1e6, surf_avg_vtor_200[0:iRLOF_200],label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=1, time delta=1',c= colors_list[3])

plt.ylabel('Surface rotational velocity [$km~s^{-1}$]',fontsize=30)
plt.xlabel('Star age [Myrs]',fontsize=30)
plt.legend(loc=1,fontsize=16)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
plt.xlim([-1,10.5])
plt.savefig(pp1_all_res_vrot, format='pdf')
pp1_all_res_vrot.close()




'''



pp1_all_res = PdfPages(paths.figures / 'age_vsurf_comparison_onestar_resolution.pdf') 

plt.figure(figsize=(10, 10))
plt.title("$\it{M}_\mathrm{ini}$ = 24 $\it{M}_\mathrm{\odot}$", fontsize=30)
plt.plot(10**log_Teff_200/1000, LOGL_200,linestyle='-',label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=1, time delta=1',lw=1,c= colors_list[3])

plt.plot(10**log_Teff_200[0:iRLOF_200]/1000, LOGL_200[0:iRLOF_200],linestyle='-',c= colors_list[3],lw=5)


plt.plot(10**log_Teff_1/1000, LOGL_1,linestyle='--',color='red',lw=1)

plt.plot(10**log_Teff_1[0:iRLOF_1]/1000, LOGL_1[0:iRLOF_1],linestyle='--',color='red',lw=8,label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=0.5, time delta=0.75')

#plt.plot(10**log_Teff_350_b, LOGL_350_b,linestyle='-',color='green',lw=2,label='$v^{ini}_{rot}$=350 $km~s^{-1}$, B+23 mesh delta=0.5, time delta=0.6')


#mesh_delta_coeff = 0.5d0
#time_delta_coeff = 0.75d0

plt.gca().invert_xaxis()

plt.legend(loc=2,fontsize=16)
plt.ylabel('log$_{10} (L/L_{\odot}$)',fontsize=30)
plt.xlabel('$T_{eff}$ [kK]',fontsize=30)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)

#plt.plot(star_age_200, center_h1_200,linestyle='-',color='green',label='$v^{ini}_{rot}$=200 $km~s^{-1}$',lw=1)
#plt.plot(star_age_250, center_h1_250,linestyle='-',color='blue',label='$v^{ini}_{rot}$=250 $km~s^{-1}$',lw=2)

plt.savefig(pp1_all_res, format='pdf')

pp1_all_res.close()

'''

'''
pp1_all_res_vrot= PdfPages('age_vsurf_comparison_onestar_resolution_vrot.pdf') 


plt.figure(figsize=(10, 10))
plt.title("$\it{M}_\mathrm{ini}$ = 24 $\it{M}_\mathrm{\odot}$", fontsize=30)
plt.plot(star_age_1[0:iRLOF_1]/1e6, surf_avg_vtor_1[0:iRLOF_1],color='red',label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=0.5, time delta=0.75',lw=8,linestyle='--')

plt.plot(star_age_200[0:iRLOF_200]/1e6, surf_avg_vtor_200[0:iRLOF_200],label='$\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, mesh delta=1, time delta=1',c= colors_list[3])

plt.ylabel('Surface rotational velocity [$km~s^{-1}$]',fontsize=30)
plt.xlabel('Star age [Myrs]',fontsize=30)
plt.legend(loc=1,fontsize=16)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
plt.xlim([-1,10.5])
plt.savefig(pp1_all_res_vrot, format='pdf')
pp1_all_res_vrot.close()

'''






