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





m_200=mp.MESA()
m2_200=mp.MESA()
m3_200=mp.MESA()

m2_200_pos=mp.MESA()
m3_200_pos=mp.MESA()

m2_200_pos_old=mp.MESA()
m3_200_pos_old=mp.MESA()


mHD_200=mp.MESA()
m2HD_200=mp.MESA()

m3HD_200=mp.MESA()
m3HD_pos_200=mp.MESA()


m3HD25_200=mp.MESA()
m3HD25_pos_200=mp.MESA()

m3HD19_200=mp.MESA()
m3HD19_pos_200=mp.MESA()

p1=mp.MESA()
p2=mp.MESA()

p1_pos=mp.MESA()
p2_pos=mp.MESA()

profiles_1_1_mesa_46_pos = mp.MESA()
profiles_1_1_mesa_46=mp.MESA()
profiles_1_1_mesa_25_pos=mp.MESA()
profiles_1_1_mesa_25=mp.MESA()
profiles_1_1_mesa_19=mp.MESA()
profiles_1_1_mesa_19_pos=mp.MESA()

m3HD_200.log_fold=os.path.join(paths.data,'CE/LOGS3_HD46485')
m3HD_200.loadHistory()
m3HD_pos_200.log_fold=os.path.join(paths.data,'CE/LOGS3_HD46485_pos')
m3HD_pos_200.loadHistory()

m3HD25_200.log_fold=os.path.join(paths.data,'CE/LOGS3_HD25631')
m3HD25_200.loadHistory()
m3HD25_pos_200.log_fold=os.path.join(paths.data,'CE/LOGS3_HD25631_pos')
m3HD25_pos_200.loadHistory()

m3HD19_200.log_fold=os.path.join(paths.data,'CE/LOGS3_HD191495')
m3HD19_200.loadHistory()
m3HD19_pos_200.log_fold=os.path.join(paths.data,'CE/LOGS3_HD191495_pos')
m3HD19_pos_200.loadHistory()


'''
profiles_1_1_mesa_46 = MesaData(os.path.join(paths.data,'CE/LOGS1_HD46485/profile8.data'))
profiles_1_1_mesa_46_pos = MesaData(os.path.join(paths.data,'CE/LOGS1_HD46485_pos/profile8.data'))


#profiles_1_1_mesa_25 = MesaData('LOGS1_HD25631/profile14.data')
#profiles_1_1_mesa_25_pos = MesaData('LOGS1_HD25631_pos/profile26.data')

profiles_1_1_mesa_25 = MesaData(os.path.join(paths.data,'CE/LOGS1_HD25631/profile3.data'))
profiles_1_1_mesa_25_pos = MesaData(os.path.join(paths.data,'CE/LOGS1_HD25631_pos/profile3.data'))



profiles_1_1_mesa_19 = MesaData(os.path.join(paths.data,'CE/LOGS1_HD191495/profile10.data'))
profiles_1_1_mesa_19_pos = MesaData(os.path.join(paths.data,'CE/LOGS1_HD191495_pos/profile4.data'))
'''

#profiles_1_1_mesa_46.log_fold=os.path.join(paths.data,'CE/LOGS1_HD46485/profile8.data')
profiles_1_1_mesa_46.loadProfile(os.path.join(paths.data,'CE/LOGS1_HD46485/profile8.data'))
#profiles_1_1_mesa_46_pos.log_fold=os.path.join(paths.data,'CE/LOGS1_HD46485_pos/profile8.data')
profiles_1_1_mesa_46_pos.loadProfile(os.path.join(paths.data,'CE/LOGS1_HD46485_pos/profile8.data'))

#profiles_1_1_mesa_25 = MesaData('LOGS1_HD25631/profile14.data')
#profiles_1_1_mesa_25_pos = MesaData('LOGS1_HD25631_pos/profile26.data')

#profiles_1_1_mesa_25.log_fold=os.path.join(paths.data,'CE/LOGS1_HD25631/profile3.data')
profiles_1_1_mesa_25.loadProfile(os.path.join(paths.data,'CE/LOGS1_HD25631/profile3.data'))
#profiles_1_1_mesa_25_pos.log_fold=os.path.join(paths.data,'CE/LOGS1_HD25631_pos/profile3.data')
profiles_1_1_mesa_25_pos.loadProfile(os.path.join(paths.data,'CE/LOGS1_HD25631_pos/profile3.data'))
#profiles_1_1_mesa_19.log_fold=os.path.join(paths.data,'CE/LOGS1_HD191495/profile10.data')
profiles_1_1_mesa_19.loadProfile(os.path.join(paths.data,'CE/LOGS1_HD191495/profile10.data'))
#profiles_1_1_mesa_19_pos.log_fold=os.path.join(paths.data,'CE/LOGS1_HD191495_pos/profile4.data')
profiles_1_1_mesa_19_pos.loadProfile(os.path.join(paths.data,'CE/LOGS1_HD191495_pos/profile4.data'))





def get_BE_from_pfile(profile_mesa, alpha_th=1.0, alpha_rot=0.0):
    """Calculates the binding energy profile of the star. See Eq. 6 in
    Dewi & Tauris 2000 (but change sign, binding energy>0 if layer is bound).
    The binding energy is calculated integrating the potential plus a
    fraction `alpha_th` of the internal energy.
    Parameters
    ----------
    pfile    : `MESA profile*.data` file
               assumed to contain the columns mass, energy, radius, and dm
    alpha_th : `float` optional, default 1
               fraction of thermal energy to include in the BE
    Returns
    -------
    BE       : `np.array` binding energy in cgs units
    """
    # from time import time
    # get the data in cgs units


    G_cgs = 6.67430e-8  # in cgs
    mu_sun = 1.3271244e26
    Lsun = 3.828e33
    Msun = mu_sun / G_cgs
    Rsun_cm = 6.957e10  # in cm
    clight = 2.99792458e10  # cm/s

    #In [16]: vcrit_46485=np.sqrt( (mu_sun * 24) / (Rsun_cm * 11) )/100/1000

    #In [17]: vcrit_25631=np.sqrt( (mu_sun * 7.5) / (Rsun_cm * 4.2) )/100/1000

    #In [18]: vcrit_191495=np.sqrt( (mu_sun * 15) / (Rsun_cm * 6.8) )/100/1000
    '''
    m = profile_mesa.data('mass') * Msun
    dm = profile_mesa.data('dm')
    r = profile_mesa.data('radius') * Rsun_cm
    u = profile_mesa.data('energy')
    omega = profile_mesa.data('omega')
    '''


    m = profile_mesa.prof.mass * Msun
    dm = profile_mesa.prof.dm
    r = profile_mesa.prof.radius * Rsun_cm
    u = profile_mesa.prof.energy
    omega = profile_mesa.prof.omega






    #src, col = getSrcCol(pfile)
    #m = src[:, col.index("mass")] * Msun  # g
    #dm = src[:, col.index("dm")]  # g
    #r = src[:, col.index("radius")] * Rsun_cm  # cm
    #u = src[:, col.index("energy")]  # erg
    # calculate local gravitational potential
    psi = -1.0 * G_cgs * np.divide(m, r)  # erg
    if alpha_rot != 0:
        print('check,check,check')
        # calculate rotationa energy of spherical shell of
        # mass dm, outer radius r,
        # thickness dr, and rotation frequency omega
        # omega = src[:, omega]  # 1/sec
        I = (2.0 / 3.0) * np.square(r)  # specific moment of inertia cgs units
        erot = 0.5 * I * np.square(omega)
        #print('omega',omega)
        #print('I',I)
        #print('erot',erot)
        #print('psi',psi)
        #print('u',u)

    else:
        erot = np.zeros(len(psi))
    # change sign: BE is the energy (>0) to provide to unbind the star
    BE =  -1.0 * np.cumsum(np.multiply(psi + alpha_th * u + alpha_rot * erot, dm))

    return np.asarray(BE, dtype=float)


alpha_th=1.0
alpha_rot=0.0

G_cgs = 6.67430e-8  # in cgs
mu_sun = 1.3271244e26
Lsun = 3.828e33
Msun = mu_sun / G_cgs
Rsun_cm = 6.957e10  # in cm
clight = 2.99792458e10  # cm/s





#HD46485

profile_fin=profiles_1_1_mesa_46
profile_fin_pos=profiles_1_1_mesa_46_pos

print('profile:', profile_fin)
print('profile_pos:', profile_fin_pos)

BE_1 = get_BE_from_pfile(profile_fin, alpha_th=alpha_th, alpha_rot=alpha_rot)
BE_1_pos = get_BE_from_pfile(profile_fin_pos, alpha_th=alpha_th, alpha_rot=alpha_rot)

BE_0 = get_BE_from_pfile(profile_fin, alpha_th=0, alpha_rot=alpha_rot)
#BE_0_pos = get_BE_from_pfile(profile_fin_pos, alpha_th=1, alpha_rot=1)

#print('BE_0',BE_0)
#print('BE_0_pos',BE_0_pos)

#print('BE_1',BE_1)
#print('BE_1_pos',BE_1_pos)


age_HD=m3HD_200.hist.age

age_HD_pos=m3HD_pos_200.hist.age


star_age_ind=find_nearest(age_HD,7788538.73)
star_age_ind_pos=find_nearest(age_HD_pos,7920428.415)

age_fitted=age_HD[star_age_ind]
age_fitted_pos=age_HD_pos[star_age_ind_pos]

print('age_fitted',age_fitted)


star_1_mass_HD=m3HD_200.hist.star_1_mass[star_age_ind]
star_2_mass_HD=m3HD_200.hist.star_2_mass[star_age_ind]
binary_separation_HD=m3HD_200.hist.binary_separation[star_age_ind]

print('binary_sepr_HD46485',binary_separation_HD)

star_1_mass_HD_pos=m3HD_pos_200.hist.star_1_mass[star_age_ind_pos]
star_2_mass_HD_pos=m3HD_pos_200.hist.star_2_mass[star_age_ind_pos]
binary_separation_HD_pos=m3HD_pos_200.hist.binary_separation[star_age_ind_pos]



#v_1_HD=m3HD_200.hist.v_orb_1[star_age_ind]
#v_2_HD=m3HD_200.hist.v_orb_2[star_age_ind]


mass_profile = profile_fin.prof.mass
r_profile = 10**profile_fin.prof.logR
h1_profile = profile_fin.prof.h1



mass_profile_pos = profile_fin_pos.prof.mass
r_profile_pos = 10**profile_fin_pos.prof.logR
h1_profile_pos = profile_fin_pos.prof.h1




ORB_E =   G_cgs * mass_profile * star_2_mass_HD * Msun**2  / (2 * r_profile * Rsun_cm) - (G_cgs * star_1_mass_HD * star_2_mass_HD * Msun**2) / ( 2* binary_separation_HD * Rsun_cm)

ORB_E_pos =   G_cgs * mass_profile_pos * star_2_mass_HD_pos * Msun**2  / (2 * r_profile_pos * Rsun_cm) - (G_cgs * star_1_mass_HD_pos * star_2_mass_HD_pos * Msun**2) / ( 2* binary_separation_HD_pos * Rsun_cm)

pp1_be = PdfPages(paths.figures / 'r_alpha_HD46485.pdf') 

#gs = gridspec.GridSpec(2, 2)

fig = plt.figure(figsize=(10, 10))
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

#ax=fig.add_subplot(gs[0,0])
#mpl.rcParams.update({'font.size': 8})
plt.tick_params(labelsize=8)

plt.title('HD46485, $\it{M}_\mathrm{1,ini}$=24, $\it{M}_\mathrm{2,ini}$=1',fontsize=30)
#plt.scatter(profile_fin.data('radius'), profile_fin.data('energy'), c='b', marker='x', lw=5, s=15,label='BE_1')
#plt.scatter(profile_fin.data('radius'), psi_int, c='b', marker='x', lw=5, s=15,label='psi_int')
star_point=5000
print('len',len(BE_1))
print(profile_fin.prof.radius[:-10])
print(BE_1[:-10]/ORB_E[:-10])
ax1.plot(profile_fin.prof.radius[:-50], BE_1[:-50]/ORB_E[:-50], c='b', linestyle='-',label=r'$\alpha_\mathrm{th}$=1, MESA default')
ax1.plot(profile_fin.prof.radius[:-50], BE_0[:-50]/ORB_E[:-50], c='r',linestyle='-',label=r'$\alpha_\mathrm{th}$=0, MESA default')
ax1.plot(profile_fin_pos.prof.radius[:-50], BE_1_pos[:-50]/ORB_E_pos[:-50], c='green',linestyle='-',label=r'$\alpha_\mathrm{th}$=1, POSYDON')
#ax1.plot(profile_fin_pos.data('radius'), BE_0_pos/ORB_E_pos, c='black',linestyle='-', lw=2 ,label='alpha_th=1, alpha_rot=1, posydon')

ax2.plot(profile_fin.prof.radius, h1_profile, c='gray',linestyle='--', lw=1,label='h1')
#plt.gca().twinx().plot(h1_profile, color = 'gray') 
ax2.set_ylabel('$X(^{1}\mathrm{H})$',fontsize=20)
ax1.set_ylabel(r'$\alpha_\mathrm{CE}$',fontsize=20)
ax1.set_yscale('log')


left, bottom, width, height = [0.33, 0.28, 0.25, 0.25]
ax3 = fig.add_axes([left, bottom, width, height])
ax4 = ax3.twinx()
ax3.plot(profile_fin.prof.radius, BE_1/ORB_E, c='b',linestyle='-', lw=2 ,label='alpha_th=1')
ax3.plot(profile_fin.prof.radius, BE_0/ORB_E, c='r', linestyle='-', lw=2 , label='alpha_th=0')
ax3.plot(profile_fin_pos.prof.radius, BE_1_pos/ORB_E_pos, c='green',linestyle='-', lw=2 ,label='alpha_th=1, posydon')
#ax3.plot(profile_fin_pos.data('radius'), BE_0_pos/ORB_E_pos, c='black',linestyle='-', lw=2 ,label='alpha_th=1, alpha_rot=1, posydon')

ax4.plot(profile_fin.prof.radius, h1_profile, c='gray',linestyle='--', lw=1,label='h1')
#ax3.set_yscale('log')
ax3.tick_params(labelsize=12)
ax4.tick_params(labelsize=12)


ax3.set_ylim([0,30])
ax3.set_xlim([1,4])

ax1.set_xlabel('Radius [$\it{R}_{\odot}$]',fontsize=20)
ax1.legend(loc=4,fontsize=15,facecolor='white', framealpha=1,frameon=True)
ax2.legend(loc=1,fontsize=15,facecolor='white', framealpha=1,frameon=True)
ax1.tick_params(labelsize=18)
ax2.tick_params(labelsize=18)


plt.savefig(pp1_be, format='pdf')

pp1_be.close()




















#HD191495

profile_fin=profiles_1_1_mesa_19
profile_fin_pos=profiles_1_1_mesa_19_pos

print('profile:', profile_fin)
print('profile_pos:', profile_fin_pos)

BE_1 = get_BE_from_pfile(profile_fin, alpha_th=alpha_th, alpha_rot=alpha_rot)
BE_1_pos = get_BE_from_pfile(profile_fin_pos, alpha_th=alpha_th, alpha_rot=alpha_rot)

BE_0 = get_BE_from_pfile(profile_fin, alpha_th=0, alpha_rot=alpha_rot)
#BE_0_pos = get_BE_from_pfile(profile_fin_pos, alpha_th=1, alpha_rot=1)

#print('BE_0',BE_0)
#print('BE_0_pos',BE_0_pos)

#print('BE_1',BE_1)
#print('BE_1_pos',BE_1_pos)


age_HD=m3HD19_200.hist.age
age_HD_pos=m3HD19_pos_200.hist.age


star_age_ind=find_nearest(age_HD,5855037.917)
star_age_ind_pos=find_nearest(age_HD_pos,12421615.4)

age_fitted=age_HD[star_age_ind]
age_fitted_pos=age_HD_pos[star_age_ind_pos]

print('age_fitted',age_fitted)


star_1_mass_HD=m3HD19_200.hist.star_1_mass[star_age_ind]
star_2_mass_HD=m3HD19_200.hist.star_2_mass[star_age_ind]
binary_separation_HD=m3HD19_200.hist.binary_separation[star_age_ind]

print('binary_sepr_HD191495',binary_separation_HD)
star_1_mass_HD_pos=m3HD19_pos_200.hist.star_1_mass[star_age_ind_pos]
star_2_mass_HD_pos=m3HD19_pos_200.hist.star_2_mass[star_age_ind_pos]
binary_separation_HD_pos=m3HD19_pos_200.hist.binary_separation[star_age_ind_pos]



#v_1_HD=m3HD_200.hist.v_orb_1[star_age_ind]
#v_2_HD=m3HD_200.hist.v_orb_2[star_age_ind]


mass_profile = profile_fin.prof.mass
r_profile = 10**profile_fin.prof.logR
h1_profile = profile_fin.prof.h1



mass_profile_pos = profile_fin_pos.prof.mass
r_profile_pos = 10**profile_fin_pos.prof.logR
h1_profile_pos = profile_fin_pos.prof.h1




ORB_E =   G_cgs * mass_profile * star_2_mass_HD * Msun**2  / (2 * r_profile * Rsun_cm) - (G_cgs * star_1_mass_HD * star_2_mass_HD * Msun**2) / ( 2* binary_separation_HD * Rsun_cm)

ORB_E_pos =   G_cgs * mass_profile_pos * star_2_mass_HD_pos * Msun**2  / (2 * r_profile_pos * Rsun_cm) - (G_cgs * star_1_mass_HD_pos * star_2_mass_HD_pos * Msun**2) / ( 2* binary_separation_HD_pos * Rsun_cm)

pp1_be = PdfPages(paths.figures / 'r_alpha_HD191495.pdf') 

#gs = gridspec.GridSpec(2, 2)

fig = plt.figure(figsize=(10, 10))
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

#ax=fig.add_subplot(gs[0,0])
#mpl.rcParams.update({'font.size': 8})
plt.tick_params(labelsize=8)

plt.title('HD191495, $\it{M}_\mathrm{1,ini}$=15, $\it{M}_\mathrm{2,ini}$=1.5',fontsize=30)
#plt.scatter(profile_fin.data('radius'), profile_fin.data('energy'), c='b', marker='x', lw=5, s=15,label='BE_1')
#plt.scatter(profile_fin.data('radius'), psi_int, c='b', marker='x', lw=5, s=15,label='psi_int')
star_point=5000
print('len',len(BE_1))
print(profile_fin.prof.radius[:-10])
print(BE_1[:-10]/ORB_E[:-10])
ax1.plot(profile_fin.prof.radius[:-50], BE_1[:-50]/ORB_E[:-50], c='b', linestyle='-',label=r'$\alpha_\mathrm{th}$=1, MESA default')
ax1.plot(profile_fin.prof.radius[:-50], BE_0[:-50]/ORB_E[:-50], c='r',linestyle='-',label=r'$\alpha_\mathrm{th}$=0, MESA default')
ax1.plot(profile_fin_pos.prof.radius[:-50], BE_1_pos[:-50]/ORB_E_pos[:-50], c='green',linestyle='-',label=r'$\alpha_\mathrm{th}$=1, POSYDON')
#ax1.plot(profile_fin_pos.data('radius'), BE_0_pos/ORB_E_pos, c='black',linestyle='-', lw=2 ,label='alpha_th=1, alpha_rot=1, POSYDON')

ax2.plot(profile_fin.prof.radius, h1_profile, c='gray',linestyle='--', lw=1,label='h1')
#plt.gca().twinx().plot(h1_profile, color = 'gray') 
ax2.set_ylabel('$X(^{1}\mathrm{H})$',fontsize=20)
ax1.set_ylabel(r'$\alpha_\mathrm{CE}$',fontsize=20)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax1.set_yscale('log')


left, bottom, width, height = [0.35, 0.28, 0.25, 0.25]
ax3 = fig.add_axes([left, bottom, width, height])
ax4 = ax3.twinx()
ax3.plot(profile_fin.prof.radius, BE_1/ORB_E, c='b',linestyle='-', lw=2 ,label='alpha_th=1, MESA default')
ax3.plot(profile_fin.prof.radius, BE_0/ORB_E, c='r', linestyle='-', lw=2 , label=r'$\alpha_\mathrm{th}$=0, MESA default')
ax3.plot(profile_fin_pos.prof.radius, BE_1_pos/ORB_E_pos, c='green',linestyle='-', lw=2 ,label=r'$\alpha_\mathrm{th}$, POSYDON')
#ax3.plot(profile_fin_pos.data('radius'), BE_0_pos/ORB_E_pos, c='black',linestyle='-', lw=2 ,label='alpha_th=1, alpha_rot=1, posydon')

ax4.plot(profile_fin.prof.radius, h1_profile, c='gray',linestyle='--', lw=1,label='h1')
#ax3.set_yscale('log')
ax3.tick_params(labelsize=12)
ax4.tick_params(labelsize=12)


ax3.set_ylim([0,30])
ax3.set_xlim([1,4])

ax1.set_xlabel('Radius [$\it{R}_{\odot}$]',fontsize=20)
ax1.legend(loc=4,fontsize=15,facecolor='white', framealpha=1,frameon=True)
ax2.legend(loc=1,fontsize=15,facecolor='white', framealpha=1,frameon=True)
ax1.tick_params(labelsize=18)
ax2.tick_params(labelsize=18)

plt.savefig(pp1_be, format='pdf')

pp1_be.close()













#HD25631

profile_fin=profiles_1_1_mesa_25
profile_fin_pos=profiles_1_1_mesa_25_pos

print('profile:', profile_fin)
print('profile_pos:', profile_fin_pos)

BE_1 = get_BE_from_pfile(profile_fin, alpha_th=alpha_th, alpha_rot=alpha_rot)
BE_1_pos = get_BE_from_pfile(profile_fin_pos, alpha_th=alpha_th, alpha_rot=alpha_rot)

BE_0 = get_BE_from_pfile(profile_fin, alpha_th=0, alpha_rot=alpha_rot)
#BE_0_pos = get_BE_from_pfile(profile_fin_pos, alpha_th=1, alpha_rot=1)

#print('BE_0',BE_0)
#print('BE_0_pos',BE_0_pos)

#print('BE_1',BE_1)
#print('BE_1_pos',BE_1_pos)


age_HD=m3HD25_200.hist.age
age_HD_pos=m3HD25_pos_200.hist.age


#star_age_ind=find_nearest(age_HD,56000632.84) #last profile
star_age_ind=find_nearest(age_HD,55600000.00) # TAMS
#star_age_ind_pos=find_nearest(age_HD_pos,55983531.15)
star_age_ind_pos=find_nearest(age_HD_pos,55600000.00)

age_fitted=age_HD[star_age_ind]
age_fitted_pos=age_HD_pos[star_age_ind_pos]

print('age_fitted',age_fitted)
print('age_fitted_pos',age_fitted_pos)


star_1_mass_HD=m3HD25_200.hist.star_1_mass[star_age_ind]
star_2_mass_HD=m3HD25_200.hist.star_2_mass[star_age_ind]
binary_separation_HD=m3HD25_200.hist.binary_separation[star_age_ind]


print('binary_sepr_HD25631',binary_separation_HD)


star_1_mass_HD_pos=m3HD25_pos_200.hist.star_1_mass[star_age_ind_pos]
star_2_mass_HD_pos=m3HD25_pos_200.hist.star_2_mass[star_age_ind_pos]
binary_separation_HD_pos=m3HD25_pos_200.hist.binary_separation[star_age_ind_pos]



#v_1_HD=m3HD_200.hist.v_orb_1[star_age_ind]
#v_2_HD=m3HD_200.hist.v_orb_2[star_age_ind]


mass_profile = profile_fin.prof.mass
r_profile = 10**profile_fin.prof.logR
h1_profile = profile_fin.prof.h1



mass_profile_pos = profile_fin_pos.prof.mass
r_profile_pos = 10**profile_fin_pos.prof.logR
h1_profile_pos = profile_fin_pos.prof.h1




ORB_E =   G_cgs * mass_profile * star_2_mass_HD * Msun**2  / (2 * r_profile * Rsun_cm) - (G_cgs * star_1_mass_HD * star_2_mass_HD * Msun**2) / ( 2* binary_separation_HD * Rsun_cm)

ORB_E_pos =   G_cgs * mass_profile_pos * star_2_mass_HD_pos * Msun**2  / (2 * r_profile_pos * Rsun_cm) - (G_cgs * star_1_mass_HD_pos * star_2_mass_HD_pos * Msun**2) / ( 2* binary_separation_HD_pos * Rsun_cm)








pp1_be = PdfPages(paths.figures / 'r_alpha_HD25631.pdf') 

#gs = gridspec.GridSpec(2, 2)

fig = plt.figure(figsize=(10, 10))
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

#ax=fig.add_subplot(gs[0,0])
#mpl.rcParams.update({'font.size': 8})
plt.tick_params(labelsize=8)

plt.title("HD25631, $\it{M}_\mathrm{1,ini}$=7, $\it{M}_\mathrm{2,ini}$=1", fontsize=30)
#plt.scatter(profile_fin.data('radius'), profile_fin.data('energy'), c='b', marker='x', lw=5, s=15,label='BE_1')
#plt.scatter(profile_fin.data('radius'), psi_int, c='b', marker='x', lw=5, s=15,label='psi_int')
star_point=5000
print('len',len(BE_1))
print(profile_fin.prof.radius[:-10])
print(BE_1[:-10]/ORB_E[:-10])
ax1.plot(profile_fin.prof.radius[:-50], BE_1[:-50]/ORB_E[:-50], c='b', linestyle='-',label=r'$\alpha_\mathrm{th}$=1, MESA default')
ax1.plot(profile_fin.prof.radius[:-50], BE_0[:-50]/ORB_E[:-50], c='r',linestyle='-',label=r'$\alpha_\mathrm{th}$=0, MESA default')
ax1.plot(profile_fin_pos.prof.radius[:-50], BE_1_pos[:-50]/ORB_E_pos[:-50], c='green',linestyle='-',label=r'$\alpha_\mathrm{th}$=1, POSYDON')
#ax1.plot(profile_fin_pos.data('radius'), BE_0_pos/ORB_E_pos, c='black',linestyle='-', lw=2 ,label='alpha_th=1, alpha_rot=1, POSYDON')

ax2.plot(profile_fin.prof.radius, h1_profile, c='gray',linestyle='--',label='h1',lw=1)
#plt.gca().twinx().plot(h1_profile, color = 'gray') 
ax2.set_ylabel('$X(^{1}\mathrm{H})$',fontsize=20)
ax1.set_ylabel(r'$\alpha_\mathrm{CE}$',fontsize=20)
ax1.set_yscale('log')


left, bottom, width, height = [0.30, 0.29, 0.25, 0.25]
ax3 = fig.add_axes([left, bottom, width, height])
ax4 = ax3.twinx()
ax3.plot(profile_fin.prof.radius, BE_1/ORB_E, c='b',linestyle='-',label='alpha_th=1')
ax3.plot(profile_fin.prof.radius, BE_0/ORB_E, c='r', linestyle='-', label='alpha_th=0')
ax3.plot(profile_fin_pos.prof.radius, BE_1_pos/ORB_E_pos, c='green',linestyle='-',label='alpha_th=1, POSYDON')
#ax3.plot(profile_fin_pos.data('radius'), BE_0_pos/ORB_E_pos, c='black',linestyle='-', lw=2 ,label='alpha_th=1, alpha_rot=1, posydon')

ax4.plot(profile_fin.prof.radius, h1_profile, c='gray',linestyle='--', lw=1,label='h1')
#ax3.set_yscale('log')
ax3.tick_params(labelsize=12)
ax4.tick_params(labelsize=12)


ax3.set_ylim([0,15])
ax3.set_xlim([0,2])

ax1.set_xlabel('Radius [$\it{R}_{\odot}$]',fontsize=20)
ax1.legend(loc=4,fontsize=15,facecolor='white', framealpha=1,frameon=True)
ax2.legend(loc=1,fontsize=15,facecolor='white', framealpha=1,frameon=True)
ax1.tick_params(labelsize=18)
ax2.tick_params(labelsize=18)
plt.savefig(pp1_be, format='pdf')

pp1_be.close()

