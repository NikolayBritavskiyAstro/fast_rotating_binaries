import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import mesaPlot as mp
from showyourwork.paths import user as Paths

paths = Paths()
plt.style.use(paths.scripts / "matplotlibrc")

if os.path.exists(os.path.join(paths.data, 'eccentricity/HD191495/LOGS3/history.data')):
    pass
else:
    os.system(f'python {os.path.join(paths.scripts / "unzip_MESA_output.py")}')

m_200_191 = mp.MESA()
m2_200_191 = mp.MESA()
m3_200_191 = mp.MESA()

m_200_191_newtides = mp.MESA()
m2_200_191_newtides = mp.MESA()
m3_200_191_newtides = mp.MESA()

m_200_nowind_191 = mp.MESA()
m2_200_nowind_191 = mp.MESA()
m3_200_nowind_191 = mp.MESA()

m_200_46 = mp.MESA()
m2_200_46 = mp.MESA()
m3_200_46 = mp.MESA()

m_200_46_newtides = mp.MESA()
m2_200_46_newtides = mp.MESA()
m3_200_46_newtides = mp.MESA()

m_200_nowind_46 = mp.MESA()
m2_200_nowind_46 = mp.MESA()
m3_200_nowind_46 = mp.MESA()

m_200_25 = mp.MESA()
m2_200_25 = mp.MESA()
m3_200_25 = mp.MESA()

m_200_25_newtides = mp.MESA()
m2_200_25_newtides = mp.MESA()
m3_200_25_newtides = mp.MESA()

m_200_nowind_25 = mp.MESA()
m2_200_nowind_25 = mp.MESA()
m3_200_nowind_25 = mp.MESA()

m3_200_191.log_fold = os.path.join(paths.data, 'eccentricity/HD191495/LOGS3')
m3_200_191.loadHistory()

m3_200_46.log_fold = os.path.join(paths.data, 'eccentricity/HD46485/LOGS3')
m3_200_46.loadHistory()

m3_200_25.log_fold = os.path.join(paths.data, 'eccentricity/HD25631/LOGS3')
m3_200_25.loadHistory()

rsun = 696000  # km

age_200_191 = m3_200_191.hist.age
eccentricity_191 = m3_200_191.hist.eccentricity

age_200_46 = m3_200_46.hist.age
eccentricity_46 = m3_200_46.hist.eccentricity

age_200_25 = m3_200_25.hist.age
eccentricity_25 = m3_200_25.hist.eccentricity

pp1_all = PdfPages(paths.figures / 'age_e_comparison_panel.pdf')

fig = plt.figure(figsize=(10, 10))

ax1 = fig.add_subplot(311)
ax1.set_title(
    'HD46485, $\it{v}^\mathrm{ini}_\mathrm{rot}$=350 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=24, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=6.9 d',
    fontdict={'fontsize': 17, 'fontweight': 'medium'})

plt.plot(age_200_46 / 1e6, eccentricity_46, linestyle='-', color='red', label='MESA default')

plt.legend(loc=3, fontsize=16)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylim([-0.05, 0.55])
plt.ylabel('$\it{e}$', fontsize=24)

ax2 = fig.add_subplot(312)

ax2.set_title(
    'HD191495, $\it{v}^\mathrm{ini}_\mathrm{rot}$=200 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=15, $\it{M}_\mathrm{2,ini}$=1.5, $\it{P}_\mathrm{ini}$=3.6 d',
    fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(age_200_191 / 1e6, eccentricity_191, linestyle='-', color='red', label='MESA default')

plt.legend(loc=3, fontsize=16)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$\it{e}$', fontsize=24)

plt.ylim([-0.05, 0.55])

ax3 = fig.add_subplot(313)

ax3.set_title(
    'HD25631, $\it{v}^\mathrm{ini}_\mathrm{rot}$=220 $km~s^{-1}$, $\it{M}_\mathrm{1,ini}$=7, $\it{M}_\mathrm{2,ini}$=1, $\it{P}_\mathrm{ini}$=5.2 d',
    fontdict={'fontsize': 17, 'fontweight': 'medium'})
plt.plot(age_200_25 / 1e6, eccentricity_25, linestyle='-', color='red', label='MESA default')

plt.legend(loc=3, fontsize=16)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Star age [Myrs]', fontsize=24)
plt.ylabel('$\it{e}$', fontsize=24)
plt.ylim([0.3, 0.55])

plt.savefig(pp1_all, format='pdf')

pp1_all.close()
