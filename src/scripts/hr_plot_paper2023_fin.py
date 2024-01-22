#!/usr/bin/env python3

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import mesaPlot as mp
from showyourwork.paths import user as Paths

paths = Paths()

if os.path.exists(os.path.join(paths.data, 'FINAL_TABLE_TO_PLOTS_mar2023.txt')):
    pass
else:
    os.system(f'python {os.path.join(paths.scripts / "unzip_MESA_output.py")}')

plt.style.use(paths.scripts / "matplotlibrc")

m_7 = mp.MESA()
m_8 = mp.MESA()
m_10 = mp.MESA()
m_15 = mp.MESA()
m_20 = mp.MESA()
m_25 = mp.MESA()
m_30 = mp.MESA()
m_40 = mp.MESA()
m_50 = mp.MESA()
m_60 = mp.MESA()

m_7.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_7')
m_7.loadHistory()

m_8.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_8')
m_8.loadHistory()

m_10.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_10')
m_10.loadHistory()

m_15.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_15')
m_15.loadHistory()

m_20.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_20')
m_20.loadHistory()

m_25.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_25')
m_25.loadHistory()

m_30.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_30')
m_30.loadHistory()

m_40.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_40')
m_40.loadHistory()

m_50.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_50')
m_50.loadHistory()

m_60.log_fold = os.path.join(paths.data, 'sample_tracks/LOGS_60')
m_60.loadHistory()

center_h7_1 = m_7.hist.center_h1
LOGL_7 = m_7.hist.log_L
log_Teff_7 = m_7.hist.log_Teff
iRLOF_7 = np.argmax(center_h7_1 < 1e-2)
Lnuc_7 = m_7.hist.log_Lnuc
logl_7_ms = LOGL_7[0:iRLOF_7]
log_Teff_7_ms = log_Teff_7[0:iRLOF_7]
Lnuc_7_ms = Lnuc_7[0:iRLOF_7]
i_zams_7 = np.argmin(abs(logl_7_ms - Lnuc_7_ms) > 1e-3)
print('i_zams_7', i_zams_7)

center_h8_1 = m_8.hist.center_h1
LOGL_8 = m_8.hist.log_L
log_Teff_8 = m_8.hist.log_Teff
iRLOF_8 = np.argmax(center_h8_1 < 1e-2)

Lnuc_8 = m_8.hist.log_Lnuc
logl_8_ms = LOGL_8[0:iRLOF_8]
log_Teff_8_ms = log_Teff_8[0:iRLOF_8]
Lnuc_8_ms = Lnuc_8[0:iRLOF_8]
i_zams_8 = np.argmin(abs(logl_8_ms - Lnuc_8_ms) > 1e-3)
print('i_zams_8', i_zams_8)

center_h10_1 = m_10.hist.center_h1
LOGL_10 = m_10.hist.log_L
log_Teff_10 = m_10.hist.log_Teff
iRLOF_10 = np.argmax(center_h10_1 < 1e-2)

Lnuc_10 = m_10.hist.log_Lnuc
logl_10_ms = LOGL_10[0:iRLOF_10]
log_Teff_10_ms = log_Teff_10[0:iRLOF_10]
Lnuc_10_ms = Lnuc_10[0:iRLOF_10]
i_zams_10 = np.argmin(abs(logl_10_ms - Lnuc_10_ms) > 1e-3)
print('i_zams_10', i_zams_10)

center_h15_1 = m_15.hist.center_h1
LOGL_15 = m_15.hist.log_L
log_Teff_15 = m_15.hist.log_Teff
iRLOF_15 = np.argmax(center_h15_1 < 1e-2)

Lnuc_15 = m_15.hist.log_Lnuc
logl_15_ms = LOGL_15[0:iRLOF_15]
log_Teff_15_ms = log_Teff_10[0:iRLOF_15]
Lnuc_15_ms = Lnuc_15[0:iRLOF_15]
i_zams_15 = np.argmin(abs(logl_15_ms - Lnuc_15_ms) > 4e-4)
print('i_zams_15', i_zams_15)

center_h20_1 = m_20.hist.center_h1
LOGL_20 = m_20.hist.log_L
log_Teff_20 = m_20.hist.log_Teff
iRLOF_20 = np.argmax(center_h20_1 < 1e-2)

Lnuc_20 = m_20.hist.log_Lnuc
logl_20_ms = LOGL_20[0:iRLOF_20]
log_Teff_20_ms = log_Teff_20[0:iRLOF_20]
Lnuc_20_ms = Lnuc_20[0:iRLOF_20]
i_zams_20 = np.argmin(abs(logl_20_ms - Lnuc_20_ms) > 1e-3)
print('i_zams_20', i_zams_20)

center_h25_1 = m_25.hist.center_h1
LOGL_25 = m_25.hist.log_L
log_Teff_25 = m_25.hist.log_Teff
iRLOF_25 = np.argmax(center_h25_1 < 1e-2)

Lnuc_25 = m_25.hist.log_Lnuc
logl_25_ms = LOGL_25[0:iRLOF_25]
log_Teff_25_ms = log_Teff_25[0:iRLOF_25]
Lnuc_25_ms = Lnuc_25[0:iRLOF_25]
i_zams_25 = np.argmin(abs(logl_25_ms - Lnuc_25_ms) > 1e-3)
print('i_zams_25', i_zams_25)

center_h30_1 = m_30.hist.center_h1
LOGL_30 = m_30.hist.log_L
log_Teff_30 = m_30.hist.log_Teff
iRLOF_30 = np.argmax(center_h30_1 < 1e-2)

Lnuc_30 = m_30.hist.log_Lnuc
logl_30_ms = LOGL_30[0:iRLOF_30]
log_Teff_30_ms = log_Teff_30[0:iRLOF_30]
Lnuc_30_ms = Lnuc_30[0:iRLOF_30]
i_zams_30 = np.argmin(abs(logl_30_ms - Lnuc_30_ms) > 1e-3)
print('i_zams_30', i_zams_30)

center_h40_1 = m_40.hist.center_h1
LOGL_40 = m_40.hist.log_L
log_Teff_40 = m_40.hist.log_Teff
iRLOF_40 = np.argmax(center_h40_1 < 1e-2)

Lnuc_40 = m_40.hist.log_Lnuc
logl_40_ms = LOGL_40[0:iRLOF_40]
log_Teff_40_ms = log_Teff_40[0:iRLOF_40]
Lnuc_40_ms = Lnuc_40[0:iRLOF_40]
i_zams_40 = np.argmin(abs(logl_40_ms - Lnuc_40_ms) > 1e-3)
print('i_zams_40', i_zams_40)

center_h50_1 = m_50.hist.center_h1
LOGL_50 = m_50.hist.log_L
log_Teff_50 = m_50.hist.log_Teff
iRLOF_50 = np.argmax(center_h50_1 < 5e-2)

Lnuc_50 = m_50.hist.log_Lnuc
logl_50_ms = LOGL_50[0:iRLOF_50]
log_Teff_50_ms = log_Teff_10[0:iRLOF_50]
Lnuc_50_ms = Lnuc_50[0:iRLOF_50]
i_zams_50 = np.argmin(abs(logl_50_ms - Lnuc_50_ms) > 1e-3)
print('i_zams_50', i_zams_50)

center_h60_1 = m_60.hist.center_h1
LOGL_60 = m_60.hist.log_L
log_Teff_60 = m_60.hist.log_Teff
iRLOF_60 = np.argmax(center_h60_1 < 6e-2)

Lnuc_60 = m_60.hist.log_Lnuc
logl_60_ms = LOGL_60[0:iRLOF_60]
log_Teff_60_ms = log_Teff_60[0:iRLOF_60]
Lnuc_60_ms = Lnuc_60[0:iRLOF_60]
i_zams_60 = np.argmin(abs(logl_60_ms - Lnuc_60_ms) > 1e-3)
print('i_zams_60', i_zams_60)

arr = []
arr_rot = []

gon_rot = []
inp = open(os.path.join(paths.data, 'FINAL_TABLE_TO_PLOTS_mar2023.txt'), "r")
lines = inp.readlines()[1:]
g = 0
for line in lines:
    numbers = line.split()
    gon_rot.append(numbers)
    g = g + 1

Teff_Sol = 5778

# Slow Rotators parameters
vsin_gon_low = []
pp_gon_low = []

vsin_gon_lot = []
pp_gon_lot = []

pp_s_gon_total = []
vsini_sb1_slow = []
pp_s_sb1_slow = []
vsini_sb1q_slow = []
pp_s_sb1q_slow = []
vsini_sb1p_slow = []
pp_s_sb1p_slow = []
v_sini_gon_total = []
dv_sini_gon_total = []
n_sp = []
teff_gon_total = []
dteff_gon_total = []
l_gon_total = []
dl_gon_total = []

teff_gon_total_sb1 = []
teff_gon_total_sb1p = []
teff_gon_total_sb1q = []

vsini_sb2p_slow = []
pp_s_sb2p_slow = []
teff_gon_total_sb2p = []

# Fast rotators parameters:

PPtot = []
err_PPtot = []
PP_cl = []
err_PPcl = []

PPfit_cl = []

vsini = []
vsini_err = []

PPtot_norm = []
err_PPtot_norm = []

name_star = []
PPcl_norm = []
err_PPcl_norm = []
err_PPcl_global_norm = []
pp_s = []
err_pp_s = []
spt = []
notes = []
fast_rot = []
teff_work = []

pp_s_sb1 = []
err_pp_s_sb1 = []
spt_sb1 = []
notes_sb1 = []
vsini_sb1 = []

vsini_sb1p = []
pp_s_sb1p = []
err_pp_s_sb1p = []
spt_sb1p = []
notes_sb1p = []

vsini_sb1q = []
pp_s_sb1q = []
err_pp_s_sb1q = []
spt_sb1q = []
notes_sb1q = []

vsini_sb1f = []
pp_s_sb1f = []
err_pp_s_sb1f = []
spt_sb1f = []
notes_sb1f = []
teff_sb1 = []
teff_sb1p = []
teff_sb1q = []

vsini_fast_rot = []
vsini_fast_rot_he_gon = []

pp_s_fast_rot = []
teff_run = []

vsini_fast_rotq = []
pp_s_fast_rotq = []

vsini_eb = []
pp_s_eb = []
teff_eb = []

vsini_sb2p_fast = []

pp_s_sb2p_fast = []

teff_gon_total_sb2p_fast = []

vsin_gon_low_fast = []
pp_gon_low_fast = []
vsin_gon_lot_fast = []
pp_gon_lot_fast = []

l_gon_total_sb2p_fast = []
l_work = []
l_sb1 = []
l_sb1q = []
l_sb1p = []
l_run = []
l_eb = []
lumc = []

vsini_abun = []

pp_s_abun = []
c_abun = []
n_abun = []
o_abun = []
n_o_abun = []
c_n_abun = []
he_abun = []

vsin_gon_bh = []
pp_gon_bh = []
teff_gon_bh = []
dteff_gon_bh = []
l_gon_bh = []
dl_gon_bh = []

vsin_gon_bh1 = []
pp_gon_bh1 = []
teff_gon_bh1 = []
dteff_gon_bh1 = []
l_gon_bh1 = []
dl_gon_bh1 = []

he_gon = []
he_gon_err = []

he_gon_total = []
he_gon_total_err = []

he_gon_total_sb1 = []
he_gon_total_sb1q = []
teff_eb_fast = []
l_eb_fast = []
err_pp_s_sb2p_fast = []

vsini_el = []

pp_s_el = []
teff_el = []
l_el = []
vsini_el_slow = []

pp_s_el_slow = []
v_sini_gon_red = []
pp_s_gon_red = []

v_sini_gon_red_sb1 = []
pp_s_gon_red_sb1 = []

vsini_eleb = []

pp_s_eleb = []
teff_eleb = []
l_eleb = []

vsin_gon_low_fast_new = []
pp_gon_low_fast_new = []

l_gon_total_sb1 = []
l_gon_el = []
teff_gon_el = []

teff_eleb = []
l_eleb = []

teff_run_slow = []
l_run_slow = []

teff_eb_slow = []
l_eb_slow = []

l_real = []
l_real_sb1 = []
l_real_sb2 = []
l_real_run = []
l_real_eb = []
l_real_ev = []

teff_gon_total_real = []
l_gon_total_real = []

teff_eleb_real = []
l_eleb_real = []
teff_gon_el_real = []
l_gon_el_real = []
teff_eb_slow_real = []
l_eb_slow_real = []

teff_gon_total_sb1_real = []
l_gon_total_sb1_real = []

for ii in range(0, g):

    if (gon_rot[ii][9] == 'nan'):

        n_sp.append(float(gon_rot[ii][5]))
        if not math.isnan(float(gon_rot[ii][15])):
            v_sini_gon_total.append(float(gon_rot[ii][8]))

            dv_sini_gon_total.append(float(gon_rot[ii][9]))

            teff_gon_total.append(float(gon_rot[ii][12]))
            dteff_gon_total.append(float(gon_rot[ii][13]))

            l_gon_total.append(float(gon_rot[ii][10]))
            dl_gon_total.append(float(gon_rot[ii][11]))
            pp_s_gon_total.append(float(gon_rot[ii][15]))
            # print('SPECTRA',gon_rot[ii][1],gon_rot[ii][27],gon_rot[ii][12])
            if not gon_rot[ii][27] == "_":
                teff_gon_total_real.append(float(gon_rot[ii][12]))
                l_gon_total_real.append(float(gon_rot[ii][27]))

            if (float(gon_rot[ii][15]) > 30) and (float(gon_rot[ii][8]) < 200):

                v_sini_gon_red.append(float(gon_rot[ii][8]))
                pp_s_gon_red.append(float(gon_rot[ii][15]))

                if gon_rot[ii][19] == 'EL':
                    vsini_el_slow.append(float(gon_rot[ii][8]))

                    pp_s_el_slow.append(float(gon_rot[ii][15]))

                    l_gon_el.append(float(gon_rot[ii][10]))
                    teff_gon_el.append(float(gon_rot[ii][12]))
                    if not gon_rot[ii][27] == "_":
                        teff_gon_el_real.append(float(gon_rot[ii][12]))
                        l_gon_el_real.append(float(gon_rot[ii][27]))

            if (float(gon_rot[ii][15]) > 20) and (float(gon_rot[ii][15]) < 30) and (float(gon_rot[ii][8]) < 200) and (
                    gon_rot[ii][17] == 'SB1'):
                # print(gon_rot[ii][1])
                v_sini_gon_red_sb1.append(float(gon_rot[ii][8]))
                pp_s_gon_red_sb1.append(float(gon_rot[ii][15]))

            if float(gon_rot[ii][5]) < 3:
                vsin_gon_low_fast_new.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_low_fast_new.append(float(gon_rot[ii][15]))

            # if math.isnan((float(gon_rot[ii][15]))) and (float(gon_rot[ii][8]) < 200):
            if (float(gon_rot[ii][15]) > 50) and (float(gon_rot[ii][15]) < 100):
                # print(gon_rot[ii][1],gon_rot[ii][6],float(gon_rot[ii][8]),float(gon_rot[ii][15]))

                print(gon_rot[ii][1], float(gon_rot[ii][8]), float(gon_rot[ii][15]))

            if (float(gon_rot[ii][15]) > 20) and (float(gon_rot[ii][8]) < 200) and (gon_rot[ii][18] == 'M19'):
                # print(gon_rot[ii][1],gon_rot[ii][6],float(gon_rot[ii][8]),float(gon_rot[ii][15]))

                print(gon_rot[ii][1])

            # print('ALL SB1 slow',len(v_sini_gon_red))

            if float(gon_rot[ii][5]) < 3:
                vsin_gon_low.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_low.append(float(gon_rot[ii][15]))
            else:
                vsin_gon_lot.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_lot.append(float(gon_rot[ii][15]))

            if gon_rot[ii][1] == 'HD226868':
                vsin_gon_bh.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_bh.append(float(gon_rot[ii][15]))
                teff_gon_bh.append(float(gon_rot[ii][12]))
                dteff_gon_bh.append(float(gon_rot[ii][13]))

                l_gon_bh.append(float(gon_rot[ii][10]))
                dl_gon_bh.append(float(gon_rot[ii][11]))

            if gon_rot[ii][1] == 'HD130298' or gon_rot[ii][1] == 'HD12323' or gon_rot[ii][1] == 'HD94024':
                vsin_gon_bh1.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_bh1.append(float(gon_rot[ii][15]))
                teff_gon_bh1.append(float(gon_rot[ii][12]))
                dteff_gon_bh1.append(float(gon_rot[ii][13]))

                l_gon_bh1.append(float(gon_rot[ii][10]))
                dl_gon_bh1.append(float(gon_rot[ii][11]))

            if gon_rot[ii][17] == 'SB1':

                vsini_sb1_slow.append(float(gon_rot[ii][8]))

                pp_s_sb1_slow.append(float(gon_rot[ii][15]))

                teff_gon_total_sb1.append(float(gon_rot[ii][12]))
                l_gon_total_sb1.append(float(gon_rot[ii][10]))
                if not gon_rot[ii][27] == "_":
                    teff_gon_total_sb1_real.append(float(gon_rot[ii][12]))
                    l_gon_total_sb1_real.append(float(gon_rot[ii][27]))


            elif gon_rot[ii][17] == 'LPV/SB1?':

                vsini_sb1q_slow.append(float(gon_rot[ii][8]))

                pp_s_sb1q_slow.append(float(gon_rot[ii][15]))

                teff_gon_total_sb1q.append(float(gon_rot[ii][12]))



            elif (gon_rot[ii][17] == 'SB2') or (gon_rot[ii][17] == 'SB2?'):
                vsini_sb1p_slow.append(float(gon_rot[ii][8]))

                pp_s_sb1p_slow.append(float(gon_rot[ii][15]))

                teff_gon_total_sb1p.append(float(gon_rot[ii][12]))


            elif (gon_rot[ii][17] == 'LPV/SB2?'):
                vsini_sb2p_slow.append(float(gon_rot[ii][8]))

                pp_s_sb2p_slow.append(float(gon_rot[ii][15]))

                teff_gon_total_sb2p.append(float(gon_rot[ii][12]))

            if gon_rot[ii][19] == 'EL/EB':
                print('EL/EB slow')
                print(gon_rot[ii][1], float(gon_rot[ii][15]))
                vsini_eleb.append(float(gon_rot[ii][8]))

                pp_s_eleb.append(float(gon_rot[ii][15]))
                teff_eleb.append(float(gon_rot[ii][12]))
                l_eleb.append(float(gon_rot[ii][10]))
                if not gon_rot[ii][27] == "_":
                    teff_eleb_real.append(float(gon_rot[ii][12]))
                    l_eleb_real.append(float(gon_rot[ii][27]))

            if gon_rot[ii][19] == 'EB':
                print('EB slow')
                print(gon_rot[ii][1], float(gon_rot[ii][15]))

                teff_eb_slow.append(float(gon_rot[ii][12]))
                l_eb_slow.append(float(gon_rot[ii][10]))

                if not gon_rot[ii][27] == "_":
                    teff_eb_slow_real.append(float(gon_rot[ii][12]))
                    l_eb_slow_real.append(float(gon_rot[ii][27]))


    else:

        if not math.isnan(float(gon_rot[ii][15])):
            name_star.append(gon_rot[ii][1])

            vsini.append(float(gon_rot[ii][8]))
            vsini_err.append(float(gon_rot[ii][9]))

            pp_s.append(float(gon_rot[ii][15]))
            err_pp_s.append(float(gon_rot[ii][16]))
            spt.append(gon_rot[ii][2])
            lumc.append(gon_rot[ii][3])

            teff_work.append(float(gon_rot[ii][12]))
            l_work.append(float(gon_rot[ii][10]))
            l_real.append(float(gon_rot[ii][27]))

            # print(gon_rot[ii][1])
            he_gon_total.append(float(gon_rot[ii][25]) / 100)
            he_gon_total_err.append(float(gon_rot[ii][26]) / 100)

            if float(gon_rot[ii][6]) < 3:
                vsin_gon_low_fast.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_low_fast.append(float(gon_rot[ii][15]))
            else:
                vsin_gon_lot_fast.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_lot_fast.append(float(gon_rot[ii][15]))

            if (float(gon_rot[ii][15]) > 50) and (float(gon_rot[ii][15]) < 100):
                # print(gon_rot[ii][1],gon_rot[ii][6],float(gon_rot[ii][8]),float(gon_rot[ii][15]))

                print(gon_rot[ii][1], float(gon_rot[ii][8]), float(gon_rot[ii][15]))

            if (gon_rot[ii][17] == 'SB1'):
                vsini_sb1.append(float(gon_rot[ii][8]))

                pp_s_sb1.append(float(gon_rot[ii][15]))
                err_pp_s_sb1.append(float(gon_rot[ii][16]))
                spt_sb1.append(gon_rot[ii][2])
                teff_sb1.append(float(gon_rot[ii][12]))
                l_sb1.append(float(gon_rot[ii][10]))
                he_gon_total_sb1.append(float(gon_rot[ii][25]) / 100)
                l_real_sb1.append(float(gon_rot[ii][27]))


            elif gon_rot[ii][17] == 'LPV/SB1?':

                vsini_sb1q.append(float(gon_rot[ii][8]))

                pp_s_sb1q.append(float(gon_rot[ii][15]))
                err_pp_s_sb1q.append(float(gon_rot[ii][16]))
                spt_sb1q.append(gon_rot[ii][2])
                teff_sb1q.append(float(gon_rot[ii][12]))
                l_sb1q.append(float(gon_rot[ii][10]))
                he_gon_total_sb1q.append(float(gon_rot[ii][25]) / 100)


            elif (gon_rot[ii][17] == 'SB2') or (gon_rot[ii][17] == 'SB2?'):
                vsini_sb1p.append(float(gon_rot[ii][8]))

                pp_s_sb1p.append(float(gon_rot[ii][15]))
                err_pp_s_sb1p.append(float(gon_rot[ii][16]))
                spt_sb1p.append(gon_rot[ii][2])
                teff_sb1p.append(float(gon_rot[ii][12]))
                l_sb1p.append(float(gon_rot[ii][10]))


            elif (gon_rot[ii][17] == 'LPV/SB2?'):
                vsini_sb2p_fast.append(float(gon_rot[ii][8]))
                err_pp_s_sb2p_fast.append(float(gon_rot[ii][16]))

                pp_s_sb2p_fast.append(float(gon_rot[ii][15]))

                teff_gon_total_sb2p_fast.append(float(gon_rot[ii][12]))
                l_gon_total_sb2p_fast.append(float(gon_rot[ii][10]))
                l_real_sb2.append(float(gon_rot[ii][27]))

            if gon_rot[ii][19] == 'EB':
                teff_eb_fast.append(float(gon_rot[ii][12]))
                l_eb_fast.append(float(gon_rot[ii][10]))
                l_real_eb.append(float(gon_rot[ii][27]))

            if gon_rot[ii][19] == 'EL':
                vsini_el.append(float(gon_rot[ii][8]))

                pp_s_el.append(float(gon_rot[ii][15]))
                teff_el.append(float(gon_rot[ii][12]))
                l_el.append(float(gon_rot[ii][10]))
                l_real_ev.append(float(gon_rot[ii][27]))

            if gon_rot[ii][1] == 'HD165174':
                vsin_gon_bh1.append(float(gon_rot[ii][8]))
                # print(gon_rot[ii][8])
                pp_gon_bh1.append(float(gon_rot[ii][15]))
                teff_gon_bh1.append(float(gon_rot[ii][12]))
                dteff_gon_bh1.append(float(gon_rot[ii][13]))

                l_gon_bh1.append(float(gon_rot[ii][10]))
                dl_gon_bh1.append(float(gon_rot[ii][11]))

vsini_fast_rot_slow = []
pp_s_fast_rot_slow = []

vsini_fast_rot_slow_sample = []

pp_s_fast_rot_slow_smple = []
teff_run_slow_sample = []
l_run_slow_sample = []

teff_run_slow_sample_real = []
l_run_slow_sample_real = []

teff_run_slow_real = []
l_run_slow_real = []
teff_sb1_2023 = []
l_real_sb1_2023 = []

for ii in range(0, g):

    if not math.isnan(float(gon_rot[ii][15])):
        if gon_rot[ii][18] == 'M18':
            vsini_fast_rot.append(float(gon_rot[ii][8]))

            pp_s_fast_rot.append(float(gon_rot[ii][15]))
            teff_run.append(float(gon_rot[ii][12]))
            l_run.append(float(gon_rot[ii][10]))
            vsini_fast_rot_he_gon.append(float(gon_rot[ii][25]) / 100)
            l_real_run.append(float(gon_rot[ii][27]))

        if gon_rot[ii][18] == 'M19':
            vsini_fast_rot_slow_sample.append(float(gon_rot[ii][8]))

            pp_s_fast_rot_slow_smple.append(float(gon_rot[ii][15]))
            teff_run_slow_sample.append(float(gon_rot[ii][12]))
            l_run_slow_sample.append(float(gon_rot[ii][10]))

            if not gon_rot[ii][27] == "_":
                teff_run_slow_sample_real.append(float(gon_rot[ii][12]))
                l_run_slow_sample_real.append(float(gon_rot[ii][27]))

        if gon_rot[ii][18] == 'M20':
            vsini_fast_rot_slow.append(float(gon_rot[ii][8]))

            pp_s_fast_rot_slow.append(float(gon_rot[ii][15]))
            teff_run_slow.append(float(gon_rot[ii][12]))
            l_run_slow.append(float(gon_rot[ii][10]))

            if not gon_rot[ii][27] == "_":
                teff_run_slow_real.append(float(gon_rot[ii][12]))
                l_run_slow_real.append(float(gon_rot[ii][27]))

        if gon_rot[ii][1] == "HD46485":
            teff_sb1_2023.append(float(gon_rot[ii][12]))
            l_real_sb1_2023.append(float(gon_rot[ii][27]))
        if gon_rot[ii][1] == "HD25631":
            teff_sb1_2023.append(float(gon_rot[ii][12]))
            l_real_sb1_2023.append(float(gon_rot[ii][27]))
        if gon_rot[ii][1] == "HD191495":
            teff_sb1_2023.append(float(gon_rot[ii][12]))
            l_real_sb1_2023.append(float(gon_rot[ii][27]))

l_real_abun = []
teff_abun = []
for ii in range(0, g):

    if not math.isnan(float(gon_rot[ii][15])):
        if gon_rot[ii][19] == 'EB':
            vsini_eb.append(float(gon_rot[ii][8]))

            pp_s_eb.append(float(gon_rot[ii][15]))
            teff_eb.append(float(gon_rot[ii][12]))
            l_eb.append(float(gon_rot[ii][10]))

for ii in range(0, g):

    if not math.isnan(float(gon_rot[ii][15])):
        if gon_rot[ii][20] != '_':
            vsini_abun.append(float(gon_rot[ii][8]))

            pp_s_abun.append(float(gon_rot[ii][15]))
            c_abun.append(float(gon_rot[ii][20]))
            n_abun.append(float(gon_rot[ii][21]))
            o_abun.append(float(gon_rot[ii][22]))
            n_o_abun.append(float(gon_rot[ii][21]) - float(gon_rot[ii][22]))
            c_n_abun.append(float(gon_rot[ii][20]) - float(gon_rot[ii][21]))
            he_abun.append(float(gon_rot[ii][24]))
            l_real_abun.append(float(gon_rot[ii][27]))
            teff_abun.append(float(gon_rot[ii][12]))

for ii in range(0, g):
    if not math.isnan(float(gon_rot[ii][15])):

        if gon_rot[ii][24] != '_':
            he_gon.append(float(gon_rot[ii][25]) / 100)
            he_gon_err.append(float(gon_rot[ii][26]) / 100)

n_o_abun_sb1 = []
c_n_abun_sb1 = []
he_abun_sb1 = []

n_o_abun_sb1q = []
c_n_abun_sb1q = []
he_abun_sb1q = []

n_o_abun_run = []
c_n_abun_run = []
he_abun_run = []

n_o_abun_single = []
c_n_abun_single = []
he_abun_single = []

vsini_abun_sb1 = []
vsini_abun_run = []
vsini_abun_sb1q = []
vsini_abun_single = []

l_abun_sb1 = []
l_abun_run = []
l_abun_sb1q = []
l_abun_single = []
he_composite_sb1 = []
he_composite_sb1q = []

he_composite_single = []
he_composite_run = []

he_composite_total = []
he_composite_total_err = []
n_o_total = []
n_total = []
n_abun_sb1 = []
n_abun_single = []
n_abun_run = []
names_n = []
for ii in range(0, g):

    if not math.isnan(float(gon_rot[ii][15])):
        if gon_rot[ii][20] != '_':

            he_composite_total.append(float(gon_rot[ii][25]) / 100)
            he_composite_total_err.append(float(gon_rot[ii][26]) / 100)
            n_o_total.append(float(gon_rot[ii][21]) - float(gon_rot[ii][22]))
            n_total.append(float(gon_rot[ii][21]))
            names_n.append((gon_rot[ii][1]))
            if (gon_rot[ii][17] == 'SB1'):
                l_abun_sb1.append(float(gon_rot[ii][10]))
                vsini_abun_sb1.append(float(gon_rot[ii][8]))
                n_abun_sb1.append(float(gon_rot[ii][21]))

                n_o_abun_sb1.append(float(gon_rot[ii][21]) - float(gon_rot[ii][22]))
                c_n_abun_sb1.append(float(gon_rot[ii][20]) - float(gon_rot[ii][21]))
                he_abun_sb1.append(float(gon_rot[ii][24]))
                he_composite_sb1.append(float(gon_rot[ii][25]) / 100)


            elif gon_rot[ii][17] == 'LPV/SB1?':
                l_abun_sb1q.append(float(gon_rot[ii][10]))

                vsini_abun_sb1q.append(float(gon_rot[ii][8]))

                n_o_abun_sb1q.append(float(gon_rot[ii][21]) - float(gon_rot[ii][22]))
                c_n_abun_sb1q.append(float(gon_rot[ii][20]) - float(gon_rot[ii][21]))
                he_abun_sb1q.append(float(gon_rot[ii][24]))
                he_composite_sb1q.append(float(gon_rot[ii][25]) / 100)



            else:
                l_abun_single.append(float(gon_rot[ii][10]))

                vsini_abun_single.append(float(gon_rot[ii][8]))
                n_abun_single.append(float(gon_rot[ii][21]))

                n_o_abun_single.append(float(gon_rot[ii][21]) - float(gon_rot[ii][22]))
                c_n_abun_single.append(float(gon_rot[ii][20]) - float(gon_rot[ii][21]))
                he_abun_single.append(float(gon_rot[ii][24]))

                he_composite_single.append(float(gon_rot[ii][25]) / 100)

for ii in range(0, g):

    if gon_rot[ii][20] != '_':
        if not math.isnan(float(gon_rot[ii][15])):

            if gon_rot[ii][18] == 'M18':
                vsini_abun_run.append(float(gon_rot[ii][8]))
                l_abun_run.append(float(gon_rot[ii][10]))
                n_abun_run.append(float(gon_rot[ii][21]))
                n_o_abun_run.append(float(gon_rot[ii][21]) - float(gon_rot[ii][22]))
                c_n_abun_run.append(float(gon_rot[ii][20]) - float(gon_rot[ii][21]))
                he_abun_run.append(float(gon_rot[ii][24]))
                he_composite_run.append(float(gon_rot[ii][25]) / 100)

            else:
                pass

print('--Slow Rotators--')
print('Total: ', len(v_sini_gon_total))
print('SB1: ', len(vsini_sb1_slow))
print('LPV/SB1?: ', len(vsini_sb1q_slow))
print('SB2-SB2?:', len(vsini_sb1p_slow))
print('LPV/SB2?: ', len(vsini_sb2p_slow))
print('Cyg X1: ', len(vsin_gon_bh))

print('--Fast Rotators--')
print('Total: ', len(vsini))
print('SB1: ', len(vsini_sb1))
print('LPV/SB1?: ', len(vsini_sb1q))
print("SB2-SB2?:", len(vsini_sb1p))
print("LPV/SB2?:", len(vsini_sb2p_fast))

print("EB", len(vsini_eb))
print("Runaways", len(vsini_fast_rot))

pp0 = PdfPages(paths.figures / 'HR_REAL_ROT_Paper2023.pdf')
plt.figure(figsize=(10, 10))

plt.subplots_adjust(bottom=0.08, left=0.09, top=0.92, right=0.91)

ZAMS = [[], []]

ZAMS.append([[10 ** log_Teff_7[i_zams_7] / 1000], [LOGL_7[i_zams_7]]])
ZAMS.append([[10 ** log_Teff_8[i_zams_8] / 1000], [LOGL_8[i_zams_8]]])
ZAMS.append([[10 ** log_Teff_10[i_zams_10] / 1000], [LOGL_10[i_zams_10]]])
ZAMS.append([[10 ** log_Teff_15[i_zams_15] / 1000], [LOGL_15[i_zams_15]]])
ZAMS.append([[10 ** log_Teff_20[i_zams_20] / 1000], [LOGL_20[i_zams_20]]])
ZAMS.append([[10 ** log_Teff_25[i_zams_25] / 1000], [LOGL_25[i_zams_25]]])
ZAMS.append([[10 ** log_Teff_30[i_zams_30] / 1000], [LOGL_30[i_zams_30]]])
ZAMS.append([[10 ** log_Teff_40[i_zams_40] / 1000], [LOGL_40[i_zams_40]]])
ZAMS.append([[10 ** log_Teff_50[i_zams_50] / 1000], [LOGL_50[i_zams_50]]])
ZAMS.append([[10 ** log_Teff_60[i_zams_60] / 1000], [LOGL_60[i_zams_60]]])

listaMasas = ['7', '8', '10', '15', '20', '25', '30', '40', '50', '60']
print(len(listaMasas))
print(ZAMS)
print(ZAMS[0])

for i in range(len(listaMasas)):
    print(ZAMS[i + 2][0][0])
    plt.text(float(ZAMS[i + 2][0][0]) + 1.9, float(ZAMS[i + 2][1][0]) - 0.03, listaMasas[i], fontsize=20)
    if i == 2:
        plt.text(float(ZAMS[i + 2][0][0]) - 2.5, float(ZAMS[i + 2][1][0]) - 0.67, '$\it{M}_\mathrm{\odot}$',
                 fontsize=20)

# plt.plot(10**log_Teff_7[70:iRLOF_7]/1000, LOGL_7[70:iRLOF_7],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_7[i_zams_7:iRLOF_7] / 1000, LOGL_7[i_zams_7:iRLOF_7], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_8[70:iRLOF_8]/1000, LOGL_8[70:iRLOF_8],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_8[i_zams_8:iRLOF_8] / 1000, LOGL_8[i_zams_8:iRLOF_8], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_10[70:iRLOF_10]/1000, LOGL_10[70:iRLOF_10],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_10[i_zams_10:iRLOF_10] / 1000, LOGL_10[i_zams_10:iRLOF_10], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_15[70:iRLOF_15]/1000, LOGL_15[70:iRLOF_15],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_15[i_zams_15:iRLOF_15] / 1000, LOGL_15[i_zams_15:iRLOF_15], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_20[70:iRLOF_20]/1000, LOGL_20[70:iRLOF_20],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_20[i_zams_20:iRLOF_20] / 1000, LOGL_20[i_zams_20:iRLOF_20], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_25[70:iRLOF_25]/1000, LOGL_25[70:iRLOF_25],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_25[i_zams_25:iRLOF_25] / 1000, LOGL_25[i_zams_25:iRLOF_25], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_30[70:iRLOF_30]/1000, LOGL_30[70:iRLOF_30],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_30[i_zams_30:iRLOF_30] / 1000, LOGL_30[i_zams_30:iRLOF_30], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_40[70:iRLOF_40]/1000, LOGL_40[70:iRLOF_40],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_40[i_zams_40:iRLOF_40] / 1000, LOGL_40[i_zams_40:iRLOF_40], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_50[70:iRLOF_50]/1000, LOGL_50[70:iRLOF_50],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_50[i_zams_50:iRLOF_50] / 1000, LOGL_50[i_zams_50:iRLOF_50], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_60[70:iRLOF_60]/1000, LOGL_60[70:iRLOF_60],linestyle='--',c='gray',lw=2)
plt.plot(10 ** log_Teff_60[i_zams_60:iRLOF_60] / 1000, LOGL_60[i_zams_60:iRLOF_60], linestyle='--', c='gray', lw=2)

# plt.plot(10**log_Teff_60/1000, LOGL_60,linestyle='--',c='orange',lw=3)


plt.plot([10 ** log_Teff_60[i_zams_60] / 1000, 10 ** log_Teff_50[i_zams_50] / 1000],
         [LOGL_60[i_zams_60], LOGL_50[i_zams_50]], linestyle='-', color='gray')

plt.plot([10 ** log_Teff_50[i_zams_50] / 1000, 10 ** log_Teff_40[i_zams_40] / 1000],
         [LOGL_50[i_zams_50], LOGL_40[i_zams_40]], linestyle='-', color='gray')

plt.plot([10 ** log_Teff_40[i_zams_40] / 1000, 10 ** log_Teff_30[i_zams_30] / 1000],
         [LOGL_40[i_zams_40], LOGL_30[i_zams_30]], linestyle='-', color='gray')

plt.plot([10 ** log_Teff_30[i_zams_30] / 1000, 10 ** log_Teff_25[i_zams_25] / 1000],
         [LOGL_30[i_zams_30], LOGL_25[i_zams_25]], linestyle='-', color='gray')

plt.plot([10 ** log_Teff_25[i_zams_25] / 1000, 10 ** log_Teff_20[i_zams_20] / 1000],
         [LOGL_25[i_zams_25], LOGL_20[i_zams_20]], linestyle='-', color='gray')

plt.plot([10 ** log_Teff_20[i_zams_20] / 1000, 10 ** log_Teff_10[i_zams_10] / 1000],
         [LOGL_20[i_zams_20], LOGL_10[i_zams_10]], linestyle='-', color='gray')

plt.plot([10 ** log_Teff_10[i_zams_10] / 1000, 10 ** log_Teff_7[i_zams_7] / 1000],
         [LOGL_10[i_zams_10], LOGL_7[i_zams_7]], linestyle='-', color='gray')

# print len(teff_O), len(teff_B)


# plt.title('Fast-rotating O-type stars (>200 km/s)',fontsize=22)
# plt.title('Characterized sample of fast-rotating O-type stars',fontsize=22)

# plt.grid(color='black', linestyle='-', linewidth=0.1,alpha=0.1)
# plt.scatter(teff_gon_total,l_gon_total,color='grey',marker='o',linestyle='None',s=15,alpha=0.6)
# plt.scatter(teff_B,l_B,color='green',marker='o',label='B type SGc',linestyle='None',s=15,alpha=0.6)


# plt.scatter(teff_A,l_A,color='magenta',marker='o',label='A stars',linestyle='None',s=15)
# plt.text(54,4.37,' $\log{\mathcal{L}}=4\log{T_{eff}-\log{g}}$, $v_{rot}$=0 ($km$ $s^{-1}$)',fontsize=12)
# plt.text(54,3.64,' $\mathcal{L}=T_{eff}^{4}/g$',fontsize=16)


# plt.scatter(teff_gon_total,l_gon_total,color='black',marker='o',label='LPV (slow-rot)',linestyle='None',s=45,facecolor='black')
# plt.plot(teff_gon_total_sb1,l_gon_total_sb1,color='red',marker='o',linestyle='None',markersize=13)
# plt.plot(teff_gon_total_sb1p,l_gon_total_sb1p,color='green',marker='o',linestyle='None',markersize=13)
# plt.plot(teff_gon_total_sb1q,l_gon_total_sb1q,color='orange',marker='o',linestyle='None',markersize=13)


# plt.scatter(teff_work,l_work,color='black',marker='o',label='Fast-rotating stars',linestyle='None',s=45,facecolor='black')
# plt.plot(teff_work,l_work,color='black',marker='o',label='v sini > 200 $km$ $s^{-1}$',linestyle='None',markersize=10)
plt.plot(teff_work, l_real, color='grey', marker='o', label='Apparently single, B23+', linestyle='None', markersize=6)

# plt.scatter(teff_Ob,l_Ob,color='black',marker='o',linestyle='None',s=45,facecolor='red')
# for i in range(0,len(l_Ob)):
# plt.text(teff_Ob[i],l_Ob[i]+0.01,name_star_Ob[i],fontsize=7,rotation=30,color='red')#,bbox=dict(facecolor='white',edgecolor='None')


# plt.plot(teff_Ob_sb1,l_Ob_sb1,color='red',marker='o',linestyle='None',markersize=13)

plt.plot(teff_sb1, l_real_sb1, color='blue', marker='o', label='SB1, B23+', linestyle='None', markersize=13)

plt.plot(teff_sb1_2023, l_real_sb1_2023, color='red', marker='o', label='SB1, N23+, present work', linestyle='None',
         markersize=13)

# plt.plot(teff_O_SB1_original,l_O_SB1_original,color='red',marker='o',linestyle='None',markersize=10,label='SB1')


# plt.plot(teff_sb1q,l_sb1q,color='orange',marker='o',label='LPV/SB1?',linestyle='None',markersize=13)
# plt.plot(teff_sb1q,l_sb1q,color='orange',marker='o',linestyle='None',markersize=10,label='SB1?')

# plt.plot(teff_sb1p,l_sb1p,color='blue',marker='o',label='SB2/SB2?',linestyle='None',markersize=13)

# plt.plot(teff_gon_total_sb2p_fast,l_gon_total_sb2p_fast,color='red',marker='o',linestyle='None',markersize=10)


# plt.plot(teff_fast_rot/1000,l_Or_fast_rot,marker='o',label='Runaway stars',linestyle='None',markersize=20,markerfacecolor='none',color='red')

# plt.plot(teff_Ob_eb,l_Ob_eb,marker='X',linestyle='None',markersize=20,markerfacecolor='none',color='black')
plt.plot(teff_run, l_real_run, marker='o', linestyle='None', markersize=15, markerfacecolor='none', color='red',
         label='Runaways')
# plt.plot(teff_runq,l_runq,marker='o',linestyle='None',markersize=20,markerfacecolor='none',color='red')

# plt.plot(teff_gon_bh,l_gon_bh,color='blue',marker='x',linestyle='None',markersize=15,markerfacecolor='none')
# plt.text(30,4.15,'Cyg X-1', fontsize=16)

# plt.plot(teff_gon_bh1,l_gon_bh1,color='blue',marker='o',linestyle='None',markersize=15,markerfacecolor='none',label='Candidate BH')

# plt.plot(teff_gon_total_sb2p_fast,l_real_sb2,color='green',marker='o',label='LPV/SB2?',linestyle='None',markersize=10)

# plt.plot(teff_eb_fast,l_real_eb,marker='X',label='Eclipsing binary',linestyle='None',markersize=20,markerfacecolor='none',color='black')
# plt.plot(teff_el,l_real_ev,marker='s',label='Ellipsoidal variable',linestyle='None',markersize=16,markerfacecolor='none',color='black')


'''
####abundances + names
plt.scatter(teff_abun,l_real_abun,c=n_abun,s=700)
cbar=plt.colorbar()
cbar.set_label('log(N/H)+12 (Cazorla et al. 2017)')

for i in range(0,len(name_star)):

	plt.text(teff_work[i]+0.005, l_real[i]+0.0005,name_star[i],fontsize=5)


#####
'''
plt.legend(loc=3, fontsize=25)
plt.ylabel('log$_{10} (L/L_{\odot}$)', fontsize=30)
plt.xlabel('Teff [kK]', fontsize=30)
# plt.yticks(size=15)
# plt.xticks(size=15)
plt.gca().invert_xaxis()
# plt.xlim([51.5,17])
plt.xlim([51.5, 15])
plt.ylim([3, 6.0])
# plt.ylim([3.21,4.35])


plt.savefig(pp0, format='pdf')
pp0.close()
