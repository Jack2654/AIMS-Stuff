import BasicFunc
import Geometry
import OneTimeScripts as ots
import BasicGeo as bg
import PlottingTools as pt
import BasicBandOut as bbo
import BasicAimsOut as bao
import BasicControl as bc
import BasicFunc as bf
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import Methods2024 as M24
import scipy

splittings = [0.0440, 0.1080, 0.1370, 0.0000, 0.0820, 0.1400, 0.0330, 0.0090, 0.0120, 0.0000, 0.0920]
sigmas = []
delt_d = []
betas = []
betas_opt = []
dl_avg = []
dl_max = []
l2 = []
L_diffs = []
base_path = "../../FHI-aims/Double_Perovskites/New_Structures/Experimental/"

# 01  =  [R-4-Cl-MBA]2 PbBr4
cur = "01"
bonds = [[[1], [9], [3]],
         [[3], [11], [1, "001"]],
         [[1, "001"], [7, "010"], [3, "010"]],
         [[3, "010"], [5], [1]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 02  =  [S-1-1-NEA]2 PbBr4
cur = "02"
bonds = [[[1, "001"], [7], [2]],
         [[2], [9], [1, "011"]],
         [[1, "011"], [10, "011"], [2, "001"]],
         [[2, "001"], [8, "001"], [1, "001"]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 03  =  [S-4-NO2-MBA]2 PbBr4 H20
cur = "03"
bonds = [[[2, "001"], [5, "001"], [1]],
         [[1], [9, "001"], [2, "011"]],
         [[2, "011"], [10, "011"], [1, "001"]],
         [[1, "001"], [6, "001"], [2, "001"]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 04  =  [S-2-Me-BuA]2 PbBr4
cur = "04"
bonds = [[[1, "001"], [3, "001"], [2]],
         [[2], [4, "010"], [1, "011"]],
         [[1, "011"], [5, "001"], [2, "001"]],
         [[2, "001"], [6, "001"], [1, "001"]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 05  =  [NMA]2 PbBr4
cur = "05"
bonds = [[[2], [5], [1]],
         [[1], [17], [2, "010"]],
         [[2, "010"], [18, "011"], [1, "001"]],
         [[1, "001"], [6], [2]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 06  =  [4AMP] PbBr4
cur = "06"
bonds = [[[4], [19], [3]],
         [[3], [15], [2]],
         [[2], [14], [3, "001"]],
         [[3, "001"], [20], [4]]]
octahedras = [[[4], [8], [12], [13], [16], [20], [19]],
              [[2], [6], [10], [14], [15], [17], [18, [0, 1, 0]]],
              [[3, [0, 0, 1]], [7, [0, 0, 1]], [11, [0, 0, 1]], [15, [0, 0, 1]], [19, [0, 0, 1]], [14], [20]],
              [[1, [0, 0, 1]], [5, [0, 0, 1]], [9, [0, 0, 1]], [13, [0, 0, 1]], [16], [18], [17, [0, -1, 1]]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# for octahedra in octahedras:
#     sigmas.append(Geometry.sig_sq_general(base_path + "06_exp.in", octahedra))
#     delt_d.append(Geometry.delt_d_general(base_path + "06_exp.in", octahedra))
# for x in range(1, 19):
#     print(str(x) + ": " + bbo.band_info(base_path + "06_exp/", band=f'band10{str(x).rjust(2, "0")}.out', band_gap=False, spin_splitting=True, verbose=False))
# pt.mulliken_plot(base_path + "06_exp/settings_final.in", quiet=False)

# 07  =  [S-MBA]2 PbI4
cur = "07"
bonds = [[[4, "010"], [11, "010"], [3]],
         [[3], [10, "011"], [4, "011"]],
         [[4, "011"], [6, "011"], [3, "010"]],
         [[3, "010"], [7, "010"], [4, "010"]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 08  =  [S-4-NH3-MBA] PbI4
cur = "08"
bonds = [[[3, "001"], [7, "001"], [1]],
         [[1], [11, "011"], [3, "011"]],
         [[3, "011"], [9, "001"], [1, "001"]],
         [[1, "001"], [5, "001"], [3, "001"]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 09  =  [R-4-Cl-MBA]2 PbI4
cur = "09"
bonds = [[[2], [6], [1]],
         [[1], [4, "010"], [2, "010"]],
         [[2, "010"], [3, "010"], [1, "001"]],
         [[1, "001"], [5], [2]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 10  =  [S-MHA]2 PbI4
cur = "10"
bonds = [[[1, "001"], [7, "001"], [3]],
         [[3], [13, "001"], [1, "011"]],
         [[1, "011"], [15, "011"], [3, "001"]],
         [[3, "001"], [5, "001"], [1, "001"]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

# 11  =  [PMA]2 PbCl4
cur = "11"
bonds = [[[4], [8], [3]],
         [[3], [20], [4, "010"]],
         [[4, "010"], [19, "010"], [3, "001"]],
         [[3, "001"], [7], [4]]]
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds)
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + f'{cur}_exp.in', method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))
L_diffs.append(bg.L_diff(base_path + f'{cur}_exp.in', bonds=bonds, mode="all"))

print(L_diffs)
print(np.average(L_diffs))
print("now to AgBi:")

# note to use optimized structures, must change some of the bond index values
# delta_beta = [max(beta) - min(beta) for beta in betas]
# delta_beta[8] = 3.3256
# for beta in delta_beta:
#     print(beta)
# slope, intercept, r_value, _, _ = scipy.stats.linregress(delta_beta, splittings)
# plt.plot([min(delta_beta), max(delta_beta)], [x * slope + intercept for x in [min(delta_beta), max(delta_beta)]],
#                  color='r', label=f'$r^2$={r_value ** 2}')
# plt.scatter(delta_beta, splittings, color='r', label="_")
# plt.xlabel("Delta Beta (degrees)")
# plt.ylabel("Spin Splitting (eV)")
# plt.title("Pb-Halide Experimental Structures")
# plt.legend()
# plt.show()

# code for 2D intensity maps
# matplotlib.rc('text', usetex='true')
# plt.scatter(delta_beta, l2, c=splittings, cmap='jet', edgecolors='k')
# plt.colorbar(label='Spin Splitting (eV)')
# plt.xlabel(r'$\Delta\beta$')
# plt.ylabel(r'$L^2$')
# plt.title("Pb Experimental")
# plt.show()

# for bonds make all start from the same place in Ag and radiate out to Bi
splittings = [0.0440, 0.1080, 0.1370, 0.0000, 0.0820, 0.1400, 0.0330, 0.0090, 0.0120, 0.0000, 0.0920]
sigmas = []
delt_d = []
betas = []
betas_opt = []
dl_avg = []
dl_max = []
l2 = []
base_path = "../../FHI-aims/Double_Perovskites/Real_Structures/real-systems/"

# 01  =  ijayuq
cur = "ijayuq"
bonds = [[[19], [10], [2]],
         [[19, "001"], [3, "001"], [1]],
         [[19, "101"], [6], [2, "001"]],
         [[19, "101"], [11], [1]]]
print(bg.L_diff(base_path + cur + "/geometry.in", bonds=bonds, frac=True))

# 02  =  selwuf
cur = "experimental-bs-selwuf"
bonds = [[[2], [9], [1]],
         [[2, "010"], [8], [1]],
         [[2, "110"], [6], [1]],
         [[2, "100"], [3], [1]]]
print(bg.L_diff(base_path + cur + "/geometry.in", bonds=bonds, frac=True))

# 03  =  selwoz
cur = "experimental-bs-selwoz"
bonds = [[[86], [2], [1]],
         [[86], [4, "010"], [1, "010"]],
         [[86], [7, "110"], [1, "110"]],
         [[86], [3, "100"], [1, "100"]]]
print(bg.L_diff(base_path + cur + "/geometry.in", bonds=bonds, frac=True))

types = ["AgBi", "Pb"]
for type in types:
    continue
    folder = f'../../FHI-aims/Double_Perovskites/New_Structures/Delta_Beta/{type}/bands/'
    pt.plot_correlations(folder, mode="ss", band_path=3, title=r'Band Path: $M\rightarrow\Gamma$')
    pt.plot_correlations(folder, mode="ss", band_path=4, title=r'Band Path: $\Gamma\rightarrow P$')

folders = ["../../FHI-aims/Double_Perovskites/New_Structures/Maurer/AgBi/"]
folders += ["../../FHI-aims/Double_Perovskites/New_Structures/Maurer/TlBi/"]
folders += ["../../FHI-aims/Double_Perovskites/New_Structures/Maurer/Pb/"]

# pt.plot_maurer_structures(folders, path=3, mode='ddiag')
# pt.plot_maurer_structures(folders, path=3, mode='db')
# pt.plot_maurer_structures(folders, path=3, mode='dd')
# pt.plot_maurer_structures(folders, path=3, mode='ss')
#

base_path = "../../FHI-aims/Double_Perovskites/New_Structures/Random/"
folders = ["AgBi", "AgBi_recenter", "AgBi_static_Cs"]

# for folder in folders:
#     pt.plot_random_structures([base_path + folder + "/"], folder)
# pt.plot_random_structures([f'{base_path}Pb_150/'], title="Pb")
# pt.plot_random_structures([f'{base_path}AgBi_static_Cs/'], title="AgBi")

folder = "../../FHI-aims/Bi_substitution/Output/"
pt.dos_plot(folder, shift=0, limits=[-3, 3], save=False, title="pDOS", filename=folder + "dos.png")

