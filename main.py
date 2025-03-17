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
temp_dl = bg.Delta_L(base_path + f'{cur}_exp.in', bonds=bonds, flag="debug")
dl_avg.append(np.average(temp_dl))
dl_max.append(max(temp_dl))
l2.append(np.average(bg.Delta_L(base_path + f'{cur}_exp.in', mode="L2", bonds=bonds)))
betas.append(bg.random_angle_info(base_path + "06_exp.in", method="ignore", bonds=bonds))
betas_opt.append(bg.random_angle_info(base_path + f'{cur}_the.in', method="ignore", bonds=bonds))

for octahedra in octahedras:
    sigmas.append(Geometry.sig_sq_general(base_path + "06_exp.in", octahedra))
    delt_d.append(Geometry.delt_d_general(base_path + "06_exp.in", octahedra))

# for x in range(1, 19):
#     print(str(x) + ": " + bbo.band_info(base_path + "06_exp/", band=f'band10{str(x).rjust(2, "0")}.out', band_gap=False, spin_splitting=True, verbose=False))
pt.mulliken_plot(base_path + "06_exp/settings_final.in", quiet=False)

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

# note to use optimized structures, must change some of the bond index values

delta_beta = [max(beta) - min(beta) for beta in betas]
delta_beta[8] = 3.3256
# for beta in delta_beta:
#     print(beta)
slope, intercept, r_value, _, _ = scipy.stats.linregress(delta_beta, splittings)
plt.plot([min(delta_beta), max(delta_beta)], [x * slope + intercept for x in [min(delta_beta), max(delta_beta)]],
                 color='r', label=f'$r^2$={r_value ** 2}')
plt.scatter(delta_beta, splittings, color='r', label="_")
plt.xlabel("Delta Beta (degrees)")
plt.ylabel("Spin Splitting (eV)")
plt.title("Pb-Halide Experimental Structures")
plt.legend()
plt.show()

# code for 2D intensity maps
matplotlib.rc('text', usetex='true')
plt.scatter(delta_beta, l2, c=splittings, cmap='jet', edgecolors='k')
plt.colorbar(label='Spin Splitting (eV)')
plt.xlabel(r'$\Delta\beta$')
plt.ylabel(r'$L^2$')
plt.title("Pb Experimental")
plt.show()

read_folder = "../../FHI-aims/French_NMR/Ethanol_new/Double/MD_eth/geometries/geometry"
write_folder = "../../FHI-aims/French_NMR/Ethanol_new/Double/NMR_MD/"
# bf.make_MD_DMF(read_folder, write_folder, calcs=10, step=700, species=["H"])
#################### Info for 2x ethanol ####################
folder = "../../FHI-aims/French_NMR/Ethanol_new/Double/NMR_particular/"
atom_dict = [[[7, 8], [16, 17]], [[4, 5, 6], [13, 14, 15]], [[9, 18]]]
# atom_dict = [[[7], [8], [16], [17]], [[4], [5], [6], [13], [14], [15]], [[9], [18]]]
type = 'H'
xlim = [-4.25, 0.25]
# pt.NMR_average(folder, atom_dict, color_dict=['lightgreen', 'b', 'r'], width=0.005, xlim=xlim, type=type)


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

# plt.rcParams.update({
#     'axes.labelsize': 17,  # Axis title size
#     'xtick.labelsize': 15,  # X-axis tick label size
#     'ytick.labelsize': 15,  # Y-axis tick label size
#     'xtick.major.width': 5,  # X-axis major tick thickness
#     'ytick.major.width': 2,  # Y-axis major tick thickness
# })

width = 0.004
# pt.plot_chem_book(width=width)

folder = "../../FHI-aims/French_NMR/Ethanol_new/Single/NMR_static/"
# single all peaks
# pt.NMR_density(folder, atom_dict=[[4, 5], [6, 7, 8], [9]], color_dict=['r', 'g', 'b'], average=False, width=width)
# single with averaging
# pt.NMR_density(folder, atom_dict=[[4, 5], [6, 7, 8], [9]], color_dict=['r', 'g', 'b'], average=True, width=width)



folder = "../../FHI-aims/French_NMR/Ethanol_new/Double/NMR_particular/"
# double all peaks MD
# pt.NMR_density(folder, atom_dict=[[7, 8], [16, 17], [4, 5, 6], [13, 14, 15], [9, 18]],
#                color_dict=['r', 'r', 'g', 'g', 'b', 'b', 'b'], average=False, MD=True, width=width)
# double MD + averaged geometry
# pt.NMR_density(folder, atom_dict=[[7, 8, 16, 17], [4, 5, 6, 13, 14, 15], [9, 18]],
#                color_dict=['r', 'g', 'b'], average=True, MD=True, width=width)

folder = "../../FHI-aims/French_NMR/Ethanol_new/Double/NMR_static/"
# pt.NMR_density(folder, atom_dict=[[7, 8], [16, 17], [4, 5, 6], [13, 14, 15], [9], [18]],
#                color_dict=['r', 'r', 'g', 'g', 'b', 'b', 'b'], average=False, width=0.001)

folder = "../../FHI-aims/random/Shielding_tutorial/"
# all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
# all_dir.sort()
# for direct in all_dir:
#     continue
#     if os.path.isdir(folder + direct):
#         bao.create_shield_out(folder + direct + "/")

base = "../../FHI-aims/Double_Perovskites/New_Structures/Random/Pb_150/"
readfile = "../../FHI-aims/Double_Perovskites/New_Structures/Random/Pb_150/geometry.in"
disturbs = [0.05, 0.10, 0.15, 0.20, 0.25]
for disturb in disturbs:
    continue
    write_base = base + str(int(disturb * 100)).rjust(2, "0") + "/"
    for i in range(1, 11):
        writefile = write_base + str(i).rjust(2, "0") + "/geometry.in"
        directory = write_base + str(i).rjust(2, "0") + "/"
        write2 = write_base + str(i).rjust(2, "0") + "/flat.in"
        bg.flatten(writefile, write2)
        # bg.disturb_positions(readfile, writefile, max_disturbance=disturb, change_Cs=False)
        # bg.recenter(writefile, writefile, 1)

base_path = "../../FHI-aims/Double_Perovskites/New_Structures/Random/"
folders = ["AgBi", "AgBi_recenter", "AgBi_static_Cs"]
modes = ["db_avg", "db_max", "ss_avg", "ss_max", "dd_avg", "dd_max",
         "disp_diag_ag", "disp_diag_bi", "disp_diag_avg", "disp_diag_max", "disp_avg", "disp_max"]

# for mode in modes:
# for folder in folders:
#     pt.plot_random_structures([base_path + folder + "/"], folder)
pt.plot_random_structures([f'{base_path}Pb_150/'], title="Pb")
# pt.plot_random_structures([f'{base_path}AgBi_static_Cs/'], title="AgBi")


base = "../../FHI-aims/Double_Perovskites/New_Structures/Random/Pb_150/"
set = "20"
# for i in range(8, 11):
for i in [4]:
    continue
    bands = []
    dbs = []
    bands.append(bbo.band_info(base + f'{set}/{str(i).rjust(2, "0")}/', "band1001.out", band_gap=False, spin_splitting=True, verbose=False))
    bands.append(bbo.band_info(base + f'{set}/{str(i).rjust(2, "0")}/', "band1002.out", band_gap=False, spin_splitting=True, verbose=False))
    bands.append(bbo.band_info(base + f'{set}/{str(i).rjust(2, "0")}/', "band1003.out", band_gap=False, spin_splitting=True, verbose=False))
    bands.append(bbo.band_info(base + f'{set}/{str(i).rjust(2, "0")}/', "band1004.out", band_gap=False, spin_splitting=True, verbose=False))
    bands.append(bbo.band_info(base + f'{set}/{str(i).rjust(2, "0")}/', "band1005.out", band_gap=False, spin_splitting=True, verbose=False))
    angles = bg.random_angle_info(base + f'{set}/{str(i).rjust(2, "0")}/geometry.in', method="flat")
    print(angles)
    dbs.append(abs(np.average([angles[0], angles[1]]) - np.average([angles[2], angles[3]])))
    dbs.append(abs(np.average([angles[0], angles[2]]) - np.average([angles[1], angles[3]])))
    print(" ".join([str(round(float(x), 4)) for x in bands]))
    print(" ".join([str(round(x, 4)) for x in dbs]))
    # print(max(dbs))
    pt.mulliken_plot(base + f'{set}/{str(i).rjust(2, "0")}/settings_final.in', quiet=True)


# BZ plotting code for 2D
corner_file = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/2D_experimental/BZ_corners.in"
corners = np.asarray(bf.read_BZ_corners(corner_file))
adjacency = [[1, 6, 8], [7, 9], [3, 4, 10], [5, 11], [6, 5], [7], [7],
             [], [9, 10], [11], [11], []]
pathway = [[0, 0.5, 0, 0, 0, 0],
           [0, 0, 0, -0.1405892226, 0.0000000000, 0.5132880253],
           [-0.1231599653, 0.5000000000, 0.4496542069, 0, 0, 0],
           [0, 0, 0, -0.1231599653, -0.5000000000, 0.4496542069],
           [0, 0, 0, 0.5132880253, 0.0000000000, -0.0485142755]]

names = [['1', 'Γ'], ['Γ', '2'], ['3', 'Γ'], ['Γ', '4'], ['Γ', '5']]
geo_file = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/2D_experimental/geometry.in"
# pt.plot_3d_solid_with_path_and_names(geo_file, corners, adjacency, pathway, names, save=False, setup=False)

dos_folder = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/2D_experimental/"
file = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/2D_experimental/dos.png"
# pt.dos_plot(dos_folder, shift=0.9, limits=[-1, 4], save=True, title="Experimental pDOS", filename=file)

dos_folder = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/2D_theoretical/"
file = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/2D_theoretical/dos.png"
# pt.dos_plot(dos_folder, shift=0.9, limits=[-1, 4], save=True, title="Theoretical pDOS", filename=file)

base = "../../FHI-aims/Yi/Yi_1_5_D/band_plotting_folder/"
structures = ["2D_", "n_2_4_", "n_3_", "n_4_", "n_4_6_", "n_5_"]
structures = ["n_5_"]
types = "experimental", "theoretical"
for structure in structures:
    continue
    for type in types:
        folder = f'{base}{structure}{type}/'
        # print(folder)
        # pt.mulliken_plot(f'{folder}settings_final.in', save=False)
        with open(f'{folder}settings_final.in', "r") as f:
            for line in f.readlines():
                if "bands" in line:
                    bands = line.strip().split()[1:]
        # for band in bands:
        #     print(bbo.band_info(folder, f'band{band}.out', conduction_offset=0, band_gap=False, spin_splitting=True, verbose=False))
        # print(bands)
    #




