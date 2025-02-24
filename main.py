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


base = "../../FHI-aims/Double_Perovskites/New_Structures/Random/Pb_150/"
set = "15"
# for i in range(8, 11):
for i in [5]:
    # continue
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




