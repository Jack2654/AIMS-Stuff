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


folder = "../../FHI-aims/French_NMR/Ethanol_new/Single/NMR_static/"
# pt.NMR_density(folder, atom_dict=[[4, 5], [6, 7, 8], [9]], color_dict=None, average=False, width=0.01)

folder = "../../FHI-aims/French_NMR/Ethanol_new/Double/NMR_static/"
# pt.NMR_density(folder, atom_dict=[[7, 8], [16, 17], [4, 5, 6], [13, 14, 15], [9, 18]],
#               color_dict=['r', 'r', 'g', 'g', 'b', 'b'], average=True, width=0.01)

folder = "../../FHI-aims/random/Shielding_tutorial/"
# all_dir = [os.fsdecode(x) for x in os.listdir(os.fsencode(folder))]
# all_dir.sort()
# for direct in all_dir:
#     continue
#     if os.path.isdir(folder + direct):
#         bao.create_shield_out(folder + direct + "/")

base = "../../FHI-aims/Double_Perovskites/New_Structures/Random/AgBi_static_Cs/"
readfile = "../../FHI-aims/Double_Perovskites/New_Structures/Random/AgBi/geometry.in"
disturbs = [0.05, 0.10, 0.15, 0.20]
for disturb in disturbs:
    write_base = base + str(int(disturb * 100)).rjust(2, "0") + "/"
    for i in range(1, 11):
        writefile = write_base + str(i).rjust(2, "0") + "/geometry.in"
        directory = write_base + str(i).rjust(2, "0") + "/"

        bg.disturb_positions(readfile, writefile, max_disturbance=disturb, change_Cs=False)
        # bg.recenter(writefile, writefile, 1)

folders = [base]
# pt.plot_random_structures(folders=folders, mode="dbmax")
# pt.plot_random_structures(folders=folders, mode="dbavg")
# pt.plot_random_structures(folders=folders, mode="ss")
# pt.plot_random_structures(folders=folders, mode="dd")
# pt.plot_random_structures(folders=folders, mode="ddisp")
# pt.plot_random_structures(folders=folders, mode="ddiagmax")
# pt.plot_random_structures(folders=folders, mode="ddiagAg")
# pt.plot_random_structures(folders=folders, mode="ddiagBi")
# pt.plot_random_structures(folders=folders, mode="ddiagavg")
# for i in range(1, 11):
#    Geometry.ddiag_val(base + f'20/{str(i).rjust(2, "0")}/geometry.in')


