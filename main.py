import OneTimeScripts as ots
import BasicGeo as bg
import PlottingTools as pt

folder = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/setup/"
disturbances = [0.05, 0.10, 0.15, 0.20]
base = folder + "geometry_base.in"
new_folder = folder[:-6] + "structures/"
for disturb in disturbances:
    for i in range(10):
        read_file = new_folder + str(int(disturb * 100)) + "_" + str(i) + ".in"
        # bg.angle_info(read_file, 2, 17, 4)
        # bg.angle_info(read_file, 2, 14, 3)
        # bg.angle_info(read_file, 3, 15, 1)
        # bg.angle_info(read_file, 1, 13, 4)

folder = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/pyth_dist_positions/"
for disturb in disturbances:
    base = folder + str(int(disturb * 100)) + "/"
    for i in range(10):
        calculation = base + str(i) + "/"
        # print(base)

for structure in range(7, 8):
    # n_2_4 experimental:
    if structure == 0:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_2_4/experimental/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_2_4_experimental.png"
        eshift = -2.2702
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}
        labels = ('$1$', '$\Gamma$', '$2\|3$', '$\Gamma$', "$4'\|\Gamma$", '$5$')
        title = "Experimental m=4, n=2"
        equal = True
        debug = False
    # n_2_4 theoretical:
    if structure == 1:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_2_4/theoretical/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_2_4_theoretical.png"
        eshift = -2.2211
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}
        labels = ('$1$', '$\Gamma$', '$2\|3$', '$\Gamma$', "$4'\|\Gamma$", '$5$')
        title = "Theoretical m=4, n=2"
        equal = True
        debug = False

    # n_3 experimental:
    if structure == 2:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_3/experimental/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_3_experimental.png"
        eshift = -0.3635
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}
        labels = ('$1$', '$\Gamma$', '$2\|\Gamma$', '$3\|\Gamma$', '$4\|\Gamma$', "$5'\|\Gamma$", '$6$')
        title = "Experimental m=3, n=3"
        equal = True
        debug = False
    # n_3 theoretical:
    if structure == 3:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_3/theoretical/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_3_theoretical.png"
        eshift = -1.3568
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}
        labels = ('$1$', '$\Gamma$', '$2\|\Gamma$', '$3\|\Gamma$', '$4\|\Gamma$', "$5'\|\Gamma$", '$6$')
        title = "Theoretical m=3, n=3"
        equal = True
        debug = False

    # n_4 experimental:
    if structure == 4:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_4/experimental/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_4_experimental.png"
        eshift = -1.2848
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "F": "k"}
        labels = ('$1$', '$\Gamma$', "$2'\|\Gamma$", '$3\|\Gamma$', '$4\|\Gamma$', '$5$')
        title = "Experimental m=4, n=4"
        equal = True
        debug = False
    # n_4 theoretical:
    if structure == 5:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_4/theoretical/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_4_theoretical.png"
        eshift = -1.3466
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "F": "k"}
        labels = ('$1$', '$\Gamma$', "$2'\|\Gamma$", '$3\|\Gamma$', '$4\|\Gamma$', '$5$')
        title = "Theoretical m=4, n=4"
        equal = True
        debug = False

    # n_5 experimental:
    if structure == 6:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_5/experimental/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_5_experimental.png"
        eshift = -2.0580
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k", "O": "r"}
        labels = ('$1$', '$\Gamma$', '$2\|\Gamma$', '$3\|\Gamma$', '$4\|\Gamma$', '$5\|\Gamma$', "$6'$")
        title = "Experimental m=5, n=5"
        equal = True
        debug = False
    # n_5 theoretical:
    if structure == 7:
        file = "../../FHI-aims/Yi_1_5_D/Results/New_Results/n_5/theoretical/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_5_theoretical.png"
        eshift = -2.2826
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k", "O": "r"}
        labels = ('$1$', '$\Gamma$', '$2\|\Gamma$', '$3\|\Gamma$', '$4\|\Gamma$', '$5\|\Gamma$', "$6'$")
        title = "Theoretical m=5, n=5"
        equal = True
        debug = False
    pt.mulliken_plot(file, figure_loc, eshift, ymin, ymax, substate, cd, labels, title, equal, debug)
