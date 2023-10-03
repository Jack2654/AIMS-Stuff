import OneTimeScripts as ots
import BasicGeo as bg
import PlottingTools as pt
import BasicBandOut as bbo

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

for structure in range(8):
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
    # pt.mulliken_plot(file, figure_loc + ".ne.png", eshift, ymin, ymax, substate, cd, labels, title, eq=False, debug=False)

folder = "../../FHI-aims/Yi_1_5_D/Results/New_Results/"
options = ["n_2_4/", "n_3/", "n_4/", "n_5/"]
options = ["n_3/"]
for opt in options:
    continue
    base = folder + opt
    exp = base + "experimental/"
    the = base + "theoretical/"
    print(base)
    for i in range(1, 7):
        band = "band100" + str(i) + ".out"
        if i == 6 and "3" not in opt and "5" not in opt:
            continue
        dbg = False

        bbo.band_info(exp, band, steps=1, band_gap=True, k0=False, spin_splitting=True,
                      effective_mass=False, verbose=False, debug=dbg)
        bbo.band_info(the, band, steps=1, band_gap=True, k0=False, spin_splitting=True,
                      effective_mass=False, verbose=False, debug=dbg)

        continue
        temp = []
        for j in range(1, 6):
            temp.append(bbo.band_info(exp, band, steps=j, band_gap=False, k0=False, spin_splitting=False,
                                      effective_mass=True, verbose=False, debug=dbg))
        temp_cpy = temp.copy()
        temp.sort()
        print(str(temp[2]) + " " + str(temp_cpy))

        temp = []
        for j in range(1, 6):
            temp.append(bbo.band_info(the, band, steps=j, band_gap=False, k0=False, spin_splitting=False,
                                      effective_mass=True, verbose=False, debug=dbg))
        temp_cpy = temp.copy()
        temp.sort()
        print(str(temp[2]) + " " + str(temp_cpy))

folder = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/python_dist_positions/"
disturbs = ["5/", "10/", "15/", "20/"]
for dist in disturbs:
    temp = folder + dist
    for i in range(10):
        current = temp + str(i) + "/"
        for j in range(1, 7):
            band = "band100" + str(j) + ".out"
            bbo.band_info(current, band, band_gap=True, spin_splitting=True, verbose=False)
