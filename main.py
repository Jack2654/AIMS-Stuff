import OneTimeScripts as ots
import BasicGeo as bg
import PlottingTools as pt
import BasicBandOut as bbo
import BasicControl as bc

folder = "../../FHI-aims/Double_Perovskites/AgBi-Perovskites/ideal/disturbed_positions/python_dist_positions/"
disturbs = ["5/", "10/", "15/", "20/"]

if True:
    corners_ag1 = [20, 17, 14, 18, 26, 22]
    shiftmap_ag1 = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    center_ag1 = 2

    corners_ag2 = [13, 15, 16, 19, 21, 25]
    corners_sig_ag2 = [13, 16, 19, 15, 25, 21]
    shiftmap_ag2 = [[0, 0, 0], [0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, 0], [0, 0, 0]]
    center_ag2 = 1

    corners_bi1 = [14, 15, 16, 20, 23, 27]
    corners_sig_bi1 = [14, 15, 20, 16, 27, 23]
    shiftmap_bi1 = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [1, 0, 0], [0, 0, 0], [0, 0, 0]]
    center_bi1 = 3

    corners_bi2 = [13, 17, 19, 18, 24, 28]
    corners_sig_bi2 = [19, 18, 13, 17, 28, 24]
    shiftmap_bi2 = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 0]]
    center_bi2 = 4

    beta_1 = [2, 14, 3]
    shiftmap_1 = -1
    beta_2 = [3, 20, 2]
    shiftmap_2 = [[0, 0, 0], [1, 0, 0], [1, 0, 0]]
    beta_3 = [2, 17, 4]
    shiftmap_3 = -1
    beta_4 = [3, 15, 1]
    shiftmap_4 = -1
    beta_5 = [4, 13, 1]
    shiftmap_5 = -1
    beta_6 = [1, 19, 4]
    shiftmap_6 = [[0, 0, 0], [1, 0, 0], [1, 0, 0]]
    beta_7 = [4, 18, 2]
    shiftmap_7 = [[0, 0, 0], [0, 1, 0], [0, 1, 0]]
    beta_8 = [1, 16, 3]
    shiftmap_8 = [[0, 0, 0], [0, 1, 0], [0, 1, 0]]

for dist in disturbs:
    # print(dist)
    temp = folder + dist
    for i in range(10):
        current = temp + str(i) + "/"
        for j in range(1, 2):
            band = "band100" + str(j) + ".out"
            # res = bg.angle_info(current + "geometry.in", beta_8, shiftmap_8)
            # print(f'%.10f' % (res[0]))
            # print(f'%0.10f' % bg.delta_d(current + "geometry.in", center_bi2, corners_sig_bi2, debug=False))
            # print(f'%0.10f' % bg.sigma_squared(current + "geometry.in", center_bi2, corners_sig_bi2, debug=False))
            # bbo.band_info(current, band, band_gap=True, spin_splitting=True, verbose=False)

base = "../../FHI-aims/"
control = base + "Yi_1_5_D/Orbital_Contributions/n_2_4_experimental/control_tight.in"
geo = base + "Yi_1_5_D/Orbital_Contributions/n_2_4_experimental/geometry.in"
species = base + "Repositories/Repo_mins/species_defaults/defaults_2020/tight/"
# bc.write_species(control, geo, species)

for structure in range(0, 2):
    # n_2_4 experimental:
    if structure == -1:
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
    if structure == -1:
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
    if structure == 0:
        file = "../../FHI-aims/Yi_1_5_D/Diagonal_IP_Calcs/n_3/experimental/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_3_diag_experimental.png"
        eshift = -0.3635
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}
        labels = ('$X$', '$\Gamma$', '$Y$')
        title = "Experimental m=3, n=3"
        equal = True
        debug = False
    # n_3 theoretical:
    if structure == 1:
        file = "../../FHI-aims/Yi_1_5_D/Diagonal_IP_Calcs/n_3/theoretical/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_3_diag_theoretical.png"
        eshift = -1.3568
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}
        labels = ('$X$', '$\Gamma$', '$Y$')
        title = "Theoretical m=3, n=3"
        equal = True
        debug = False

    # n_4 experimental:
    if structure == 2:
        file = "../../FHI-aims/Yi_1_5_D/Diagonal_IP_Calcs/n_4/experimental/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_4_diag_experimental.png"
        eshift = -1.2848
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "F": "k"}
        labels = ('$X$', '$\Gamma$', '$Y$')
        title = "Experimental m=4, n=4"
        equal = True
        debug = False
    # n_4 theoretical:
    if structure == 3:
        file = "../../FHI-aims/Yi_1_5_D/Diagonal_IP_Calcs/n_4/theoretical/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_4_diag_theoretical.png"
        eshift = -1.3466
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "F": "k"}
        labels = ('$X$', '$\Gamma$', '$Y$')
        title = "Theoretical m=4, n=4"
        equal = True
        debug = False

    # n_5 experimental:
    if structure == 4:
        file = "../../FHI-aims/Yi_1_5_D/Diagonal_IP_Calcs/n_5/experimental/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_5_diag_experimental.png"
        eshift = -2.0580
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k", "O": "r"}
        labels = ('$X$', '$\Gamma$', '$Y$')
        title = "Experimental m=5, n=5"
        equal = True
        debug = False
    # n_5 theoretical:
    if structure == 5:
        file = "../../FHI-aims/Yi_1_5_D/Diagonal_IP_Calcs/n_5/theoretical/"
        figure_loc = "../../FHI-aims/Yi_1_5_D/Results/Figures/New_Figures/n_5_diag_theoretical.png"
        eshift = -2.2826
        ymin = -1
        ymax = 4
        substate = 0
        cd = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k", "O": "r"}
        labels = ('$X$', '$\Gamma$', '$Y$')
        title = "Theoretical m=5, n=5"
        equal = True
        debug = False
    print("Band Gap, Spin Splitting, Effective Mass for " + file + ":")
    for j in range(1, 3):
        print()
        bbo.band_info(file, steps=1, band="band100" + str(j) + ".out", band_gap=True, spin_splitting=True,
                      effective_mass=True, verbose=False)
        bbo.band_info(file, steps=2, band="band100" + str(j) + ".out", band_gap=True, spin_splitting=True,
                      effective_mass=True, verbose=False)
        bbo.band_info(file, steps=3, band="band100" + str(j) + ".out", band_gap=True, spin_splitting=True,
                      effective_mass=True, verbose=False)
        bbo.band_info(file, steps=4, band="band100" + str(j) + ".out", band_gap=True, spin_splitting=True,
                      effective_mass=True, verbose=False)
        bbo.band_info(file, steps=5, band="band100" + str(j) + ".out", band_gap=True, spin_splitting=True,
                      effective_mass=True, verbose=False)
    # pt.mulliken_plot(file, figure_loc, eshift, ymin, ymax, substate, cd, labels, title, eq=equal, debug=False)
