import PlottingTools as pt
import Geometry

b = 150
a = 170

# aang = 150
# newfile = "../../FHI-aims/TlBi-Perovskites/in_plane/"
# newfile += str(aang) + "/" + str(aang) + "_"
# for x in range(9):
#   fill = newfile + str(aang - 20 + x * 5) + "/geometry.in"
#   Geometry.centroid(fill, -10, 4, -10, 4, "Tl", False)

geo = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/Relax/light11/geometry.in.next_step"
geo_two = "../../FHI-aims/ChiralSurfaceDefects/S_AEA/S_AEA_slab.in"
source = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/AngleAnalysis/sourceA2.in"
source2 = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/AngleAnalysis/sourceT2.in"
write = "../../FHI-aims/ChiralSurfaceDefects/SurfaceDefects/S_AEA/Relax/light11/geometry.in.moved"
geo0 = "../../FHI-aims/AgBi-Perovskites/ideal/out/final_OOP/180/180_200/band1001.out"
geo1 = "../../FHI-aims/AgBi-Perovskites/ideal/out/bro_idk/180/180_200/band1001.out"

# m23.manyLocAngles(source, write, 1, "s", True)
# print()
# m23.manyAbsAngles(source, write, 1, "s", True)
# m23.moveManyPoints(source, geo, 0, write, True)
# m23.moveManyPoints(source2, write, 0, write, False)
# m23.moveManyPoints(source2, write, 1, write, True)
# mecha.Delt_Y_Transform("D", 12, 6, 14)

start = "../../FHI-aims/CationDoping/New_Results/RR_Zoom"
start2 = "../../FHI-aims/CationDoping/New_Results/DR_Zoom"
fp = "../../FHI-aims/Yi_1_5_D/Results/Band_Results/n_3/theoretical"
color_dict = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}

pt.mulliken_plot(filepath=fp, energyshift=0, ymin=-2.1, ymax=2.1, color_dict=color_dict, debug=True)
# pt.mulliken_plot(filepath=start2, energyshift=-2.10031361, ymin=-4.0, ymax=4.0, color_dict=color_dict)

# file = "../../FHI-aims/CationDoping/Final_Results/Regular_Relaxed/Zoom/band100"
# for i in range(4):
#    file_plus = file + str(i + 1) + ".out"
    #Geometry.bandProcess(file_plus, True, 0.0)

#file = "../../FHI-aims/AgBi-Perovskites/ideal/out/final_OOP/180/180_"
#for i in range(9):
#    file_plus = file + str(160 + i * 5) + "/band1002.out"
#    #print(file_plus)
#    Geometry.bandProcess(file_plus, True, 0.0)


# file1 = "../../FHI-aims/Yi_1_5_D/n_3_3/no_cation.in"
# file2 = "../../FHI-aims/Yi_1_5_D/n_3_3/cation_two.in"
# location = "../../FHI-aims/Yi_1_5_D/n_3_3/config_two.in"
# Geometry.combine(file1, file2, location, 13, 2)