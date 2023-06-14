import PlottingTools as pt
import Geometry
import ASE_tools as ASE_t
import Basic

fp = "../../FHI-aims/Yi_1_5_D/Results/Band_Results/n_3/theoretical"
color_dict = {"Pb": "m", "I": "g", "N": "b", "C": "y", "H": "c", "S": "k"}

# pt.mulliken_plot(filepath=fp, energyshift=0, ymin=-2.1, ymax=2.1, color_dict=color_dict, debug=True)
# pt.mulliken_plot(filepath=start2, energyshift=-2.10031361, ymin=-4.0, ymax=4.0, color_dict=color_dict)

file1 = "../../FHI-aims/KellerPBE/s66x10_in_copy/"
file2 = "../../FHI-aims/KellerPBE/all_geometries/xx_1.00.in"
file3 = "../../FHI-aims/KellerPBE/all_geometries/xx_1.10.in"

number = "07"
a = 2
b = 12

file2 = file2.replace("xx", number)
file3 = file3.replace("xx", number)
one = Basic.distance(file2, a, file2, b)
two = Basic.distance(file3, a, file3, b)
print(two/one)