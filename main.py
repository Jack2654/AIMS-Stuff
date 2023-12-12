import OneTimeScripts as ots
import BasicGeo as bg
import PlottingTools as pt
import BasicBandOut as bbo
import BasicControl as bc
import matplotlib.pyplot as plt
import math

settings = "../../FHI-aims/Double_Perovskites/Figures/Cs_test_OOP/"
options = ["base/", "disp_3/", "disp_4/"]
for opt in options:
    temp = f'%s%ssettings.in' % (settings, opt)
    # pt.mulliken_plot(temp, debug=False, save=True)

settings = "../../FHI-aims/Double_Perovskites/Figures/plotting/Random_plotting_20_03/settings.in"
pt.mulliken_plot(settings, debug=False, save=True)


