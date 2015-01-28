# -*- coding: utf-8 -*-
"""
This script generates all relevant figures from the experiment and
stores them in the `Figures` directory.

"""

from Modules.plotting import *

if __name__ == "__main__":
    if not os.path.isdir("Figures"):
        os.makedirs("Figures")

    plt.close("all")    
    plot_cp_curve(1.0, save=True, show=False)
    plot_perf_curves(save=True, show=False)
    plot_perf_re_dep(save=True, show=False, normalize_by="default", 
                     dual_xaxes=True)
    plot_wake_profiles(z_H=0.0, save=True)
    plot_meancontquiv(U_infty=0.4, show=False)
    plot_meancontquiv(U_infty=1.2, show=True)