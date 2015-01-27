# -*- coding: utf-8 -*-
"""
This script generates all relevant figures from the experiment and
store them in the `Figures` directory.

"""

from Modules.plotting import *

if __name__ == "__main__":
    """Choose functions to run here."""
    plt.close("all")
    if not os.path.isdir("Figures"):
        os.makedirs("Figures")
    
    plot_cp_curve(1.0, save=True, show=False)
    plot_perf_curves(save=True, show=False)
    plot_perf_re_dep(save=True, show=False, normalize_by="default", 
                     dual_xaxes=True)
    plot_wake_profiles(z_H=0.0, save=True)
    plot_meancomboquiv(U_infty=0.4, show=False)
    plot_meancomboquiv(U_infty=1.2, show=True)
#    plot_settling(1.0)