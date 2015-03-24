# -*- coding: utf-8 -*-
"""
This script generates all relevant figures from the experiment and
stores them in the `Figures` directory.

"""

from Modules.plotting import *

show = False
save = True
savetype = ".pdf"

def main():
    plot_cp_curve(1.0, save=save, savetype=savetype)
    plot_perf_curves(save=save, savetype=savetype)
    plot_perf_re_dep(save=save, savetype=savetype, normalize_by="default", 
                     dual_xaxes=True)
    plot_wake_profiles(z_H=0.0, save=save, savetype=savetype)
    plot_meancontquiv(U_infty=0.4, save=save, savetype=savetype)
    plot_meancontquiv(U_infty=1.2, save=save, savetype=savetype)
    if show:
        plt.show()

if __name__ == "__main__":
    if not os.path.isdir("Figures"):
        os.makedirs("Figures")
    main()
    