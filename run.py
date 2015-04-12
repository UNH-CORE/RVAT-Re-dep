# -*- coding: utf-8 -*-
"""
This script generates all relevant figures from the experiment and
stores them in the `Figures` directory.

"""

from Modules.plotting import *

set_sns()
show = True
save = True
savetype = ".pdf"

def main():
    plot_perf_curves(save=save, savetype=savetype)
    plot_perf_re_dep(save=save, savetype=savetype)
    plot_perf_re_dep(subplots=False, save=save, savetype=savetype)
    plot_wake_profiles(save=save, savetype=savetype)
    make_k_bar_graph(save=save, savetype=savetype)
    make_mom_bar_graph(save=save, savetype=savetype)
    wm = WakeMap(0.4)
    wm.plot_meancontquiv(save=save, savetype=savetype)
    wm.plot_k(save=save, savetype=savetype)
    wm = WakeMap(1.2)
    wm.plot_meancontquiv(save=save, savetype=savetype)
    wm.plot_k(save=save, savetype=savetype)
    if show:
        plt.show()

if __name__ == "__main__":
    if not os.path.isdir("Figures"):
        os.makedirs("Figures")
    main()
    