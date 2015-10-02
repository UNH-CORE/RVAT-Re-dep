# -*- coding: utf-8 -*-
"""
This script generates all relevant figures from the experiment and
stores them in the `Figures` directory.

"""

from py_rvat_re_dep.plotting import *

set_sns()
show = True
save = True
savetype = ".pdf"

kwargs = {"markerfacecolor": "none",
          "color": "black"}

def main():
    plot_perf_curves(save=save, savetype=savetype, **kwargs)
    plot_perf_re_dep(save=save, savetype=savetype, **kwargs)
    plot_perf_re_dep(subplots=False, save=save, savetype=savetype, **kwargs)
    plot_wake_profiles(save=save, savetype=savetype)
    make_k_bar_graph(save=save, savetype=savetype)
    make_mom_bar_graph(save=save, savetype=savetype)
    plot_all_meancontquiv(save=save, savetype=savetype)
    plot_all_kcont(save=save, savetype=savetype)
    plot_multi_spec(save=save, savetype=savetype)
    plot_wake_trans_totals(save=save, savetype=savetype)
    if show:
        plt.show()

if __name__ == "__main__":
    if not os.path.isdir("Figures"):
        os.makedirs("Figures")
    main()
