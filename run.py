# -*- coding: utf-8 -*-
"""
This script generates all relevant figures from the experiment and
stores them in the `Figures` directory.

"""

from Modules.plotting import *

sns.set(style="white", context="paper", font_scale=1.75,
        rc={"lines.markersize": 9, "lines.markeredgewidth": 1.25,
            "legend.fontsize": 12})
show = True
save = True
savetype = ".pdf"

def main():
    plot_cp_curve(1.0, save=save, savetype=savetype)
    plot_perf_curves(save=save, savetype=savetype)
    plot_perf_re_dep(save=save, savetype=savetype)
    plot_wake_profiles(z_H=0.0, save=save, savetype=savetype)
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
    