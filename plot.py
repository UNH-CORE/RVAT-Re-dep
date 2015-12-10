#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script generates all relevant figures from the experiment and
stores them in the `Figures` directory.
"""

from pyrvatrd.plotting import *
import argparse

kwargs = {"markerfacecolor": "none",
          "color": "black"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create figures from the "
                                     "UNH-RVAT Reynolds number dependence "
                                     "experiment")
    parser.add_argument("figures", nargs="?", help="Which figures to create",
                        choices=["perf_curves", "perf_re_dep", "wake_profiles",
                        "k_bar_graph", "mom_bar_graph", "all_meancontquiv",
                        "all_kcont", "multi_spec", "wake_trans_totals",
                        "vel_unc_table"],
                        default=[])
    parser.add_argument("--all", "-A", help="Plot all figures",
                        action="store_true", default=False)
    parser.add_argument("--no-errorbars", "-e", help="Do not plot error bars",
                        action="store_true", default=False)
    parser.add_argument("--save", "-s", help="Save figures",
                        action="store_true", default=False)
    parser.add_argument("--savetype", help="Format to save figures",
                        default=".pdf")
    parser.add_argument("--style", help="Matplotlib stylesheet")
    parser.add_argument("--noshow", help="Do not show figures",
                        action="store_true", default=False)
    args = parser.parse_args()

    if args.style:
        plt.style.use(args.style)
    else:
        from pxl.styleplot import set_sns
        set_sns()
    savetype = args.savetype
    save = args.save
    errorbars = not args.no_errorbars
    if save:
        if not os.path.isdir("Figures"):
            os.makedirs("Figures")

    if "perf_curves" in args.figures or args.all:
        plot_perf_curves(subplots=False, save=save, savetype=savetype, **kwargs)
        plot_perf_curves(save=save, savetype=savetype, **kwargs)
    if "perf_re_dep" in args.figures or args.all:
        plot_perf_re_dep(errorbars=errorbars, save=save, savetype=savetype,
                         **kwargs)
        plot_perf_re_dep(subplots=False, errorbars=errorbars, save=save,
                         savetype=savetype, **kwargs)
    if "wake_profiles" in args.figures or args.all:
        plot_wake_profiles(save=save, savetype=savetype)
    if "k_bar_graph" in args.figures or args.all:
        make_k_bar_graph(save=save, savetype=savetype)
    if "mom_bar_graph" in args.figures or args.all:
        make_mom_bar_graph(save=save, savetype=savetype)
    if "all_meancontquiv" in args.figures:
        plot_all_meancontquiv(save=save, savetype=savetype)
    if "all_kcont" in args.figures:
        plot_all_kcont(save=save, savetype=savetype)
    if "multi_spec" in args.figures or args.all:
        plot_multi_spec(plot_conf_int=errorbars, save=save, savetype=savetype)
    if "wake_trans_totals" in args.figures or args.all:
        plot_wake_trans_totals(save=save, savetype=savetype)
    if "vel_unc_table" in args.figures or args.all:
        make_velocity_unc_table(save=save)

    if not args.noshow:
        plt.show()
