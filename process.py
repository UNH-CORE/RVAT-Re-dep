#!/usr/bin/env python
"""Process data."""

from __future__ import print_function
from pyrvatrd.processing import *
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process data from the "
                                     "UNH-RVAT Reynolds number dependence "
                                     "experiment.")
    parser.add_argument("--single-run", "-r", nargs=2, help="Process a single "
                        "run; args=[section, nrun]")
    parser.add_argument("--single-tare-drag", help="Process a tare drag run")
    parser.add_argument("--single-tare-torque",
                        help="Process a tare torque run")
    parser.add_argument("--section", help="Process a test matrix section")
    parser.add_argument("--tare-drag", help="Process tare drag runs",
                        action="store_true")
    parser.add_argument("--tare-torque", help="Process tare torque runs",
                        action="store_true")
    parser.add_argument("--all", "-a", action="store_true",
                        help="Process all data")
    parser.add_argument("--plot", action="store_true", default=False,
                        help="Create plots (if applicable)")
    args = parser.parse_args()

    if args.plot:
        from pxl.styleplot import set_sns
        set_sns()

    if args.single_run:
        # Dealing with individual runs
        section = args.single_run[0]
        nrun = int(args.single_run[1])
        print("Processing {} run {}".format(section, nrun))
        r = Run(section, nrun)
        if args.plot:
            r.plot_perf("cp")
            r.plot_wake()
        print(r.summary)

    if args.single_tare_drag:
        process_tare_drag(args.single_tare_drag, plot=args.plot)

    if args.tare_drag:
        batch_process_tare_drag(plot=args.plot)

    if args.single_tare_torque:
        process_tare_torque(args.single_tare_torque, plot=args.plot)

    if args.tare_torque:
        batch_process_tare_torque(plot=args.plot)

    if args.section:
        print("Processing {}".format(args.section))
        batch_process_section(args.section)

    if args.all:
        batch_process_all()

    if args.plot:
        plt.show()
