# -*- coding: utf-8 -*-
"""
This script runs processing functions from `Modules.processing`.

"""

from __future__ import print_function
import os
if os.getcwd()[-7:] == "Scripts":
    print("Changing working directory to experiment root directory")
    os.chdir("../")
from Modules.processing import *

if __name__ == "__main__":
    print("Running processing script")

    # Dealing with individual runs
    r = Run("Wake-0.4", 20)
    r.plot_perf("cp")
    r.plot_wake()
    r.print_perf_stats()
    r.print_wake_stats()

    # Tare drag and torque"""
#    process_tare_torque(2, plot=True)
#    batch_process_tare_torque(plot=True)
#    process_tare_drag(5, plot=True)
#    batch_process_tare_drag(plot=True)
    
    # Batch processing
#    batch_process_section("Perf-0.3")
#    batch_process_all()