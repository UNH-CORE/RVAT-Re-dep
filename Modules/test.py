# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 08:51:05 2014

@author: Pete
"""
import os
if os.getcwd()[-7:] == "Modules":
    print("Changing working directory to experiment root directory")
    os.chdir("../")
from Modules.processing import *

def test_run():
    print("Testing Run class")
    run = Run("Perf-0.4", 20)
    run.calcperf()
    run.calcwake()
    print("PASS")
    
def test_perf_curve():
    print("Testing PerfCurve class")
    pc = PerfCurve(0.6)
    pc.plotcp()
    print("PASS")
    
def test_wake_profile():
    print("Testing WakeProfile class")
    wp = WakeProfile(0.6, 0.25, "horizontal")
    wp.plot("mean_u")
    print("PASS")
    
def test_all():
    test_run()
    test_perf_curve()
    print("Testing plot_perf_re_dep")
    plot_perf_re_dep()
    plot_perf_re_dep(dual_xaxes=True)
    print("PASS")
    print("Testing plot_perf_curves")
    plot_perf_curves()
    print("PASS")
    print("Testing plot_trans_wake_profile")
    plot_trans_wake_profile("mean_u")
    print("PASS")
    print("Testing plot_wake_profiles")
    plot_wake_profiles(z_H=0.0, save=False)
    print("PASS")
    test_wake_profile()
    print("Testing process_tare_torque")
    process_tare_torque(2, plot=False)
    print("PASS")
    print("Testing process_tare_drag")
    process_tare_drag(5, plot=False)
    print("PASS")
    print("All tests passed")
    
if __name__ == "__main__":
    test_all()
#    test_wake_profile()