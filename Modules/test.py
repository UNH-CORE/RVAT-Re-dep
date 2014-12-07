# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 08:51:05 2014

@author: Pete
"""
import os
import time
if os.getcwd()[-7:] == "Modules":
    print("Changing working directory to experiment root directory")
    os.chdir("../")
from Modules.processing import *

def test_run():
    print("Testing Run class")
    run = Run("Wake-1.0", 20)
    print(run.cp_per_rev)
    print(run.std_cp_per_rev)
    print(run.cp_conf_interval)
    print(run.mean_cp)
    print(run.unc_cp)
    print(run.exp_unc_cp)
    run.print_perf_stats()
    print("PASS")
    
def test_section():
    print("Testing Section class")
    section = Section("Wake-1.0")
    print("PASS")
    
def test_batch_process_section():
    print("Testing batch_process_section")
    batch_process_section("Perf-1.0")
    df = pd.read_csv("Data/Processed/Perf-1.0.csv")
    print(df)
    plt.figure()
    plt.plot(df.mean_tsr, df.mean_cp)
    plt.show()
    
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
    
def test_wake_map():
    print("Testing WakeMap class")
    wm = WakeMap(0.4)
    wm.plot_meancomboquiv()
    wm2 = WakeMap(1.2)
    wm2.plot_meancomboquiv()
#    wm.plot_diff(quantity="mean_w", U_infty_diff=0.6)
#    wm.plot_meancomboquiv_diff(0.8, percent=False)
    print("PASS")
    
def test_process_section_parallel():
    nproc = 4
    nruns = 32
    t0 = time.time()
    s = Section("Wake-1.0")
    s.process_parallel(nproc=nproc, nruns=nruns)
    print("Parallel elapsed time: {} seconds".format(time.time() - t0))
    t0 = time.time()
    df = pd.DataFrame()
    for n in range(nruns):
        r = Run(s.name, n)
        df = df.append(r.summary, ignore_index=True)
    print("Serial elapsed time: {} seconds".format(time.time() - t0))
    assert(np.all(s.data.run == df.run))
    assert(np.all(s.data.mean_cp == df.mean_cp))
    assert(np.all(s.data.mean_cd == df.mean_cd))
    print("PASS")
    
def test_all():
    test_run()
    test_section()
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
    test_wake_map()
    plt.show()
    print("All tests passed")
    
if __name__ == "__main__":
#    test_run()
#    test_all()
#    test_wake_profile()
#    test_wake_map()
#    plot_perf_curves()
#    test_section()
#    test_batch_process_section()
    test_process_section_parallel()