# -*- coding: utf-8 -*-
"""
This module contains test functions.
"""

import os
import time
from .processing import *
from .plotting import *
from nose.tools import assert_almost_equal


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
    wm.plot_meancontquiv()
    wm2 = WakeMap(1.2)
    wm2.plot_meancontquiv()
#    wm.plot_diff(quantity="mean_w", U_infty_diff=0.6)
#    wm.plot_meancontquiv_diff(0.8, percent=False)
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

def test_batch_process_section_vs_parallel():
    name = "Perf-1.0"
    t0 = time.time()
    batch_process_section_old(name)
    print(time.time() - t0)
    t0 = time.time()
    Section(name).process()
    print(time.time() - t0)

def test_download_raw():
    """Tests the `processing.download_raw` function."""
    print("Testing processing.download_raw")
    # First rename target file
    fpath = "Data/Raw/Perf-1.0/0/metadata.json"
    fpath_temp = "Data/Raw/Perf-1.0/0/metadata-temp.json"
    exists = False
    if os.path.isfile(fpath):
        exists = True
        os.rename(fpath, fpath_temp)
    try:
        download_raw("Perf-1.0", 0, "metadata")
        # Check that file contents are equal
        with open(fpath) as f:
            content_new = f.read()
        if exists:
            with open(fpath_temp) as f:
                content_old = f.read()
            assert(content_new == content_old)
    except ValueError as e:
        print(e)
    os.remove(fpath)
    if exists:
        os.rename(fpath_temp, fpath)
    print("PASS")

def test_plot_settling():
    print("Testing plotting.plot_settling")
    plot_settling(1.0)
    print("PASS")

def test_calc_mom_transport():
    print("Testing WakeMap.calc_mom_transport")
    wm = WakeMap(1.0)
    wm.calc_mom_transport()
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
    test_download_raw()
    test_plot_settling()
    test_calc_mom_transport()
    print("All tests passed")

if __name__ == "__main__":
    pass
