# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:03:29 2014

@author: Pete
"""

from Modules.processing import *
from Modules import test

def plot_meancomboquiv(U_infty=1.0):
    wm = WakeMap(U_infty)
    wm.plot_meancomboquiv()

if __name__ == "__main__":
    """Choose functions to run here."""
    plt.close("all")
    p = "Google Drive/Research/Papers/CFT Re dependence/figures"
    p = os.path.join(os.path.expanduser("~"), p)

    """Tests"""
#    test.test_all()

    """Dealing with individual runs"""
#    r = Run("Wake-1.0", 50)
#    r.calcperf()
#    r.calcwake()
#    r.plotperf()
#    r.plotwake()

    """Tare drag and torque"""
#    process_tare_torque(2, plot=True)
#    batch_process_tare_torque(plot=True)

#    process_tare_drag(5, plot=True)
#    batch_process_tare_drag(plot=True)
#    plot_tare_drag()
    
    """Batch processing"""
#    batch_process_section("Perf-1.0", reprocess=True)
#    batch_process_all()
    
    """Plotting"""
    plot_perf_curves(save=False, savepath=p)
#    plot_perf_re_dep(save=False, cfd=False, savepath=p, normalize_by="default",
#                     dual_xaxes=True)
#    plot_wake_profiles(z_H=0.0, save=True, savepath=p)
#    plot_settling(1.0)
#    plot_meancomboquiv(U_infty=0.4)
#    plot_meancomboquiv(U_infty=1.2)
