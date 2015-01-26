# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 16:03:29 2014

@author: Pete
"""

from Modules.processing import *
from Modules.plotting import *

if __name__ == "__main__":
    """Choose functions to run here."""
    plt.close("all")
    p = "Figures"
    if not os.path.isdir(p):
        os.makedirs(p)

    """Dealing with individual runs"""
    r = Run("Wake-0.4", 20)
    r.plot_perf("cp")
    r.plot_wake()
    
    """Performance curves"""
#    pc = PerfCurve(1.0)
#    pc.plotcp(save=True, savepath=p)

    """Tare drag and torque"""
#    process_tare_torque(2, plot=True)
#    batch_process_tare_torque(plot=True)
#    process_tare_drag(5, plot=True)
#    batch_process_tare_drag(plot=True)
#    plot_tare_drag()
    
    """Batch processing"""
#    batch_process_section("Perf-0.3")
#    batch_process_all()
    
    """Plotting"""
#    plot_perf_curves(save=False, savepath=p)
#    plot_perf_re_dep(save=False, cfd=False, savepath=p, normalize_by="default",
#                     dual_xaxes=True)
#    plot_wake_profiles(z_H=0.0, save=True, savepath=p)
#    plot_settling(1.0)
#    plot_meancomboquiv(U_infty=0.4)
#    plot_meancomboquiv(U_infty=1.2)
