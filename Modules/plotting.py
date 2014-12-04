# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 09:15:20 2014

@author: Pete
"""

import matplotlib
import os

style_file = "Config/ggplot-mod.mplstyle"

try:
    matplotlib.style.use(style_file)
except ValueError:
    matplotlib.style.use(os.path.join("../", style_file))

ylabels = {"mean_u" : r"$U/U_\infty$",
           "std_u" : r"$\sigma_u/U_\infty$",
           "mean_v" : r"$V/U_\infty$",
           "mean_w" : r"$W/U_\infty$",
           "mean_upvp" : r"$\overline{u'v'}/U_\infty^2$",
           "mean_u_diff" : r"$\Delta U$ (\%)",
           "mean_v_diff" : r"$\Delta V$ (\%)",
           "mean_w_diff" : r"$\Delta W$ (\%)"}