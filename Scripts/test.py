# -*- coding: utf-8 -*-
"""
This script will run tests from `Modules.tests`.

"""
from __future__ import print_function
import os
if os.path.split(os.getcwd())[-1] == "Scripts":
    print("Changing working directory to experiment root directory")
    os.chdir("../")
from Modules import tests


if __name__ == "__main__":
    tests.test_all()
#    tests.test_download_raw()
#    tests.test_plot_settling()
