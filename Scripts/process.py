# -*- coding: utf-8 -*-
"""
This script runs processing functions from `Modules.processing`.

"""

from __future__ import print_function
import os
if os.path.split(os.getcwd())[-1] == "Scripts":
    print("Changing working directory to experiment root directory")
    os.chdir("../")
from Modules.processing import *

if __name__ == "__main__":
    print("Running processing script")
    batch_process_all()