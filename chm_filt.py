#!/usr/bin/env python
#
# Utility to filter the HRSI CHMs
# 

import sys
import os
import subprocess
import glob
import argparse
import shutil

import numpy as np

from pygeotools.lib import iolib
from pygeotools.lib import malib
from pygeotools.lib import geolib
from pygeotools.lib import filtlib
from pygeotools.lib import warplib

import matplotlib
##https://stackoverflow.com/questions/37604289/tkinter-tclerror-no-display-name-and-no-display-environment-variable
matplotlib.use('Agg')
import matplotlib.pyplot, matplotlib.mlab, math
import scipy.stats

def main():
    pass

"""
First, merge 2 versions of dem_control

Then, write a chm_filt.py

import dem_control

# Basic logic
Divide HRSI CHM into forest and non-forest;
to estimate max canopy height, within the forest mask run a 'max'filter (filtlib)
to remove spurious 'heights' in the non-forest using a 'min' filter (filtlib)


(1) use roughmask (as is now in dem_control) to get 'non-forest' valid pixels
    -run a 'min' filter,

(2) invert roughmask to get 'forest'
    -run a 'max' filter,

(3) then mask the result with the toamask (this removes water and other dark (shadow) areas

(4) for later: then mask with the slopemask, toatrimask
"""
if __name__ == '__main__':
    main()
