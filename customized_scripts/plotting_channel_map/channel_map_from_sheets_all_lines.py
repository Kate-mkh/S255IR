import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import aplpy
import matplotlib as mpl
import matplotlib.patches as patches
from astropy import units as u
import pyexcel_ods
import os
import warnings

from my_functions import get_ra, get_dec, find_all_lines, plot_channel_map

warnings.filterwarnings("ignore") # working, but not with all the warnings
#You can run without warnings in terminal:
#python3 channel_map_from_sheets_all_lines.py 2> /dev/null


file_path = "/home/mikh_kate/kalenskii/CASA/lines_for_fixing.ods"


cols = 3

for i, (version, molecule, ref_freq, file, first_chan, num_of_chans) in enumerate(find_all_lines(file_path)):
    #if i >= 2: #for tests
    #    break  #for tests
    plot_channel_map(version, molecule, ref_freq, file, first_chan, num_of_chans, cols)
