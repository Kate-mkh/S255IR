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

def find_lines(file_path):
    data = pyexcel_ods.get_data(file_path)
    sheet = data[list(data.keys())[0]]  # Assuming there's only one sheet

    version = 1

    for row in sheet:
        if len(row) >= 6 and (row[5] == 'Yes' or row[5] == 'Almost'):  # Check if row has at least 4 columns
            molecule = row[1]
            file = row[13] 
            round_freq = row[3]
            path = '/home/mikh_kate/kalenskii/CASA/integrated_mol/' + molecule + '/'
            imagename = molecule + '_' + str(round_freq) + '_integrated_' + str(file) + '.fits'
            logfile = molecule + '_' + str(round_freq) + '_integrated_' + str(file) + '.txt'
            yield path, imagename, logfile
            


with open('commands.sh', 'w') as f:
    for path, imagename, logfile in find_lines('S255IR_lines.ods'):
        casa_command = f"casa -c \"imfit(imagename='{path + imagename}', chans='0', logfile='{path + logfile}')\""
        f.write(casa_command + '\n')
