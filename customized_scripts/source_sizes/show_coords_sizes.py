import re
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

# Script to get data from a text file

def parse_file(input_file, output_file):
    ra = dec = ra_err = dec_err = majaxis = minaxis = PA = majaxis_err = minaxis_err = PA_err = None
    with open(input_file, 'r') as f:
        lines = f.readlines()

    ra_pattern = re.compile(r'ra:\s+([0-9:.]+)\s+\+/-\s+([0-9.]+)')
    dec_pattern = re.compile(r'dec:\s+([+\d:.]+)\s+\+/-\s+([0-9.]+)')
    
    for line in lines:
        ra_match = ra_pattern.search(line)
        if ra_match:
            ra = ra_match.group(1)
            ra_err = ra_match.group(2)
            break

    for line in lines:
        dec_match = dec_pattern.search(line)
        if dec_match:
            dec = dec_match.group(1).lstrip('+')  # I do not need "+"
            dec_err = dec_match.group(2)
            break

    # major, minor and position angle
    for i, line in enumerate(lines):
        if 'Image component size (deconvolved from beam) ---' in line:
            next_line = lines[i + 1].strip()
            if "Could not deconvolve source from beam" in next_line:
                majaxis = "Could not deconvolve source from beam"
                minaxis = PA = ""
                majaxis_err = minaxis_err = PA_err = None
            else:
                majaxis_pattern = re.compile(r'major axis FWHM:\s+([0-9.]+)\s+\+/-\s+([0-9.]+)')
                minaxis_pattern = re.compile(r'minor axis FWHM:\s+([0-9.]+)\s+\+/-\s+([0-9.]+)')
                pa_pattern = re.compile(r'position angle:\s+([0-9.]+)\s+\+/-\s+([0-9.]+)')
                
                # major axis
                majaxis_match = majaxis_pattern.search(lines[i+1])
                if majaxis_match:
                    majaxis = majaxis_match.group(1)
                    majaxis_err = majaxis_match.group(2)
                
                # minor axis
                minaxis_match = minaxis_pattern.search(lines[i+2])
                if minaxis_match:
                    minaxis = minaxis_match.group(1)
                    minaxis_err = minaxis_match.group(2)
                
                # position angle
                pa_match = pa_pattern.search(lines[i+3])
                if pa_match:
                    PA = pa_match.group(1)
                    PA_err = pa_match.group(2)
    
    result = f'{ra}\t{ra_err}\t{dec}\t{dec_err}\t{majaxis}\t{majaxis_err}\t{minaxis}\t{minaxis_err}\t{PA}\t{PA_err}\n'
    return result
    
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
            yield path, imagename, logfile, molecule, round_freq


input_file = '/home/mikh_kate/kalenskii/CASA/integrated_mol/CH3CHO/CH3CHO_234826_integrated_14337_2.txt'
output_file = 'sorces.txt'

with open(output_file, 'w') as f:
    f.write("Molecule\tFreq\tRA\tRA_err\tDec\tDec_err\tMajAxis\tMajAxis_err\tMinAxis\tMinAxis_err\tPA\tPA_err\n")    
    for path, imagename, logfile, molecule, freq in find_lines('S255IR_lines.ods'):
        f.write(molecule+ '\t' + str(freq) + '\t' + parse_file(path+logfile, output_file))

