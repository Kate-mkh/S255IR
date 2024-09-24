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

from my_functions import get_ra, get_dec, find_lines, find_all_lines, plot_map

warnings.filterwarnings("ignore") # working, but not with all the warnings


file_path = "/home/mikh_kate/kalenskii/CASA/lines.ods"
path_save = '/home/mikh_kate/kalenskii/CASA/integrated_mol/'

for version_f, molecule_f, ref_freq_f, file_f, first_chan_f, num_of_chans_f, round_freq in find_all_lines(file_path):
    version = version_f
    #if version == 3:
    #    break
    molecule = molecule_f
    dir_ = '/home/mikh_kate/kalenskii/CASA/s255all/s255ch' + str(file_f) + '/'
    fits_file_name = 's255ch' + str(file_f) + '.fits'
    first_chan = first_chan_f
    num_of_chans = num_of_chans_f
    ref_freq = ref_freq_f
    print(f"Version: {version}, Molecule: {molecule}, File: {fits_file_name}, first channal: {first_chan}, {dir_}")

    #input data
    fits_file = dir_ + fits_file_name
    data = fits.getdata(fits_file)
    header = fits.getheader(fits_file)
    beam_data = fits.getdata(fits_file, ext=1)
    beam_header = fits.getheader(fits_file, ext=1)
    
    
    #INTEGRATING
    slices_to_keep = range(first_chan, first_chan + num_of_chans)
    channels = len(slices_to_keep)
    selected_slices = data[slices_to_keep]
    integrated = [selected_slices[0]]
    for i in range(1, channels):
        integrated[0] += selected_slices[i]
        
    new_header = header.copy()
    bmaj_value = beam_data['BMAJ'][first_chan]  
    bmin_value = beam_data['BMIN'][first_chan]   
    bpa_value = beam_data['BPA'][first_chan]       
    new_header.set('BMAJ', bmaj_value, 'Beam major axis FWHM [arcsec]')# for futher mapping
    new_header.set('BMIN', bmin_value, 'Beam minor axis FWHM [arcsec]')  
    new_header.set('BPA', bpa_value, 'Beam Position Angle [degrees]')
        
    beam_data_new = beam_data[first_chan:first_chan+1].copy()
    beam_data_new['CHAN'][0] = 0
    
    
    # new HDU from BEAMS
    hdu_beams = fits.BinTableHDU(data=beam_data_new, header=beam_header, name='BEAMS')
    
    hdu_beams.header.set('NCHAN', 1, 'Number of channels')
    
    dir_integrated = path_save + molecule + '/'
    new_fits_file_name = dir_integrated + molecule + '_' + str(round_freq) + '_integrated_' + str(file_f) + '.fits'
    print(new_fits_file_name)
    
    os.makedirs(dir_integrated , exist_ok=True) 
    
    # HDUList with PrimaryHDU and BEAMS HDU
    hdul = fits.HDUList([fits.PrimaryHDU(data=integrated, header=new_header), hdu_beams])   

    hdul.writeto(new_fits_file_name, overwrite=True)
    
    #PLOTTING 
    plot_map(new_fits_file_name, molecule, ref_freq, round_freq, dir_integrated, file_f)
    
