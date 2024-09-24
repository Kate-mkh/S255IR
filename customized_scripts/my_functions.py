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
#from astropy.utils.exceptions import AstropyWarning, FITSFixedWarning


def get_ra(hours, minutes, seconds):
    return (hours + ((minutes + (seconds / 60)) / 60))/24*360
def get_dec(degrees, minutes, seconds):
    return degrees+((minutes+(seconds/60))/60)

def find_lines(file_path, value):
    data = pyexcel_ods.get_data(file_path)
    sheet = data[list(data.keys())[0]]  # Assuming there's only one sheet

    version = 1

    for row in sheet:
        if len(row) >= 6 and row[5] == value:  # Check if row has at least 4 columns
            molecule = row[5] if len(row) >= 6 else None  # column E is the 5th column (index 4)
            file = row[0] 
            first_chan = row[7]
            num_of_chans = row[9]
            ref_freq = row[4]
            round_freq = row[3]
            yield version, molecule, ref_freq, file, first_chan, num_of_chans, round_freq
            version += 1

def find_all_lines(file_path):
    data = pyexcel_ods.get_data(file_path)
    sheet = data[list(data.keys())[0]]  # there's only one sheet

    version = 1

    for row in sheet[1:]:
        if len(row) >= 6:  # Check if row has at least 4 columns
            molecule = row[5] if len(row) >= 6 else None  # Assuming column E is the 5th column (index 4)
            ref_freq = row[4]
            file = row[0] 
            first_chan = row[7]
            num_of_chans = row[9]
            round_freq = row[3]
            yield version, molecule, ref_freq, file, first_chan, num_of_chans, round_freq
            version += 1
  
def plot_channel_map(version_f, molecule_f, ref_freq_f, file_f, first_chan_f, num_of_chans_f, round_freq_f, cols):
    
    warnings.filterwarnings("ignore", category=UserWarning, module="aplpy")
    warnings.filterwarnings("ignore", category=UserWarning, module="astropy.wcs.wcs")
    
    
    version = str(file_f)
    molecule = molecule_f
    ref_freq = ref_freq_f
    round_freq = round_freq_f
    dir_ = '/home/mikh_kate/kalenskii/CASA/s255all/s255ch' + str(file_f) + '/'
    fits_file_name = 's255ch' + str(file_f) + '.fits'
    first_chan = first_chan_f
    num_of_chans = num_of_chans_f
    print(f"Version: {version}, Molecule: {molecule}, File: {fits_file_name}, first channel: {first_chan}")
    

    warnings.filterwarnings("ignore")   
    

    
    beam = True #True - if you need the beam, False - if you don't
    
     
    
    # Modify if you want to change the design
    
    cmap = 'YlOrRd' # Define colormap
    #cmap = 'Greys' # Define colormap
    num_of_chan = num_of_chans
    rows = (num_of_chans + 1) // cols
        
    
    
    cmapmin = -0.2 #minimum value of colormap
    #cmapmax = 10 #maximum value of colormap, to use it, change vmax in ax.show_colorscale to cmapmax
    contours_percent = [0.25, 0.5, 0.75] #percentage of the max value for contours
    shift_x = 1.5 # number in absolute value
    shift_y = 1.2 # number in absolute value
    size_x = 2 #size of the cells
    size_y = 2 #size of the cells
    gap = 0.01 #the gap between the cells
    
    
    # Calculations
    c = 299792 #speed of light
    figsize_x = cols*size_x + shift_x + 0.25 + 1.4 #9
    figsize_y = rows*size_y + shift_y + 0.8 #8
    cell_size_x = size_x / figsize_x 
    cell_size_y = size_y / figsize_y 
    gap_x = cell_size_x * gap
    gap_y = cell_size_y * gap
    shift_x_norm = shift_x / figsize_x 
    shift_y_norm = shift_y / figsize_y 
    
    
    fits_file = dir_ + fits_file_name
    data = fits.getdata(fits_file)
    header = fits.getheader(fits_file)
    beam_data = fits.getdata(fits_file, ext=1)
    BMAJ = beam_data['BMAJ'][first_chan]  # Assuming the first line is line 0
    BMIN = beam_data['BMIN'][first_chan]   
    BPA = beam_data['BPA'][first_chan]    
    
    channels = range(first_chan, first_chan + num_of_chan)
    
    
    
    # Calculate frequency using header information
    n_channels, n_y, n_x = data.shape
    channel_width = header['CDELT3']  # Assuming frequency increment per channel
    channel_start = header['CRVAL3']  # Starting frequency
    frequencies = [channel_start + (ch * channel_width) for ch in range(n_channels)]
    
    fig = plt.figure(figsize=(figsize_x, figsize_y))
    
    
    for i, channel in enumerate(channels, start=0):  # Start from 0 to match array indexing
        row = i // cols  # Calculate row index
        col = i % cols   # Calculate column index
        freq = frequencies[channel]  # Indexing starts from 0
        x_coord = (cell_size_x + gap_x) * col  + shift_x_norm # x-coordinate of cell left edge
        y_coord = (cell_size_y + gap_y) * (rows - 1 - row) + shift_y_norm # y-coordinate of cell bottom edge
        x_coord_cbar = (cell_size_x + gap_x) * (col + 1) + shift_x_norm # x-coordinate of colorbar left edge
        if i != len(channels)-1 :
            warnings.filterwarnings("ignore")
            ax = aplpy.FITSFigure(fits_file, slices=[channel], figure=fig, subplot=[x_coord, y_coord, cell_size_x, cell_size_y])
            warnings.filterwarnings("ignore")
        if i == len(channels)-1 :
            warnings.filterwarnings("ignore")
            ax = aplpy.FITSFigure(fits_file, slices=[channel], figure=fig, subplot=[x_coord, y_coord, cell_size_x + 0.25/figsize_x , cell_size_y]) #we need to increase the x-size because of the colorbar
            warnings.filterwarnings("ignore")
         
    
        # Add contours (percentage of the max value)
        max_value = np.max(data[first_chan: first_chan + num_of_chan])
        #max_value = 4
        levels = max_value * np.array(contours_percent)
        ax.show_contour(data[channel], levels=levels, colors='black', linewidths=0.5, zorder=1)#  cmap = 'Greys_r' 
        #print(max_value, end = '\n')
        
        # Show the data with the specified colormap
        ax.show_colorscale(vmin = cmapmin, vmax = max_value, cmap=cmap)      
    
        # Add veloсity on each map
        veloсity = (ref_freq - freq/1e6) * c / ref_freq
        ax.add_label(0.3, 0.9, f'{veloсity:.2f} km s$^{{-1}}$', relative=True, color='black', size='large') #fontweight='bold')
      
        ax.ticks.set_tick_direction('in')
        #ax.ticks.set_linewidth(2)
        #ax.ticks.set_color('white')
        ax.ticks.set_minor_frequency(1)
        #ax.frame.set_linewidth(2)
        #ax.frame.set_color('white')
        
        #Hide axis labels inside the figure
        if (col > 0):
            ax.tick_labels.hide_y()
            ax.axis_labels.hide_y()
        if (row < rows - 1):
            ax.tick_labels.hide_x()
            ax.axis_labels.hide_x()    
        
        # Define a celestial coordinates for sources 
        ra = [93.225, get_ra(6, 12, 53.77500), get_ra(6, 12, 53.84300)]
        dec = [17.989755555555554, get_dec(17, 59, 26.17), get_dec(17, 59, 23.6200)]
    
        ax.show_markers(ra,dec,marker='^',c='black',s=20)
        #ax.add_colorbar()
        
        if (row == rows - 1) & (col == 0) & (beam == True):
            # Define rectangle properties in pixel coordinates
            rect_x = 2.5  
            rect_y = 3  
            rect_width = BMAJ*3       
            rect_height = BMAJ*3      
            
            # Create the rectangle patch for the beam
            rectangle = patches.Rectangle((rect_x, rect_y), rect_width, rect_height, linewidth=1, edgecolor='black', facecolor='white')
            
            # Add the rectangle patch to the plot
            ax.ax.add_patch(rectangle) #The function built into aplpy doesn't work for our data, so we use matplotlib directly
            
            #ax.add_beam() #it works if beam parametrs are in the header and if comment error lines in popup file :) But just ellipse is better.
            #ax.beam.set_color('black')
            #ax.beam.set_hatch('+') # What does it even do?..
            #ax.add_scalebar(0.2) #doesn't work with my data 
            
            #show beam
            ax.show_ellipses(rect_x + rect_width/2, rect_y + rect_height/2, BMIN*2, BMAJ*2, BPA, coords_frame='pixel', zorder = 3, facecolor = 'black')
            
    ax.ticks.set_tick_direction('out')
       
       
    
    
    
    ax.add_label(-0.5, rows + 0.2, molecule + ' ' + str(round_freq), relative=True, color='black', size='large', fontweight='bold')
    #ax.add_label(-0.5, 2.2, molecule + ' ' + str(round_freq), relative=True, color='black', size='large', fontweight='bold')
    
    
    # Add colorbar near the last map
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    ax.add_colorbar()
    ax.colorbar.set_axis_label_text("Jy/beam")
    ax.colorbar.set_box([x_coord_cbar, y_coord, 0.2/ figsize_x ,cell_size_y])
    #ax.colorbar.set_ticks([0, 5, 10]) #if you need specific values displayed
    
    # Save the plot    
    path_save = '/home/mikh_kate/kalenskii/CASA/channel_maps/'+ molecule + '/'
    os.makedirs(path_save , exist_ok=True)    
    
    output_file_2 = path_save + molecule + '_' + str(round_freq) + '_' + str(version) + '.pdf'
    plt.savefig(output_file_2)
    
    # Show the plot
    #plt.show    ()
    
def plot_map(new_fits_file_name, molecule, ref_freq, round_freq, dir_integrated, version):
    label = True
    colorbar_ = True
    levels_ = True
    markers = True
    beam = True #True - if you need the beam, False - if you don't.
    
    contours_percent = [0.25, 0.5, 0.75] #percentage of the max value for contours
    
    # Open the FITS file
    with fits.open(new_fits_file_name) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    
    fig = plt.figure(figsize=(9, 9))
    
    f = aplpy.FITSFigure(new_fits_file_name, figure=fig, subplot=[0.15, 0.1, 0.75, 0.75])
    f.show_colorscale(cmap= 'YlOrRd') 
    
    if (label == True):
        f.add_label(0.45, 1.1, molecule + ' ' + str(ref_freq), relative=True, color='black', size='large', fontweight='bold')
    
    if (colorbar_ == True):
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        f.add_colorbar()
        f.colorbar.set_axis_label_text("Jy/beam")
        mpl.rcParams['xtick.direction'] = 'out'
        mpl.rcParams['ytick.direction'] = 'out'        
    
    if (levels_ == True):
        max_value = np.max(data)
        levels = max_value * np.array(contours_percent)
        f.show_contour(data[0], levels=levels, colors='black', linewidths=0.5, zorder=1)#  cmap = 'Greys_r'
    
    if (markers == True):
        # Define a celestial coordinates for sources 
        ra = [93.225, get_ra(6, 12, 53.77500), get_ra(6, 12, 53.84300)]
        dec = [17.989755555555554, get_dec(17, 59, 26.17), get_dec(17, 59, 23.6200)]
        f.show_markers(ra,dec,marker='^',c='black',s=80)
    
    if (beam == True):
        #beam parameters
        BMAJ = header['BMAJ']
        BMIN = header['BMIN']
        BPA = header['BPA']
        
        # Define rectangle properties in pixel coordinates
        rect_x = 2.5  
        rect_y = 3  
        rect_width = BMAJ*3       
        rect_height = BMAJ*3      
        
        # Create the rectangle patch for the beam
        rectangle = patches.Rectangle((rect_x, rect_y), rect_width, rect_height, linewidth=1, edgecolor='black', facecolor='white')
        
        # Add the rectangle patch to the plot
        f.ax.add_patch(rectangle) #The function built into aplpy doesn't work for our data, so we use matplotlib directly
        
        #show beam
        f.show_ellipses(rect_x + rect_width/2, rect_y + rect_height/2, BMIN*2, BMAJ*2, BPA, coords_frame='pixel', zorder = 3, facecolor = 'black')
        
        
    
    output_file = dir_integrated  + molecule + '_' + str(round_freq) + '_integrated_' + str(version) + '.pdf'
    os.makedirs(dir_integrated , exist_ok=True) 
    plt.savefig(output_file)
    #plt.show()
