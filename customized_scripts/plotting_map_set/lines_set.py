import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import aplpy
import matplotlib as mpl
import matplotlib.patches as patches
from astropy import units as u
import pyexcel_ods
import os

# Directory and filename
paths_fits = []
paths_fits.append("/home/mikh_kate/kalenskii/CASA/integrated_mol/CH3OH/CH3OH_241904_integrated_32.fits")
paths_fits.append("/home/mikh_kate/kalenskii/CASA/integrated_mol/CH3OH/CH3OH_243915_integrated_34.fits")
#paths_fits.append("/home/mikh_kate/kalenskii/CASA/integrated_mol/CH3OH/CH3OH_241700_integrated_22.fits")
#paths_fits.append("/home/mikh_kate/kalenskii/CASA/integrated_mol/CH3OH/CH3OH_241767_integrated_23.fits")
#paths_fits.append("/home/mikh_kate/kalenskii/CASA/integrated_mol/CH3OH/CH3OH_241791_integrated_24.fits")
#paths_fits.append("/home/mikh_kate/kalenskii/CASA/integrated_mol/CH3OH/CH3OH_241879_integrated_30.fits")

names = []
names.append(("CH3OH", '-'))
names.append(("CH3OH", '-'))
#names.append(("CH3OH", '-'))
#names.append(("CH3OH", '-'))
#names.append(("CH3OH", '-'))

latex_names = []
latex_names.append(("CH$_{{3}}$OH", 241904.643))
latex_names.append(("CH$_{{3}}$OH", 243915.811))
#latex_names.append(("CH$_{{3}}$OH", 241700.168))
#latex_names.append(("CH$_{{3}}$OH", 241767.247))
#latex_names.append(("CH$_{{3}}$OH", 241791.367))
#latex_names.append(("CH$_{{3}}$OH", 241879.038))



#latex_names.append(("CN", ))
#latex_names.append(("CH$_{{3}}$OH", 241879.038))

cols = 2
rows = (len(latex_names) + 1) // cols

beam = True 
markers = True
levels_ = True

 

# Modify if you want to change the design

cmap = 'YlOrRd' # Define colormap
#cmap = 'Greys' # Define colormap
#num_of_chan = rows * cols
cmapmin = -0.2 #minimum value of colormap
#cmapmax = 10 #maximum value of colormap, to use it, change vmax in ax.show_colorscale to cmapmax
contours_percent = [0.15, 0.45, 0.75] #percentage of the max value for contours
shift_x = 1.2 # number in absolute value
shift_y = 0.8 # number in absolute value
size_x = 3 #size of the cells
size_y = 3 #size of the cells
gap = 0.01 #the gap between the cells

def get_ra(hours, minutes, seconds):
    return (hours + ((minutes + (seconds / 60)) / 60))/24*360
def get_dec(degrees, minutes, seconds):
    return degrees+((minutes+(seconds/60))/60)

# Calculations
c = 299792 #speed of light
figsize_x = cols*size_x + shift_x + 0.45 + 1.4 #9
figsize_y = rows*size_y + shift_y + 0.8 #8
cell_size_x = size_x / figsize_x 
cell_size_y = size_y / figsize_y 
gap_x = cell_size_x * (gap + 0.35) #0.4 for colorbar
gap_y = cell_size_y * (gap + 0.1)
shift_x_norm = shift_x / figsize_x 
shift_y_norm = shift_y / figsize_y 

fig = plt.figure(figsize=(figsize_x, figsize_y))


for i, fits_file in enumerate(paths_fits, start=0):  # Start from 0 to match array indexing
    row = i // cols  # Calculate row index
    col = i % cols   # Calculate column index
    #freq = frequencies[channel]  # Indexing starts from 0
    x_coord = (cell_size_x + gap_x) * col  + shift_x_norm # x-coordinate of cell left edge
    y_coord = (cell_size_y + gap_y) * (rows - 1 - row) + shift_y_norm # y-coordinate of cell bottom edge
    #x_coord_cbar = (cell_size_x + gap_x) * (col + 1) + shift_x_norm # x-coordinate of colorbar left edge
    
    data = fits.getdata(fits_file)
    header = fits.getheader(fits_file)    
    
    ax = aplpy.FITSFigure(fits_file, slices=[0], figure=fig, subplot=[x_coord, y_coord, cell_size_x, cell_size_y])

    if (levels_ == True):
        max_value = np.max(data)
        levels = max_value * np.array(contours_percent)
        ax.show_contour(data[0], levels=levels, colors='black', linewidths=0.5, zorder=1)#  cmap = 'Greys_r'
    
    # Show the data with the specified colormap
    ax.show_colorscale(cmap=cmap)      #vmin = cmapmin, vmax = max_value, 

    # Add title on each map
    #veloсity = (ref_freq - freq/1e6) * c / ref_freq
    ax.add_label(0.5, 1.03, f'{latex_names[i][0]} {str(latex_names[i][1])} МГц', relative=True, color='black', size='large') #fontweight='bold')
  
    ax.ticks.set_tick_direction('in')
    ax.ticks.set_minor_frequency(1)

    
    #Hide axis labels inside the figure
    if (col > 0):
    #    ax.tick_labels.hide_y()
        ax.axis_labels.hide_y()
    if (row < rows - 1):
    #    ax.tick_labels.hide_x()
        ax.axis_labels.hide_x()    
    
    if (markers == True):
        # Define a celestial coordinates for sources 
        ra = [93.225, get_ra(6, 12, 53.77500), get_ra(6, 12, 53.84300)]
        dec = [17.989755555555554, get_dec(17, 59, 26.17), get_dec(17, 59, 23.6200)]
        ax.show_markers(ra,dec,marker='^',c='black',s=40)
    
    if (beam == True):
        
        #beam parameters
        print(fits_file)
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
        ax.ax.add_patch(rectangle) #The function built into aplpy doesn't work for our data, so we use matplotlib directly
        
        #show beam
        ax.show_ellipses(rect_x + rect_width/2, rect_y + rect_height/2, BMIN*2, BMAJ*2, BPA, coords_frame='pixel', zorder = 3, facecolor = 'black')

    # Add colorbar near each map
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    ax.add_colorbar()
    
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'    
    if (col == cols-1):
        ax.colorbar.set_axis_label_text("Jy/beam")
    #ax.colorbar.set_box([x_coord_cbar, y_coord, 0.2/ figsize_x ,cell_size_y])

# Save the plot
name_for_file = ""
for i in range(len(names)):
    name_for_file = name_for_file + names[i][0] + '_'
    
print(name_for_file)

path_save = '/home/mikh_kate/kalenskii/CASA/map_sets/'
os.makedirs(path_save, exist_ok=True)
output_file = path_save + name_for_file + '.pdf'
plt.savefig(output_file)

# Show the plot
plt.show()