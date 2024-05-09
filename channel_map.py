import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import aplpy
import matplotlib as mpl
import matplotlib.patches as patches
from astropy import units as u

#input data

dir_ = '/home/mikh_kate/kalenskii/CASA/s255all/s255ch18433/'
fits_file_name = 's255ch18433.fits'
first_chan = 1375
ref_freq = 247228.737 #line frequency 
molecule = 'CH3OH'
version = 1

rows = 3
cols = 3
beam = True #True - if you need the beam, False - if you don't.
#beam parameters
BMAJ    =   4.318158    # you can find it in the header or 
BMIN    =   3.294571    # in the fv software in [1] in "all" 
BPA     =   -0.3356704  #

 

# Modify if you want to change the design

cmap = 'YlOrRd' # Define colormap
#cmap = 'Greys' # Define colormap
num_of_chan = rows * cols
cmapmin = -0.2 #minimum value of colormap
#cmapmax = 10 #maximum value of colormap, to use it, change vmax in ax.show_colorscale to cmapmax
contours_percent = [0.15, 0.45, 0.75] #percentage of the max value for contours
shift_x = 1.5 # number in absolute value
shift_y = 1.2 # number in absolute value
size_x = 2 #size of the cells
size_y = 2 #size of the cells
gap = 0.01 #the gap between the cells

def get_ra(hours, minutes, seconds):
    return (hours + ((minutes + (seconds / 60)) / 60))/24*360
def get_dec(degrees, minutes, seconds):
    return degrees+((minutes+(seconds/60))/60)

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
        ax = aplpy.FITSFigure(fits_file, slices=[channel], figure=fig, subplot=[x_coord, y_coord, cell_size_x, cell_size_y])
    if i == len(channels)-1 :
        ax = aplpy.FITSFigure(fits_file, slices=[channel], figure=fig, subplot=[x_coord, y_coord, cell_size_x + 0.25/figsize_x , cell_size_y]) #we need to increase the x-size because of the colorbar
     

    # Add contours (percentage of the max value)
    max_value = np.max(data[first_chan: first_chan + num_of_chan])
    #max_value = 4.5
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



ax.add_label(-0.5, rows + 0.2, molecule + ' ' + str(ref_freq), relative=True, color='black', size='large', fontweight='bold')
#ax.add_label(-0.5, 2.2, molecule + ' ' + str(ref_freq), relative=True, color='black', size='large', fontweight='bold')


# Add colorbar near the last map
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
ax.add_colorbar()
ax.colorbar.set_axis_label_text("Jy/beam")
ax.colorbar.set_box([x_coord_cbar, y_coord, 0.2/ figsize_x ,cell_size_y])
#ax.colorbar.set_ticks([0, 5, 10]) #if you need specific values displayed

# Save the plot
output_file = dir_ + molecule + '_' + str(int(ref_freq//1)) +'.pdf'
plt.savefig(output_file)

output_file_2 = '/home/mikh_kate/kalenskii/CASA/All_the_maps/' + molecule + '_' + str(int(ref_freq//1)) + '_' + str(version) + '.pdf'
plt.savefig(output_file_2)

# Show the plot
plt.show()
