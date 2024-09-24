from astropy.io import fits

def add_lines_to_header(fits_file_path):
    # Define the new header lines
    new_header = {'OBJECT': 'S255IR SMA2',
                  'EQUINOX': 2000.0,
                  'BUNIT': 'Jy/beam',
                  'TELESCOP': 'SMA',
                  'OBSERVER': 'syliu',
                  'DATE-OBS': '2020-10-06',
                  'UT': '14:57:22.2'}

    # Open the FITS file
    with fits.open(fits_file_path, mode='update') as hdul:
        # Get the primary header
        header = hdul[0].header
        
        
        # Update the header with new lines
        for key, value in new_header.items():
            header[key] = value        
        
        # Convert comma-separated numeric values to use dot as decimal separator
        #for key in header:
        #    print(header[key])
        #    if isinstance(header[key], str) and ',' in header[key]:
        #        print('bkbkjbkj')
        #        header[key] = float(header[key].replace(',', '.'))
        
    
        # Save the changes
        hdul.flush()
    
    print('Success')	
    
set_of_files = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

for i in set_of_files:
    fits_file_path = '/home/mikh_kate/kalenskii/CASA/SMA2/s255ch' + str(i) + '_2_SMA2.fits'
    print(fits_file_path)
    add_lines_to_header(fits_file_path)