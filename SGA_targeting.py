import os
import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.widgets import TextBox, Button
from PIL import Image
import json
import io
from urllib.request import urlopen

from astropy.table import Table



pix_scale = 0.25 # arcsec/pixel
pix_scale_arcmin = pix_scale/60

fiber_diameter = 1.52 # arcsec
fiber_diameter_pixels = fiber_diameter/pix_scale

patrol_radius = 1.4 # arcmin
patrol_radius_pixels = patrol_radius/pix_scale_arcmin



class GalaxyMarker(object):
   
    def __init__(self, input_filename):

        # Get index number of current file
        filename_parts = input_filename.split('.')
        filename_subparts = filename_parts[0].split('_')
        file_index = filename_subparts[-1]
       
        self.out_filename = "SGA_off-axis_targets_" + file_index + ".txt"
        
        # Initialize output target dictionary
        self.file_output = {}

        self.input_table = Table.read(input_filename, format='fits')
       
        # Get all the objects in the specified file
        self.objects = self.input_table['GALAXY']

        # Initialize list of image urls
        self.files = []
       
        for i in range(len(self.input_table)):

            # Extract sky coordinates for object
            ra = self.input_table['SGA_RA'][i]
            dec = self.input_table['SGA_DEC'][i]

            ####################################################################
            # Determine size of image needed
            #-------------------------------------------------------------------
            major_axis = self.input_table['DIAM'][i]

            major_axis_pixels = major_axis/pix_scale_arcmin

            img_size = int(major_axis_pixels + 100)
            ####################################################################
           
            # Build HTML address for image
            img_url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={}&dec={}&%22/pix={}&layer=ls-dr9&size={}'.format(ra, dec, pix_scale, img_size)
               
            self.files.append(img_url)

            # Extract galaxy name
            gal_name = self.objects[i]
               
            self.file_output[gal_name] = []
               
       
       
       
        self.curr_truth_display = []
       
        self.fig = plt.figure(figsize=(14,8.75))
                   
        self.display_axes = self.fig.add_axes([.29,.05,.7,.85], picker=True)
       
       
        self.next_axes = self.fig.add_axes([.08, .7, .1, .05])
       
        self.prev_axes = self.fig.add_axes([.08, .6, .1, .05])
       
        self.save_axes = self.fig.add_axes([.08, .5, .1, .05])


        #self.fiber_axes = self.fig.add_axes([0.08, 0.05, 0.2, 0.2])
       
       
       
        self.next_button = Button(self.next_axes, 'Next')
       
        self.next_button.on_clicked(self.next_button_func)
       
       
        self.prev_button = Button(self.prev_axes, "Prev")
       
        self.prev_button.on_clicked(self.prev_button_func)
       
       
        self.save_button = Button(self.save_axes, 'Save')
       
        self.save_button.on_clicked(self.save_button_func)
       
       
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
       
        self.seek_to_index(0)
       
        plt.show()
       
    def next_button_func(self, event):
       
        self.seek_to_index(self.curr_index+1)
       
    def prev_button_func(self, event):
       
        self.seek_to_index(self.curr_index-1)
       
    def save_button_func(self, event):
       
        outfile = open(self.out_filename, 'w')
       
        json.dump(self.file_output, outfile)
       
        outfile.close()
       
    def seek_to_index(self, index):
       
        self.curr_index = index
       
        if index < 0 or index >= len(self.files):
           
            print("No galaxies left!")
           
            self.seek_to_index(0)
       
       
       
       
        for artist, x_pix, y_pix in self.curr_truth_display:
           
            artist.remove()
           
        self.curr_truth_display = []
       
       
        galaxy_img_page = urlopen(self.files[index])
        galaxy_img_byte = io.BytesIO(galaxy_img_page.read())
        curr_frame = numpy.array(Image.open(galaxy_img_byte))
       
       
        self.display_axes.clear()
       
        self.display_axes.imshow(curr_frame, interpolation='nearest')

        self.display_axes.set_title(self.objects[index])

        ########################################################################
        # Plot the SGA ellipse footprint
        #-----------------------------------------------------------------------
        major_axis = self.input_table['DIAM'][index]
        axis_ratio = self.input_table['BA'][index]
        phi = self.input_table['PA'][index]

        # Convert major axis units from arcminutes to pixels
        major_axis_pixels = major_axis/pix_scale_arcmin

        img_size = curr_frame.shape

        center_row = 0.5*img_size[0]
        center_col = 0.5*img_size[1]

        SGA_ellipse = Ellipse((center_row, center_col), 
                              major_axis_pixels, 
                              major_axis_pixels*axis_ratio, 
                              angle=90 - phi, 
                              color='#03A9FC', 
                              fill=False)

        self.display_axes.add_artist(SGA_ellipse)
        ########################################################################

        '''
        ########################################################################
        # Plot the fiber patrol radius
        #-----------------------------------------------------------------------
        patrol_circle = plt.Circle((0, 0), 
                                   patrol_radius_pixels, 
                                   color='#000000', 
                                   fill=False)

        self.fiber_axes.add_artist(patrol_circle)
        ########################################################################
        '''
   
        self.add_existing_truth()
   
       
        plt.draw()
    


    def add_existing_truth(self):
       
       
        truth_data = self.file_output[self.objects[self.curr_index]]
       
        for y_pixel, x_pixel in truth_data:
           
            #new_circle = plt.Circle((x_pixel, y_pixel), 2.0, color='#00CC00', fill=False)
            new_circle = plt.Circle((x_pixel, y_pixel), 
                                    fiber_diameter_pixels, 
                                    color='#CC0000', 
                                    fill=False)
           
            self.curr_truth_display.append((new_circle, x_pixel, y_pixel))
           
            self.display_axes.add_artist(new_circle)
           
               
               
               
   
    def onpick(self, event):
       
        if event.artist == self.display_axes:
           
            x_pixel = event.mouseevent.xdata
           
            y_pixel = event.mouseevent.ydata
           
            mouse_button_pressed = event.mouseevent.button
           
            doubleclick = event.mouseevent.dblclick
           
            #remove a truth value
            if mouse_button_pressed == 3 or doubleclick:
               
                for idx, (curr_circle, x_pix, y_pix) in enumerate(self.curr_truth_display):
                   
                    if (x_pix - x_pixel)**2 + (y_pix - y_pixel)**2 < 0.5*fiber_diameter_pixels:
                       
                        curr_circle.remove()
                       
                       
                       
                        del self.curr_truth_display[idx]
                       
                       
                        del self.file_output[self.objects[self.curr_index]][idx]
                       
                        print("Removing object")
                       
                        #self.display_axes.draw(None)
                       
                        break
               
            #add a truth value
            else:
               
                new_circle = plt.Circle((x_pixel, y_pixel), 
                                        fiber_diameter_pixels, 
                                        color='#CC0000', 
                                        fill=False)
               
                self.curr_truth_display.append((new_circle, x_pixel, y_pixel))
               
                self.display_axes.add_artist(new_circle)
               
                #self.curr_truth_coords.append((y_pixel, x_pixel)) #row, col format
                self.file_output[self.objects[self.curr_index]].append((y_pixel, x_pixel))
               
            plt.draw()
               
               
               
               
if __name__ == "__main__":
   
    # Change this file name
    input_filename = "SGA_large_galaxies_0.fits"
   
    GalaxyMarker(input_filename)