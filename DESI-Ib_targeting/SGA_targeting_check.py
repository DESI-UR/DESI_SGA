import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.widgets import TextBox, Button
from PIL import Image
import json
import io
from urllib.request import urlopen

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

# Custom functions
from SGA_targeting import radec_to_xy



pix_scale = 0.25 # arcsec/pixel
pix_scale_arcmin = pix_scale/60
pix_scale_degree = pix_scale_arcmin/60

fiber_diameter = 1.52 # arcsec
fiber_diameter_pixels = fiber_diameter/pix_scale

patrol_radius = 1.4 # arcmin
patrol_radius_pixels = patrol_radius/pix_scale_arcmin



class GalaxyChecker(object):
   
    def __init__(self, SGA_filename):
       
        self.out_filename = "target_files/SGA2025_off-axis_targets_cleaned.txt"
        # self.out_filename = "target_files/test.txt"

        self.input_table = Table.read(SGA_filename, format='fits')
       
        # Get all the objects in the source file
        self.objects = self.input_table['SGAID']

        # Read in current target list
        # targets_filename = 'target_files/SGA2025_off-axis_targets_old.txt'
        targets_filename = 'target_files/SGA2025_off-axis_targets_cleaned.txt'
        infile = open(targets_filename, 'r')
        self.curr_target_list = json.load(infile)
        infile.close()

        # print([type(k) for k in self.curr_target_list.keys()])
        '''
        # Read in earlier cleaned target file
        # (This file was generated using all the large SGA-2025 galaxies, not 
        # just those that are not part of the SGA-2020.  We're just going to use 
        # entries in here if the galaxy doesn't already have targets selected.)
        incomplete_targets_filename = 'target_files/SGA2025_off-axis_targets_cleaned-ALL_do-not-use.txt'
        infile = open(incomplete_targets_filename, 'r')
        incomplete_target_list = json.load(infile)
        infile.close()
        '''
        # Initialize list of image urls
        self.files = []
       
        for i in range(len(self.input_table)):

            # Extract sky coordinates for object
            ra = self.input_table['RA'][i]
            dec = self.input_table['DEC'][i]

            ####################################################################
            # Determine size of image needed
            #-------------------------------------------------------------------
            major_axis = self.input_table['D26'][i]

            major_axis_pixels = major_axis/pix_scale_arcmin

            img_size = int(major_axis_pixels + 100)
            ####################################################################
           
            # Build HTML address for image
            img_url = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={}&dec={}&%22/pix={}&layer=ls-dr11&size={}'.format(ra, dec, pix_scale, img_size)
            
            self.files.append(img_url)
            '''
            # Extract galaxy name
            gal_name = self.objects[i]

            # Check to see if galaxy is already in target list; if not, add it
            if str(gal_name) not in self.curr_target_list.keys():
                # print(gal_name, "not yet targeted")

                # Check to see if I already "cleaned" the galaxy's targets
                if str(gal_name) in incomplete_target_list.keys():
                    # print("    I targeted", gal_name, "earlier:", len(incomplete_target_list[str(int(gal_name))]))

                    self.curr_target_list[str(int(gal_name))] = incomplete_target_list[str(int(gal_name))]

                else:

                    self.curr_target_list[str(int(gal_name))] = []
            '''
        self.curr_truth_display = []
       
        self.fig = plt.figure(figsize=(14,8.75))
                   
        self.display_axes = self.fig.add_axes([.29,.05,.7,.85], picker=True)
       
       
        self.next_axes = self.fig.add_axes([.08, .7, .1, .05])
       
        self.prev_axes = self.fig.add_axes([.08, .6, .1, .05])
       
        self.save_axes = self.fig.add_axes([.08, .8, .1, .05])

        self.clear_axes = self.fig.add_axes([.08, .1, .07, .05])

        self.zoom_axes = self.fig.add_axes([.08, .4, .1, .05])

        self.reset_axes = self.fig.add_axes([.08, .3, .1, .05])
       
       
        self.next_button = Button(self.next_axes, 'Next')
       
        self.next_button.on_clicked(self.next_button_func)
       
       
        self.prev_button = Button(self.prev_axes, "Prev")
       
        self.prev_button.on_clicked(self.prev_button_func)
       
       
        self.save_button = Button(self.save_axes, 'Save')
       
        self.save_button.on_clicked(self.save_button_func)

        self.clear_button = Button(self.clear_axes, 'ClearAll')
        self.clear_button.on_clicked(self.clear_existing_truth)

        self.zoom_button = Button(self.zoom_axes, 'TopLeft')
        self.zoom_button.on_clicked(self.topleft_zoom)

        self.reset_button = Button(self.reset_axes, 'ResetView')
        self.reset_button.on_clicked(self.reset_view)
       
       
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
       
        self.seek_to_index(0) # Change this value to start at a specific galaxy
       
        plt.show()
       
    def next_button_func(self, event):
       
        self.seek_to_index(self.curr_index+1)
       
    def prev_button_func(self, event):
       
        self.seek_to_index(self.curr_index-1)
       
    def save_button_func(self, event):
       
        outfile = open(self.out_filename, 'w')
       
        json.dump(self.curr_target_list, outfile)
       
        outfile.close()

    def topleft_zoom(self, event):

        # print('Clicked TopLeft')

        self.display_axes.set_xlim(0, 600)

        self.display_axes.set_ylim(600, 0)

        plt.draw()

    def reset_view(self, event):

        rows = self.curr_frame.shape[0]
        cols = self.curr_frame.shape[1]

        self.display_axes.set_xlim(0, cols)
        self.display_axes.set_ylim(rows, 0)

        plt.draw()
       
    def seek_to_index(self, index):
       
        self.curr_index = index
       
        if index < 0 or index >= len(self.files):
           
            print("No galaxies left!")
           
            self.seek_to_index(0)

            return
       
        ########################################################################
        # Only target galaxies at dec > -35 degrees
        #-----------------------------------------------------------------------
        curr_dec = self.input_table['DEC'][self.curr_index]

        if curr_dec < -35:

            self.seek_to_index(self.curr_index + 1)

            return
        ########################################################################
       
        '''
        for artist, x_pix, y_pix in self.curr_truth_display:
           
            artist.remove()
        '''
        self.curr_truth_display = []
       
        galaxy_img_page = urlopen(self.files[index])
        galaxy_img_byte = io.BytesIO(galaxy_img_page.read())
        curr_frame = np.array(Image.open(galaxy_img_byte))
       
        self.display_axes.clear()
       
        self.display_axes.imshow(curr_frame, interpolation='nearest')

        self.display_axes.set_title(str(index) + ' - ' + str(self.objects[index]) + ' ({:.3f}, {:.3f})'.format(self.input_table['RA'][self.curr_index], self.input_table['DEC'][self.curr_index]))

        self.curr_frame = curr_frame

        ########################################################################
        # Plot the SGA ellipse footprint
        #-----------------------------------------------------------------------
        major_axis = self.input_table['D26'][index]
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


        ########################################################################
        # Plot the higher-priority TF targets
        #-----------------------------------------------------------------------
        # Extract sky coordinates for object
        ra = self.input_table['RA'][index]
        dec = self.input_table['DEC'][index]
        center_sky = SkyCoord(ra*u.deg, dec*u.deg)

        #-----------------------------------------------------------------------
        # Center fiber
        #-----------------------------------------------------------------------
        center_fiber = plt.Circle((center_row, center_col), 
                                  fiber_diameter_pixels, 
                                  color='#ff80ff', 
                                  fill=False)
        self.display_axes.add_artist(center_fiber)
        #-----------------------------------------------------------------------


        #-----------------------------------------------------------------------
        # Minor axis fibers
        #-----------------------------------------------------------------------
        delta_b = 0.4*(0.5*major_axis*u.arcmin)*axis_ratio

        fiber1 = center_sky.directional_offset_by((phi + 90)*u.deg, delta_b)
        fiber2 = center_sky.directional_offset_by((phi - 90)*u.deg, delta_b)

        for fiber in [fiber1, fiber2]:
            fiber_xy = radec_to_xy(fiber.ra.value, fiber.dec.value, 
                                   ra0=ra, dec0=dec, 
                                   x0=center_row, y0=center_col, 
                                   xscale=pix_scale_degree, 
                                   yscale=pix_scale_degree)
            fiber_circ = plt.Circle(fiber_xy, 
                                    fiber_diameter_pixels, 
                                    color='#ff80ff', 
                                    fill=False)
            self.display_axes.add_artist(fiber_circ)
        #-----------------------------------------------------------------------


        #-----------------------------------------------------------------------
        # Major axis fibers
        #-----------------------------------------------------------------------
        x = np.arange(0.2,1.2,0.2)

        # Distances along the semi-major axis from the center coordinate
        delta_a = 0.5*(major_axis*u.arcmin)*x

        # Target positions
        fiber1 = center_sky.directional_offset_by(phi*u.deg, delta_a)
        fiber2 = center_sky.directional_offset_by((phi + 180)*u.deg, delta_a)

        for fiber in [fiber1, fiber2]:
            for i in range(len(fiber)):
                fiber_xy = radec_to_xy(fiber[i].ra.value, fiber[i].dec.value, 
                                       ra0=ra, dec0=dec, 
                                       x0=center_row, y0=center_col, 
                                       xscale=pix_scale_degree, 
                                       yscale=pix_scale_degree)
                fiber_circ = plt.Circle(fiber_xy, 
                                        fiber_diameter_pixels, 
                                        color='#ff80ff', 
                                        fill=False)
                self.display_axes.add_artist(fiber_circ)
        #-----------------------------------------------------------------------
        ########################################################################
   
        self.add_existing_truth()
   
       
        plt.draw()
    


    def add_existing_truth(self):
       
        truth_data = self.curr_target_list[str(self.objects[self.curr_index])]
        # print(self.curr_index, self.objects[self.curr_index], len(truth_data))
       
        for y_pixel, x_pixel in truth_data:
           
            #new_circle = plt.Circle((x_pixel, y_pixel), 2.0, color='#00CC00', fill=False)
            new_circle = plt.Circle((x_pixel, y_pixel), 
                                    fiber_diameter_pixels, 
                                    color='#CC0000', 
                                    # color='#bbbbbb',
                                    # edgegapcolor='#CC0000',
                                    # linestyle='dashed',
                                    fill=False)
           
            self.curr_truth_display.append((new_circle, x_pixel, y_pixel))
           
            self.display_axes.add_artist(new_circle)
           
               
    def clear_existing_truth(self, event):

        # print('Clicked ClearAll')

        for idx, (curr_circle, x_pix, y_pix) in enumerate(self.curr_truth_display):
   
            curr_circle.remove()

        self.curr_truth_display = []

        self.curr_target_list[str(self.objects[self.curr_index])] = []

        plt.draw()
               
               
   
    def onpick(self, event):
       
        if event.artist == self.display_axes:
           
            x_pixel = event.mouseevent.xdata
           
            y_pixel = event.mouseevent.ydata
           
            mouse_button_pressed = event.mouseevent.button

            #print(x_pixel, y_pixel)
           
            doubleclick = event.mouseevent.dblclick
            
            # add a truth value
            if mouse_button_pressed == 3 or doubleclick:
               
                new_circle = plt.Circle((x_pixel, y_pixel), 
                                        fiber_diameter_pixels, 
                                        color='#CC0000', 
                                        # color='#bbbbbb',
                                        fill=False)
               
                self.curr_truth_display.append((new_circle, x_pixel, y_pixel))
               
                self.display_axes.add_artist(new_circle)
               
                #self.curr_truth_coords.append((y_pixel, x_pixel)) #row, col format
                self.curr_target_list[str(self.objects[self.curr_index])].append((y_pixel, x_pixel))

            # remove a truth value
            else:

                for idx, (curr_circle, x_pix, y_pix) in enumerate(self.curr_truth_display):
               
                    if (x_pix - x_pixel)**2 + (y_pix - y_pixel)**2 < fiber_diameter_pixels**2:
                       
                        curr_circle.remove()
                       
                        del self.curr_truth_display[idx]
                       
                        del self.curr_target_list[str(self.objects[self.curr_index])][idx]
                       
                        #print("Removing object")
                       
                        #self.display_axes.draw(None)
                       
                        break
            
            plt.draw()
               
               
               
               
if __name__ == "__main__":
   
    # Change this file name
    target_file = 'SGA2025_large_galaxies_new.fits'
   
    GalaxyChecker(target_file)


