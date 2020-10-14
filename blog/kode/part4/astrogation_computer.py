# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from PIL import Image
import logging
import sys
import os

class astrogation_computer:
    def __init__(self, seed=33382, data_path='himmelkule_data',image='himmelkule.npy', projection_file='projections.npy'):

        # Configures logging
        self.logger = logging.getLogger('astrogation_computer')
        logging.basicConfig(
            filename='astrogation_computer.log',
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s : %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
        )

        self.system = SolarSystem(seed)
        self.mission = SpaceMission(seed)
        # Holds data for the spherical sky
        self.sky_image = np.load(image)
        self.logger.info(f'Sky image file {image} loaded with shape {self.sky_image.shape}')

        if (os.path.exists(f'{data_path}/{projection_file}')):
            # Holds data for the stereographic projection of the sky
            self.projections = np.load(f'{data_path}/{projection_file}')
            self.logger.info(f'Projection file {data_path}/{projection_file} loaded with shape {self.projections.shape}')

        self.log_arrays = False
        self.show_plots = False

    def field_of_view(self, theta, alpha):
        '''
        Caclulates the x_min/max
        and the y_min/max
        '''
        self.logger.info(f'Parameters theta={theta}, alpha={alpha}')
        x_numerator = 2*np.sin(alpha/2)
        x_denomerator = 1 + np.cos(alpha/2)

        y_numerator = 2*np.sin(theta/2)
        y_denomerator = 1 + np.cos(theta/2)

        x_max = x_numerator/x_denomerator
        x_min = -x_max

        y_max = y_numerator/y_denomerator
        y_min = -y_max
        self.logger.info(f'Field of view is X = [{x_min},{x_max}]\nY = [{y_min},{y_max}]')
        return x_min, x_max, y_min, y_max

    def k(self, theta, theta0, phi, phi0):
        return 2/(1+np.cos(theta0)*np.cos(theta)+np.sin(theta0)*np.sin(theta)*np.cos(phi-phi0))

    def transform_thetaphi_to_xy(self, theta, phi, theta0, phi0):
        '''
        Returns x,y transformed from theta and phi with center theta0, phi0
        '''
        k = self.k(theta,theta0,phi,phi0)
        x = k*np.sin(theta)*np.sin(phi-phi0)
        y = k*(np.sin(theta0)*np.cos(theta)-np.cos(theta0)*np.sin(theta)*np.cos(phi-phi0))
        return x,y

    def transform_xy_to_thetaphi(self,x,y,theta0,phi0):
        '''
        Returns theta, phi transformed from x,y with center theta0, phi0
        '''
        p = np.sqrt(x**2+y**2)
        B = 2*np.arctan(p/2)
        self.logger.debug(f'p:\n{p}')

        theta = np.pi/2 - np.arcsin(np.cos(B)*np.cos(theta0)+y*np.sin(B)*np.sin(theta0)/p)

        phi = phi0 + np.arctan(x*np.sin(B)/(p*np.sin(theta0)*np.cos(B)-y*np.cos(theta0)*np.sin(B)))

        return theta, phi

    def generate_sky_image(self,width,height,phi=0,theta=90,alpha=70):
        '''
        Generates an image from the sky centered around phi, theta angles in degrees with FOV alpha x alpha
        The returned array is of sixe height x width x 3
        '''
        self.logger.info(f'Generating sky image from coordinates ({phi},{theta}) with field of view {alpha:g}x{alpha:g}')

        theta0 = np.deg2rad(theta)
        phi0 = np.deg2rad(phi)
        phi_min = np.deg2rad(phi+alpha/2)
        phi_max = np.deg2rad(phi-alpha/2)
        theta_min = np.deg2rad(theta-alpha/2)
        theta_max = np.deg2rad(theta+alpha/2)

        self.logger.info(f'Image size is height x width {height} x {width}') # width = 480, height = 640
        self.logger.info(f'phi is in range [{phi_min},{phi_max}] radians')
        self.logger.info(f'phi is in range [{np.rad2deg(phi_min)},{np.rad2deg(phi_max)}] degrees')
        self.logger.info(f'theta is in range [{theta_min},{theta_max}] radians')
        self.logger.info(f'theta is in range [{np.rad2deg(theta_min)},{np.rad2deg(theta_max)}] degrees')

        xy_image_map = np.zeros((height,width,2))               # Maps pixel to x,y coordinate
        rad_image_map = np.zeros_like(xy_image_map)             # Maps pixel to theta, phi coordinate
        new_sky_image = np.zeros((height,width,3))              # Holds the sky partial image
        self.logger.info(f'new_sky_image shape {new_sky_image.shape}')

        # Finds the maximum and minimum x and y values for the given FOV
        x_min, x_max, y_min, y_max = self.field_of_view(np.deg2rad(alpha), np.deg2rad(alpha))

        # Constructs Arrays to hold x,y coordinates in the picture
        x_coor = np.linspace(x_min,x_max,width)
        y_coor = np.linspace(y_max,y_min,height)
        x, y = np.meshgrid(x_coor,y_coor)
        xy_image_map[:,:,0] = y
        xy_image_map[:,:,1] = x

        if(self.show_plots):
            plt.subplot(421)
            plt.title('XY Grid')
            plt.xlabel('X')
            plt.ylabel('Y')
            if(width >= 480):
                plt.scatter(x[0:width*height:40],y[0:width*height:40],marker='.')
            else:
                plt.scatter(x,y,marker='.')


        # Gets the theta,phi coordinates for x,y values in xy_image_map
        theta_coord, phi_coord = self.transform_xy_to_thetaphi(x,y,theta0, phi0)

        if(self.log_arrays):
            self.logger.info(f'Theta {theta_coord[0]},{theta_coord[height-1]}')
            self.logger.info(f'Phi {phi_coord[0]},{phi_coord[width-1]}')
            self.logger.info(f'Theta {np.rad2deg(theta_coord[0])},{np.rad2deg(theta_coord[height-1])}')
            self.logger.info(f'Phi {np.rad2deg(phi_coord[0])},{np.rad2deg(phi_coord[width-1])}')

        if(self.log_arrays):
            self.logger.info(f'phi_coord: {phi_coord}')
            self.logger.info(f'theta_coord {theta_coord}')

        self.logger.info(f'Creating theta phi meshgrid {theta_coord.shape} {phi_coord.shape}')
        rad_image_map[:,:,0],rad_image_map[:,:,1] = theta_coord, phi_coord

        if(self.show_plots):
            plt.subplot(422)
            plt.ylabel('Theta radians')
            plt.xlabel('Phi radians')
            plt.title('Theta Phi Grid')
            if(width >= 480):
                plt.scatter(rad_image_map[0:height:40,0:width:40,1],rad_image_map[0:height:40,0:width:40,0],c=['#ff7f0e'],marker='.')
            else:
                plt.scatter(rad_image_map[:,:,0],rad_image_map[:,:,1],c=['#ff7f0e'],marker='.')
        if(self.log_arrays):
            self.logger.info(f'rad_image_map: {np.rad2deg(rad_image_map)}')

        self.logger.info(f'Getting pixels')

        # Loop that generetas the stereographic projection
        for i in range(width):
            for j in range(height):
                t,p = rad_image_map[j,i,0], rad_image_map[j,i,1]
                 # Finds the index of the pixel we are looking at
                index = self.mission.get_sky_image_pixel(t,p)
                # Store the pixel
                pixel = self.sky_image[index,2:]
                new_sky_image[j,i] = pixel

        self.logger.info(f'new_sky_image.shape {new_sky_image.shape}')
        if(self.log_arrays):
            self.logger.info(f'new_sky_image {new_sky_image}')
        if(self.show_plots):
            plt.subplot(4,2,(5,8))
            plt.title(f'Generated picture FOV={alpha:g} rad.')

            plt.imshow(new_sky_image.astype('uint8'))
            plt.tight_layout()
            plt.show()

        return new_sky_image

    def find_orientation_angle(self,filename):
        '''
        Finds the best match for the angle the input file is oriented towards
        '''
        img = Image.open(filename)
        input = np.array(img.convert('RGB'))

        angle = 0
        # Set initial diff to infinity
        min_diff = float('inf')
        min_angle = angle
        print(f'Comparing to angle ', end='')
        for projection in self.projections:

            print(f'{angle:4d}\b\b\b\b', end = '',flush = True)
            diff = self.mean_square_diff(input,projection)
            if(diff < min_diff):
                min_diff = diff
                min_angle = angle
            angle += 1
        print('')
        self.logger.info(f'Found least difference in angle {min_angle}')

        return min_angle

    def save_array_as_png(self,img_array,filename):
        image = Image.fromarray(img_array.astype('uint8'), 'RGB')
        #image_rotated = image.rotate(90)
        image.save(filename)

    def sample(self, image='sample0000.png'):
        img = Image.open(image)
        self.logger.info(f'PNG info: {img.info}')
        pixels = np.array(img.convert('RGB'))
        self.logger.info(f'pixels array has shape {pixels.shape}')
        width = len(pixels[0, :])
        height = len(pixels[:, 0])
        #plt.imshow(pixels)
        #plt.show()
        self.logger.info(f'The sample image {image} is {height} by {width} pixels')
        self.reference_image = pixels
        return width, height

    def mean_square_diff(self,a, b):
    	err = np.sum((a.astype("float") - b.astype("float")) ** 2)
    	err /= float(a.shape[0] * b.shape[1])

    	return err
if __name__ == '__main__':
    data_path = 'himmelkule_data'

    astro_comp = astrogation_computer()

    # Find height and widt of the reference picture
    width, height = astro_comp.sample()

    # Generate a picture that should be the same as the reference and save as PNG
    generated_image = astro_comp.generate_sky_image(width,height,phi=0,theta=90,alpha=70)
    astro_comp.save_array_as_png(generated_image,'sample-generated-009070.png')

    # Generate stereographic projections of the entire sky
    angles = np.linspace(0,359,360)

    # Store projections in array
    projected_images = np.zeros((len(angles), height, width, 3))

    if(os.path.exists(data_path) == False):
        os.mkdir(data_path)

    for angle in angles:
        print(f'Generating for angle {angle}')
        generated_image = astro_comp.generate_sky_image(width,height,phi=angle,theta=90,alpha=70)
        projected_images[int(angle)] = generated_image
        img = Image.fromarray(generated_image.astype('uint8'))
        img.save(f'{data_path}/himmelkule_{int(angle)}_degs.png')

    np.save(f'{data_path}/projections.npy',projected_images)

'''
Når man kjører denne kommandoen genereres himmelkule_data/projections.npy som
brukes når man kjører manual_orientation.py
Det lages også 360 png bilder, en for hver grad i himmelkule_data/

Generer også sample-generated-009070.png som er bildet ved phi=0 og theta=90 grader og FOV 70x70

$> python astrogation_computer.py
Generating for angle 1.0
Generating for angle 2.0
Generating for angle 3.0
...
...
Generating for angle 358.0
Generating for angle 359.0
'''
