# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from PIL import Image
import logging
import os


class astrogation_computer:
    def __init__(self, seed=33382, data_path='himmelkule_data',
                 image='himmelkule.npy', projection_file='projections.npy'):

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
        k = self.k(theta, theta0, phi, phi0)
        x = k*np.sin(theta)*np.sin(phi-phi0)
        y = k*(np.sin(theta0)*np.cos(theta)-np.cos(theta0)*np.sin(theta)*np.cos(phi-phi0))
        return x, y

    def transform_xy_to_thetaphi(self, x, y, theta0, phi0):
        '''
        Returns theta, phi transformed from x,y with center theta0, phi0
        '''
        p = np.sqrt(x**2+y**2)
        B = 2*np.arctan(p/2)
        self.logger.debug(f'p:\n{p}')

        theta = np.pi/2 - np.arcsin(np.cos(B)*np.cos(theta0)+y*np.sin(B)*np.sin(theta0)/p)

        phi = phi0 + np.arctan(x*np.sin(B)/(p*np.sin(theta0)*np.cos(B)-y*np.cos(theta0)*np.sin(B)))

        return theta, phi

    def generate_sky_image(self, width, height, phi=0, theta=90, alpha=70):
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

        xy_image_map = np.zeros((height, width, 2))               # Maps pixel to x,y coordinate
        rad_image_map = np.zeros_like(xy_image_map)             # Maps pixel to theta, phi coordinate
        new_sky_image = np.zeros((height, width, 3))              # Holds the sky partial image
        self.logger.info(f'new_sky_image shape {new_sky_image.shape}')

        # Finds the maximum and minimum x and y values for the given FOV
        x_min, x_max, y_min, y_max = self.field_of_view(np.deg2rad(alpha), np.deg2rad(alpha))

        # Constructs Arrays to hold x,y coordinates in the picture
        x_coor = np.linspace(x_min, x_max, width)
        y_coor = np.linspace(y_max, y_min, height)
        x, y = np.meshgrid(x_coor, y_coor)
        xy_image_map[:, :, 0] = y
        xy_image_map[:, :, 1] = x

        if(self.show_plots):
            plt.subplot(421)
            plt.title('XY Grid')
            plt.xlabel('X')
            plt.ylabel('Y')
            if(width >= 480):
                plt.scatter(x[0:width*height:40], y[0:width*height:40], marker='.')
            else:
                plt.scatter(x, y, marker='.')


        # Gets the theta,phi coordinates for x,y values in xy_image_map
        theta_coord, phi_coord = self.transform_xy_to_thetaphi(x, y, theta0, phi0)

        if(self.log_arrays):
            self.logger.info(f'Theta {theta_coord[0]},{theta_coord[height-1]}')
            self.logger.info(f'Phi {phi_coord[0]},{phi_coord[width-1]}')
            self.logger.info(f'Theta {np.rad2deg(theta_coord[0])},{np.rad2deg(theta_coord[height-1])}')
            self.logger.info(f'Phi {np.rad2deg(phi_coord[0])},{np.rad2deg(phi_coord[width-1])}')

        if(self.log_arrays):
            self.logger.info(f'phi_coord: {phi_coord}')
            self.logger.info(f'theta_coord {theta_coord}')

        self.logger.info(f'Creating theta phi meshgrid {theta_coord.shape} {phi_coord.shape}')
        rad_image_map[:, :, 0],rad_image_map[:, :, 1] = theta_coord, phi_coord

        if(self.show_plots):
            plt.subplot(422)
            plt.ylabel('Theta radians')
            plt.xlabel('Phi radians')
            plt.title('Theta Phi Grid')
            if(width >= 480):
                plt.scatter(rad_image_map[0:height:40, 0:width:40,1],
                            rad_image_map[0:height:40,0:width:40,0],
                            c=['#ff7f0e'], marker='.')
            else:
                plt.scatter(rad_image_map[:,:,0],rad_image_map[:,:,1],c=['#ff7f0e'],marker='.')
        if(self.log_arrays):
            self.logger.info(f'rad_image_map: {np.rad2deg(rad_image_map)}')

        self.logger.info(f' Getting pixels')

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
$ python astrogation_computer.py 
Generating for angle 0.0
Generating for angle 1.0
Generating for angle 2.0
Generating for angle 3.0
Generating for angle 4.0
Generating for angle 5.0
Generating for angle 6.0
Generating for angle 7.0
Generating for angle 8.0
Generating for angle 9.0
Generating for angle 10.0
Generating for angle 11.0
Generating for angle 12.0
Generating for angle 13.0
Generating for angle 14.0
Generating for angle 15.0
Generating for angle 16.0
Generating for angle 17.0
Generating for angle 18.0
Generating for angle 19.0
Generating for angle 20.0
Generating for angle 21.0
Generating for angle 22.0
Generating for angle 23.0
Generating for angle 24.0
Generating for angle 25.0
Generating for angle 26.0
Generating for angle 27.0
Generating for angle 28.0
Generating for angle 29.0
Generating for angle 30.0
Generating for angle 31.0
Generating for angle 32.0
Generating for angle 33.0
Generating for angle 34.0
Generating for angle 35.0
Generating for angle 36.0
Generating for angle 37.0
Generating for angle 38.0
Generating for angle 39.0
Generating for angle 40.0
Generating for angle 41.0
Generating for angle 42.0
Generating for angle 43.0
Generating for angle 44.0
Generating for angle 45.0
Generating for angle 46.0
Generating for angle 47.0
Generating for angle 48.0
Generating for angle 49.0
Generating for angle 50.0
Generating for angle 51.0
Generating for angle 52.0
Generating for angle 53.0
Generating for angle 54.0
Generating for angle 55.0
Generating for angle 56.0
Generating for angle 57.0
Generating for angle 58.0
Generating for angle 59.0
Generating for angle 60.0
Generating for angle 61.0
Generating for angle 62.0
Generating for angle 63.0
Generating for angle 64.0
Generating for angle 65.0
Generating for angle 66.0
Generating for angle 67.0
Generating for angle 68.0
Generating for angle 69.0
Generating for angle 70.0
Generating for angle 71.0
Generating for angle 72.0
Generating for angle 73.0
Generating for angle 74.0
Generating for angle 75.0
Generating for angle 76.0
Generating for angle 77.0
Generating for angle 78.0
Generating for angle 79.0
Generating for angle 80.0
Generating for angle 81.0
Generating for angle 82.0
Generating for angle 83.0
Generating for angle 84.0
Generating for angle 85.0
Generating for angle 86.0
Generating for angle 87.0
Generating for angle 88.0
Generating for angle 89.0
Generating for angle 90.0
Generating for angle 91.0
Generating for angle 92.0
Generating for angle 93.0
Generating for angle 94.0
Generating for angle 95.0
Generating for angle 96.0
Generating for angle 97.0
Generating for angle 98.0
Generating for angle 99.0
Generating for angle 100.0
Generating for angle 101.0
Generating for angle 102.0
Generating for angle 103.0
Generating for angle 104.0
Generating for angle 105.0
Generating for angle 106.0
Generating for angle 107.0
Generating for angle 108.0
Generating for angle 109.0
Generating for angle 110.0
Generating for angle 111.0
Generating for angle 112.0
Generating for angle 113.0
Generating for angle 114.0
Generating for angle 115.0
Generating for angle 116.0
Generating for angle 117.0
Generating for angle 118.0
Generating for angle 119.0
Generating for angle 120.0
Generating for angle 121.0
Generating for angle 122.0
Generating for angle 123.0
Generating for angle 124.0
Generating for angle 125.0
Generating for angle 126.0
Generating for angle 127.0
Generating for angle 128.0
Generating for angle 129.0
Generating for angle 130.0
Generating for angle 131.0
Generating for angle 132.0
Generating for angle 133.0
Generating for angle 134.0
Generating for angle 135.0
Generating for angle 136.0
Generating for angle 137.0
Generating for angle 138.0
Generating for angle 139.0
Generating for angle 140.0
Generating for angle 141.0
Generating for angle 142.0
Generating for angle 143.0
Generating for angle 144.0
Generating for angle 145.0
Generating for angle 146.0
Generating for angle 147.0
Generating for angle 148.0
Generating for angle 149.0
Generating for angle 150.0
Generating for angle 151.0
Generating for angle 152.0
Generating for angle 153.0
Generating for angle 154.0
Generating for angle 155.0
Generating for angle 156.0
Generating for angle 157.0
Generating for angle 158.0
Generating for angle 159.0
Generating for angle 160.0
Generating for angle 161.0
Generating for angle 162.0
Generating for angle 163.0
Generating for angle 164.0
Generating for angle 165.0
Generating for angle 166.0
Generating for angle 167.0
Generating for angle 168.0
Generating for angle 169.0
Generating for angle 170.0
Generating for angle 171.0
Generating for angle 172.0
Generating for angle 173.0
Generating for angle 174.0
Generating for angle 175.0
Generating for angle 176.0
Generating for angle 177.0
Generating for angle 178.0
Generating for angle 179.0
Generating for angle 180.0
Generating for angle 181.0
Generating for angle 182.0
Generating for angle 183.0
Generating for angle 184.0
Generating for angle 185.0
Generating for angle 186.0
Generating for angle 187.0
Generating for angle 188.0
Generating for angle 189.0
Generating for angle 190.0
Generating for angle 191.0
Generating for angle 192.0
Generating for angle 193.0
Generating for angle 194.0
Generating for angle 195.0
Generating for angle 196.0
Generating for angle 197.0
Generating for angle 198.0
Generating for angle 199.0
Generating for angle 200.0
Generating for angle 201.0
Generating for angle 202.0
Generating for angle 203.0
Generating for angle 204.0
Generating for angle 205.0
Generating for angle 206.0
Generating for angle 207.0
Generating for angle 208.0
Generating for angle 209.0
Generating for angle 210.0
Generating for angle 211.0
Generating for angle 212.0
Generating for angle 213.0
Generating for angle 214.0
Generating for angle 215.0
Generating for angle 216.0
Generating for angle 217.0
Generating for angle 218.0
Generating for angle 219.0
Generating for angle 220.0
Generating for angle 221.0
Generating for angle 222.0
Generating for angle 223.0
Generating for angle 224.0
Generating for angle 225.0
Generating for angle 226.0
Generating for angle 227.0
Generating for angle 228.0
Generating for angle 229.0
Generating for angle 230.0
Generating for angle 231.0
Generating for angle 232.0
Generating for angle 233.0
Generating for angle 234.0
Generating for angle 235.0
Generating for angle 236.0
Generating for angle 237.0
Generating for angle 238.0
Generating for angle 239.0
Generating for angle 240.0
Generating for angle 241.0
Generating for angle 242.0
Generating for angle 243.0
Generating for angle 244.0
Generating for angle 245.0
Generating for angle 246.0
Generating for angle 247.0
Generating for angle 248.0
Generating for angle 249.0
Generating for angle 250.0
Generating for angle 251.0
Generating for angle 252.0
Generating for angle 253.0
Generating for angle 254.0
Generating for angle 255.0
Generating for angle 256.0
Generating for angle 257.0
Generating for angle 258.0
Generating for angle 259.0
Generating for angle 260.0
Generating for angle 261.0
Generating for angle 262.0
Generating for angle 263.0
Generating for angle 264.0
Generating for angle 265.0
Generating for angle 266.0
Generating for angle 267.0
Generating for angle 268.0
Generating for angle 269.0
Generating for angle 270.0
Generating for angle 271.0
Generating for angle 272.0
Generating for angle 273.0
Generating for angle 274.0
Generating for angle 275.0
Generating for angle 276.0
Generating for angle 277.0
Generating for angle 278.0
Generating for angle 279.0
Generating for angle 280.0
Generating for angle 281.0
Generating for angle 282.0
Generating for angle 283.0
Generating for angle 284.0
Generating for angle 285.0
Generating for angle 286.0
Generating for angle 287.0
Generating for angle 288.0
Generating for angle 289.0
Generating for angle 290.0
Generating for angle 291.0
Generating for angle 292.0
Generating for angle 293.0
Generating for angle 294.0
Generating for angle 295.0
Generating for angle 296.0
Generating for angle 297.0
Generating for angle 298.0
Generating for angle 299.0
Generating for angle 300.0
Generating for angle 301.0
Generating for angle 302.0
Generating for angle 303.0
Generating for angle 304.0
Generating for angle 305.0
Generating for angle 306.0
Generating for angle 307.0
Generating for angle 308.0
Generating for angle 309.0
Generating for angle 310.0
Generating for angle 311.0
Generating for angle 312.0
Generating for angle 313.0
Generating for angle 314.0
Generating for angle 315.0
Generating for angle 316.0
Generating for angle 317.0
Generating for angle 318.0
Generating for angle 319.0
Generating for angle 320.0
Generating for angle 321.0
Generating for angle 322.0
Generating for angle 323.0
Generating for angle 324.0
Generating for angle 325.0
Generating for angle 326.0
Generating for angle 327.0
Generating for angle 328.0
Generating for angle 329.0
Generating for angle 330.0
Generating for angle 331.0
Generating for angle 332.0
Generating for angle 333.0
Generating for angle 334.0
Generating for angle 335.0
Generating for angle 336.0
Generating for angle 337.0
Generating for angle 338.0
Generating for angle 339.0
Generating for angle 340.0
Generating for angle 341.0
Generating for angle 342.0
Generating for angle 343.0
Generating for angle 344.0
Generating for angle 345.0
Generating for angle 346.0
Generating for angle 347.0
Generating for angle 348.0
Generating for angle 349.0
Generating for angle 350.0
Generating for angle 351.0
Generating for angle 352.0
Generating for angle 353.0
Generating for angle 354.0
Generating for angle 355.0
Generating for angle 356.0
Generating for angle 357.0
Generating for angle 358.0
Generating for angle 359.0
'''
