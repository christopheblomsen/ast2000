from PIL import Image
import numpy as np
import sys
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission

#LOADING THE DATA
sky_sphere = np.load('himmelkule.npy')

#FIELD OF VIEW INITIAL CONDITIONS (measured in degrees)
a_phi = 70
a_theta = 70
#theta direction for camera:
theta_0 = 90
#pixel size of flat picture:
npix_x = 640   #number of pixels in x-direction
npix_y = 480   #number of pixels in y-direction

#CONVERTING DEGREES TO RADIANS
a_phi = np.deg2rad(a_phi)
a_theta = np.deg2rad(a_theta)
theta_0 = np.deg2rad(theta_0)


#you will probably not need all of the following functions, which ones you need
#depends on the way you choose to solve the problem

#Find xy-limits of the flat picture:
def find_xy_limits(a_phi, a_theta):
    x_max = (2*np.sin(a_phi/2))/(1 + np.cos(a_phi/2))
    x_min = -(2*np.sin(a_phi/2))/(1 + np.cos(a_phi/2))
    y_max = (2*np.sin(a_theta/2))/(1 + np.cos(a_theta/2))
    y_min = -(2*np.sin(a_theta/2))/(1 + np.cos(a_theta/2))
    return x_min, x_max, y_min, y_max


#convert x and y meshgrid to theta/phi angles,  all angles in radians.
def xy2ang(X,Y,theta_0,phi_0):

    rho = np.sqrt(X**2 + Y**2)                  #Absolute Distance from Origin
    beta = 2*np.arctan(rho/2)                      #A Constant
    #Calculating the azimuthal angle of every given point
    theta = np.pi/2 - np.arcsin(np.cos(beta)*np.cos(theta_0) + \
                                Y*np.sin(beta)*np.sin(theta_0)/rho)
    #Calculating the angle of every given point
    phi = phi_0 + np.arctan(X*np.sin(beta)/(rho*np.sin(theta_0)*np.cos(beta) - \
                                            Y*np.cos(theta_0)*np.sin(beta)))
    return theta, phi

#convert theta,phi to x,y:  theta_0 and phi_0 are the centre of the flat picture
def ang2xy(theta_0, phi_0, theta, phi):

    kappa = 2./(1. + cos(theta_0)*cos(theta) + sin(theta_0)*sin(theta)*cos(phi - phi_0))
    X = kappa * sin(theta)*sin(phi - phi_0)
    Y = kappa * (sin(theta_0)*cos(theta) - cos(theta_0)*sin(theta)*cos(phi - phi_0))

    return X,Y

def save_array_as_png(img_array,filename):
    image = Image.fromarray(img_array.astype('uint8'), 'RGB')
    image.save(filename)

#start the projection here:
x_min, x_max, y_min, y_max = find_xy_limits(a_phi,a_theta)
print(f'x_min {x_min} x_max {x_max} y_min {y_min} y_max {y_max}')
#insert code
mission = SpaceMission(33382)

#CREATING A MESHGRID OF PIXEL COORDINATES (you may or may not need this)
x = np.linspace(x_max, x_min, npix_x)
y = np.linspace(y_max, y_min, npix_y)
Y, X = np.meshgrid(x, y)

#insert code
angles = np.linspace(0,359,360)

projections = np.zeros((len(angles),npix_y, npix_x, 3))
print(f'projections.shape {projections.shape}')
for angle in angles:
    t, p = xy2ang(X,Y,theta_0,np.deg2rad(angle))
    #T, P = np.meshgrid(t, p)
    TP = np.zeros((npix_y,npix_x,, 2))
    TP[:,:,0]= t
    TP[:,:,1]= p
    print(f'Generating for angle {angle}')
    for i in range(npix_x):
        for j in range(npix_y):
            index = mission.get_sky_image_pixel(TP[j,i,0],TP[j,i,1])
            pixel = sky_sphere[index,2:]
            #print(f'pixel [{int(angle)},{j},{i}] {pixel}')
            projections[int(angle),j,i] = pixel

    img = Image.fromarray(projections[int(angle)].astype('uint8'))
    img.save('himmelkule_{}_degs.png'.format(int(angle)))
    sys.exit()

#SAVING THE DATA TO A NUMPY DATA FILE (.npy using "np.save(name, array)")
np.save('projections.npy', projections)

#projections = np.load('projections.npy').astype('uint8')
#CREATING AN RGB IMAGE FROM GIVEN ANGLES (in degrees)

for angle in angles:
    img = Image.fromarray(projections[int(angle)].astype('uint8'))
    img.save('himmelkule_{}_degs.png'.format(int(angle)))
