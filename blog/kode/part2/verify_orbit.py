# Egen kode
#
# Skriver ut på CSV format resultater til Part 2 del B
# Kjør med: python verify_orbit.py -h for hjelp til parametre
#
import argparse
import load_orbit_sim as los
import orbit_sim
import numpy as np
import ast2000tools.constants as c
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-rs','--run-sim',action='store_true', help='Runs the simulation of orbits from scratch and saves the result')
parser.add_argument('-d','--download',action='store_true', help='Downloads pickle file')

args = parser.parse_args()

def period(a):
    '''
    Calculates period based on semi-major axis a
    '''
    mu = c.G_sol*(orbit.system.star_mass+orbit.system.masses)
    return 2*np.pi*np.sqrt(a**3/mu)

def period_kepler(a):
    '''
    Calculates period based on semi-major axis a with Keplers constant containg only the star mass
    '''
    mu = c.G_sol*orbit.system.star_mass
    return 2*np.pi*np.sqrt(a**3/mu)

def print_table(a,b,labels,title):
    '''
    Prints a table that can be piped to a csv file
    '''
    print(title)
    for i in range(len(labels)):
        print(labels[i],end=',')
    print('')
    for i in range(len(a)):
        print(f'{i+1},{a[i]},{b[i]},{error_percentage(a[i],b[i])}%')

def error_percentage(a,b):
    '''
    Calculates the difference between a and b as percentage
    '''
    diff = abs(a-b)
    return (diff/a)*100

if __name__ == '__main__':
    filename = f"simulated_orbits.pkl"
    orbit = los.orbit_sim_factory(filename,args)

    t=0
    #Find t for perihelion of hoth
    t=orbit.hoth_perhelion_t()
    t2=7000


    initialArea=orbit.area_swiped_during_period(50,2)
    print_table(initialArea,orbit.area_swiped_during_period(t,2),['Planet #','Area at t=0',f'Area at t={t}','Error'],'Area swiped in period 2 days')
    print_table(initialArea,orbit.area_swiped_during_period(t2,2),['Planet #','Area at t=0',f'Area at t={t2}','Error'],'Area swiped in period 2 days')
    print(f'Distance travled at aphelion during 2 days is {orbit.distance_traveled_during_period(0,2)[0]} AU')
    print(f'Distance travled at perhelion during 2 days is {orbit.distance_traveled_during_period(t,2)[0]} AU')
    print(f'Planet Hoth mean velocity is {orbit.mean_velocity()[0]:g} km/s')

    aphelium = orbit.simulated_aphelion()[:,0]
    print_table(orbit.axes,aphelium,['Planet #','Real aphelion','Simulated aphelion','Error'],'Aphelion comparison')
    print_table(period(orbit.axes),period(aphelium),['Planet #','Real period','Simulated period','Error'],'Period comparison')
    print_table(period(orbit.axes),period_kepler(orbit.axes),['Planet #','Newton period','Kepler period','Error'],'Kepler vs. Newton')
