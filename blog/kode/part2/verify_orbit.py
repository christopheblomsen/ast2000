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

    orbit.verify_planet_positions()

'''
Kjøreeksempel:
janmagneandersen$ python verify_orbit.py
Area swiped in period 2 days
Planet #,Area at t=0,Area at t=5000,Error,
1,0.008826228295737558,0.008826215422738894,0.0001458493733929629%
2,0.007028825221520661,0.006796570391173079,3.3043193283063923%
3,0.0030061582483555246,0.002851326090724846,5.150499236538107%
4,0.0037170862554012922,0.0035250225055988866,5.167051195632549%
5,0.004942520366062345,0.004731210707898278,4.275342184020489%
6,0.0024354458776144213,0.0023106532502059516,5.124015629150712%
7,0.00625232609984852,0.006123854414259173,2.054782228848553%
8,0.01352647154056166,0.011985442072296097,11.392693679534217%
Area swiped in period 2 days
Planet #,Area at t=0,Area at t=7000,Error,
1,0.008826228295737558,0.008828994132020674,0.031336559518310166%
2,0.007028825221520661,0.006934360443863687,1.343962535414092%
3,0.0030061582483555246,0.0029085701703304803,3.246272150790082%
4,0.0037170862554012922,0.0035936972428370103,3.319508994040305%
5,0.004942520366062345,0.0048198446198385704,2.4820483708296615%
6,0.0024354458776144213,0.0023573546065146058,3.2064465820241437%
7,0.00625232609984852,0.00624002420621241,0.1967570699232683%
8,0.01352647154056166,0.013033307071606201,3.6459210184756055%
Distance travled at aphelion during 2 days is 0.04825532995276159 AU
Distance travled at perhelion during 2 days is 0.05081290366399219 AU
Planet Hoth mean velocity is 9.05258 km/s
Aphelion comparison
Planet #,Real aphelion,Simulated aphelion,Error,
1,0.35660598633414764,0.3658122623745753,2.5816381085092583%
2,0.5888151773912673,0.5908024532354111,0.3375041813542244%
3,3.1696784965056795,3.2023777143212433,1.0316256948965692%
4,2.310491470065132,2.2092395862699,4.38226607226461%
5,1.4291331832748377,1.3036410541769385,8.780996100750793%
6,4.466365726836965,4.686758195563071,4.934492206982457%
7,0.784031429649257,0.76686985091582,2.1888891292424733%
8,0.17134082142599272,0.16555578172716445,3.37633475238532%
Period comparison
Planet #,Real period,Simulated period,Error,
1,0.24747008708693347,0.25711484732902473,3.8973438590592853%
2,0.5250604525576981,0.5277208456137986,0.5066831910765943%
3,6.557834998474799,6.659574736005569,1.5514226502257584%
4,4.081261882456952,3.815945237856431,6.500848321960623%
5,1.9854084827328173,1.7297281819937476,12.877969594807936%
6,10.959700063942728,11.78083482031989,7.492310479177119%
7,0.8067538431277292,0.7804109061001968,3.2653004695214927%
8,0.0824197994911378,0.0782810811513506,5.021509837854209%
Kepler vs. Newton
Planet #,Newton period,Kepler period,Error,
1,0.24747008708693347,0.24747070147434885,0.0002482673452050846%
2,0.5250604525576981,0.5250604895591076,7.04707606746583e-06%
3,6.557834998474799,6.5578830227564024,0.0007323191512908767%
4,4.081261882456952,4.081289478203675,0.0006761572160075316%
5,1.9854084827328173,1.9854086067928578,6.248590229254694e-06%
6,10.959700063942728,10.969132150007768,0.08606153462238725%
7,0.8067538431277292,0.8067538643056904,2.625083393518932e-06%
8,0.0824197994911378,0.0824199054955,0.0001286151663304607%
The biggest relative deviation was for planet 7, which drifted 0.139673 % from its actual position.
Your planet trajectories were satisfyingly calculated. Well done!
*** Achievement unlocked: Well-behaved planets! ***

'''
