# Egen kode
import numpy as np
import ast2000tools.constants as c
from ast2000tools.solar_system import SolarSystem
import orbit_sim


class habitable_zone:
    '''
    A class for finding info
    about planets in the
    habitable zone
    '''
    def __init__(self, seed=33382):
        self.system = SolarSystem(seed)                         # Our system
        self.radii = self.system.radii                          # All the radii in the system
        self.star_T = self.system.star_temperature              # Temprature of the star in Kelvin
        self.star_R = self.system.star_radius*1000              # Star radius in meters
        self.planet_pos = self.system.initial_positions*c.AU    # Planet positions in meters
        self.sorting()

    def sorting(self):
        N = len(self.radii)
        s = []
        for i in range(N):
            rd = np.linalg.norm(self.planet_pos[:, i])
            s.append((rd, i))
        s.sort()
        self.s = s
        self.planet_sorted = np.zeros((2, N), float)
        for i in range(N):
            j = s[i][1]
            self.planet_sorted[:, i] = self.planet_pos[:, j]

    def flux_boltz(self, T):
        '''
        Find the flux using boltzmann
        '''
        F = c.sigma*T**4
        return F

    def luminosity(self):
        '''
        Finds the luminosity
        '''
        T = self.star_T
        R = self.star_R
        self.L = 4*np.pi*R**2 * self.flux_boltz(T)

    def flux_radii(self, r):
        '''
        Find flux from radii
        '''
        F = self.L/(4*np.pi*r**2)
        return F

    def wattage_area(self, r, W=40):
        '''
        Caculates the Area, needed to achieve a certain watt
        '''
        efficiency = 0.12                               # efficiency of the panels

        R = self.star_R
        T = self.star_T

        A = (W * r**2)/(R**2*c.sigma*T**4*efficiency)   # Area needed
        return A

    def energy_at_planet(self, planet):
        '''
        At planets initial position
        '''
        r = self.radii[planet]
        R = self.star_R
        T = self.star_T
        rd = self.planet_sorted[:, planet]

        top = (2*np.pi*r**2*c.sigma*R**2*T**4)          # Numerator

        def dE(r):
            '''
            function for the dE
            '''
            ans = top/r**2
            return ans
        E = dE(rd)  # Bit unsure of the integration

        return E

    def temp_at_planet(self, planet):
        '''
        At planets initial position
        '''
        rd = self.planet_sorted[:, planet]
        rd_norm = np.linalg.norm(rd)        # Distance to planet
        R = self.star_R                     # km
        T = self.star_T                     # K

        Tp = np.sqrt(R/rd_norm)*T           # Temp at planet in K

        return Tp

    def temp_at_sim(self, planet):
        '''
        Simulated temperature all year round
        '''
        orbit = orbit_sim()                 # Not callable for some reason
        R = self.star_R

        mu = c.G_sol*(self.M * self.system.masses[planet])

        rotational_orbit_in_years = 2*np.pi*np.sqrt(self.axes**3/mu)
        T = 45*rotational_orbit_in_years[planet]

        dt = rotational_orbit_in_years[planet]/1000

        r0 = self.planet_sorted[:, planet]
        v0 = self.system.initial_velocities[:, planet]

        r, v, a, t = orbit.leapfrog(r0, v0, T, dt)

        temp_sim = np.sqrt(R/r)*self.star_T

        return temp_sim

    def temp_at_all_planets(self):
        '''
        Table of all planets temp at initial positions
        '''
        N = len(self.radii)
        R = self.star_R                     # km
        T = self.star_T                     # K
        print(f'Star radius is {R:g} m')
        print(f'Star temprature {T} K')

        for i in range(N):
            '''
            Prints out table of tempratures
            '''
            rd = self.planet_sorted[:, i]
            rd_norm = np.linalg.norm(rd)

            print(f'Distance to planet is {rd_norm/c.AU:g}')

            Tp = np.sqrt(R/rd_norm)*T

            print(f'Planet {i+1} has temprature {Tp:.2f} K at t=0')
            if 260 <= Tp and Tp <= 390:
                '''
                Checks if current planet is withing
                habitable zone
                '''
                print(f'Planet {i+1} is in the habitable zone at t=0')
                print(f'Distance at t=0 is {rd_norm/c.AU:g} AU')
                print(f'And is a type {self.system.types[self.s[i][1]]}')

    def just_calculations(self):
        '''
        Just the calculations no print
        of temp_at_all_planets
        '''
        N = len(self.radii)
        R = self.star_R                     # km
        T = self.star_T                     # K

        self.hab_plan = []

        for i in range(N):
            rd = self.planet_sorted[:, i]
            rd_norm = np.linalg.norm(rd)


            Tp = np.sqrt(R/rd_norm)*T

            if 245 <= Tp and Tp <= 405:
                '''
                Checks if current planet is withing
                habitable zone
                '''
                temp = (self.system.types[self.s[i][1]], i, rd_norm/c.AU, Tp)
                self.hab_plan.append(temp)

    def solar_panels(self, planet, W=40):
        '''
        Calculate the size the solar panels
        needs to be in order to give enough
        wattage at the planets position
        '''
        r = self.planet_sorted[:, planet]      # Pos of planet in m
        r_norm = np.linalg.norm(r)
        A = self.wattage_area(r_norm, W)    # Area in m**2
        return A

    def print_habitable_planet(self):
        '''
        prints only the habitable planets
        '''
        self.just_calculations()
        for i in range(len(self.hab_plan)):
            print(f'Planet {1+self.hab_plan[i][1]} from the star is in the habitable zone')
            print(f'With a temprature {self.hab_plan[i][3]:g} K. And is type {self.hab_plan[i][0]}')
            print(f'And would need a solar panel of size {self.solar_panels(self.hab_plan[i][1]):g} m^2')
            print(f'To operate at 40 W\n')

    def info_normal_solar_system(self):
        '''
        Just for comparisng henceforth no
        pretty prints or fancy functions
        '''
        planets_distance = np.array([57.9, 108, 150,
                                    228, 778, 1429,
                                    2871, 4504, 5914])
        planets_distance = planets_distance*(10**9)

        planet_radii = np.array([2.44, 6.05, 6.38, 3.4,
                                71.5, 60.3, 25.6, 24.1, 1.2])*10**6

        T_sun = 5780  # K
        R_sun = c.R_sun

        Tp = np.sqrt(R_sun/planets_distance)*T_sun
        for i in range(len(Tp)):
            print(f'{Tp[i]:g}')
            print(f'{planets_distance[i]/c.AU:g}')



if __name__ == '__main__':
    find_home = habitable_zone()
    find_home.print_habitable_planet()
    find_home.temp_at_all_planets()
    #find_home.info_normal_solar_system()