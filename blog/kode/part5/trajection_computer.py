# Egen kode
import numpy as np
import matplotlib.pyplot as plt
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as c
import logging
import os
try:
    import cPickle as pickle
except:
    import pickle
from matplotlib.animation import FuncAnimation
import tools

class trajection_computer:

    def __init__(self, spacecraft_mass, seed=33382):

        # Configures logging
        self.log = logging.getLogger('trajection_computer')
        logging.basicConfig(
            filename='trajection_computer.log',
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s : %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
        )

        self.spacecraft_mass = spacecraft_mass
        # Downloads planet trajectories from URL if they do not exist on local disk
        trajectories_filename = "planet_trajectories.npz"

        url = 'https://www.uio.no/studier/emner/matnat/astro/AST2000/h20/blogger/Flukten%20fra%20Hoth/data/planet_trajectories.npz'
        if (os.path.exists(trajectories_filename) == False):
            try:
                import requests
                print(f'Downloading {trajectories_filename}')
                r = requests.get(url, allow_redirects=True)

                open(trajectories_filename, 'wb').write(r.content)
            except:
                print('You need to install requests to download file: pip install requests')

        # Load the saved planet trajectories
        trajectories = np.load(trajectories_filename,allow_pickle=True)
        self.times = trajectories[0]
        self.planet_positions = trajectories[1]
        print(f'Planet trajectories {self.planet_positions.shape[1]} planets loaded')

        self.system = SolarSystem(seed)          # Our famous system

    def planet_positions_at_t(self, t):

        positions = np.zeros((self.system.number_of_planets,2))
        for n in range(self.system.number_of_planets):
            positions[n] = self.planet_position(n, t)

        return positions

    def spacecraft_position_at_t(self, t):
        idx = (np.abs(self.rocket_trajectory_times - t)).argmin()

        spacecraft_position = self.spacecraft_positions[idx,:]
        spacecraft_velocity = self.spacecraft_velocities[idx,:]

        dt = t - self.rocket_trajectory_times[idx]

        movement = spacecraft_velocity * dt

        return spacecraft_position + movement

    def plot_planets(self,t):

        positions = self.planet_positions_at_t(t)

        plt.scatter(positions[:,0],positions[:,1])
        plt.scatter(0,0,color='r')
        plt.show()

    def planet_position(self, n, t):
        '''
        Find position of planet n at time t
        '''

        # Find index of time closest to our t
        idx = (np.abs(self.times - t)).argmin()
        self.log.info(f'Time is {t}')
        self.log.info(f'Planet position is {self.planet_positions[:,n,idx]} at index {idx} at time {self.times[idx]}')
        self.log.info(f'Planet position is {self.planet_positions[:,n,idx+1]} at index {idx+1} at time {self.times[idx+1]}')

        # Get planet position at time t
        planet_position = self.planet_positions[:,n,idx]

        # Interpolate position to exact time
        dt = t - self.times[idx]
        # Calculate orbital movement
        orbital_velocity = self.planet_orbital_velocity(n,t)

        self.log.info(f'Difference in time between position got and actual time {dt}')
        self.log.info(f'Orbital velocity is  {orbital_velocity}')
        orbital_movement = orbital_velocity*dt

        return planet_position + orbital_movement

    def planet_orbital_velocity(self,n,t):
        '''
        Calculate orbital velocity of the planet n at time t in years
        '''
        # t = 0 is special case
        if t == 0:
            return self.system.initial_velocities[:,n]

        # Find index of time closest to our t
        idx = (np.abs(self.times - t)).argmin()

        # Calculate delta t in years
        deltaT = (self.times[idx]-self.times[idx-1])

        # Calculate delta s in AU
        deltaS = (self.planet_positions[:,n,idx] - self.planet_positions[:,n,idx-1])
        # Planet orbital velocity vector in AU/year
        planet_velocity = deltaS/deltaT
        self.log.debug(f'Calculated orbital velocity {planet_velocity} = {deltaS}/{deltaT} : {np.linalg.norm(planet_velocity*c.AU/c.yr)} m/s')

        return planet_velocity

    def animate_orbits(self,t0,t,steps):
        fig, ax = plt.subplots()
        xdata, ydata = [], []
        ln, = plt.plot([], [], 'b.')
        ln_sp, = plt.plot([], [], 'r+',label='Spacecraft')
        ln_hoth, = plt.plot([], [], 'g.',label='Hoth')
        ln_pjuf, = plt.plot([],[], 'c.',label='Pjuf')

        def init():
            ax.set_xlim(-2, 2)
            ax.set_ylim(-2, 2)
            ax.set_xlabel('AU')
            ax.set_ylabel('AU')
            ax.legend(loc='upper left')
            plt.plot(0,0,'yo')
            return ln,

        def update(frame):
            pos = self.planet_positions_at_t(frame)
            sp_pos = np.array([self.spacecraft_position_at_t(frame)])
            #print(f'sp_pos: {sp_pos}')
            #all=np.append(pos,sp_pos,axis=0)

            xdata= pos[[1,2,3,4,5,7],0]
            ydata= pos[[1,2,3,4,5,7],1]
            ln.set_data(xdata, ydata)
            ln_sp.set_data(sp_pos[:,0],sp_pos[:,1])
            ln_hoth.set_data(pos[0,0],pos[0,1])
            ln_pjuf.set_data(pos[6,0],pos[6,1])
            return ln, ln_sp, ln_hoth, ln_pjuf


        anim = FuncAnimation(fig, update, frames=np.linspace(t0, t0+t, steps),
                            init_func=init, blit=True)

        anim.save('orbit_animation.gif', writer='imagemagick', fps=60)
        plt.show()

    def planet_html_table(self,t):

        print(f'<table><tr><th>Planet</th><th>Masse i sol masser</th><th>Posisjon i AU</th><th>Hastighet i AU/år</th></tr>')
        for n in range(self.system.number_of_planets):
            print(f'<tr><td>{n}</td><td>{self.system.masses[n]}</td><td>({self.planet_position(n,t)[0]},{self.planet_position(n,t)[1]})</td><td>({self.planet_orbital_velocity(n,t)[0]},{self.planet_orbital_velocity(n,t)[1]})')
        print(f'</table>')

    def simulate_trajectory(self, initial_time, initial_position, initial_velocity, simulation_time, dt):
        final_time = initial_time + simulation_time
        final_position = initial_position
        final_velocity = initial_velocity
        N = self.system.number_of_planets


        def f(r, t):
            '''
            Calculates force when we are at position r at time t
            And returns acceleration in AU/yr
            '''
            lrl = np.linalg.norm(r)     # distance in m
            #print(f'r = {r}, lrl = {lrl}, t = {t}')
            a_s = -c.G_sol*((self.system.star_mass)/lrl**3)*r # acceleration from gravitational force from star in m/s**2
            #print(f'as: {a_s} = (({self.system.star_mass*c.m_sun})/{lrl**3}*{r}')
            rp = r - self.planet_positions_at_t(initial_time+t) # Get planet positions at time t in meters
            ap = np.zeros((self.system.number_of_planets,2)) # Hold sum of gravitational force vectors for all planets
            for i in range(N):
                ap[i]  = -c.G_sol*((self.system.masses[i])/np.linalg.norm(rp[i])**3)*rp[i]
                #print(f'ap[{i}]: {ap[i]} = (({self.system.masses[i]*c.m_sun})/{np.linalg.norm(rp[i])**3})*{rp[i]}')

            a = a_s # Total accelleration from gravitational force
            for x in ap:
                a += x
            #print(f'a = {a}')
            #print(f'sum = {np.sum(a,axis=0)} at {t}')
            return a

        r0 = initial_position
        v0 = initial_velocity
        T = simulation_time

        self.spacecraft_positions, self.spacecraft_velocities, a, self.rocket_trajectory_times = tools.leapfrog(r0, v0, T, dt, f)

        return self.rocket_trajectory_times, self.spacecraft_positions, self.spacecraft_velocities, a

if __name__ == '__main__':

    spacecraft_mass = 1579
    tcomp = trajection_computer(spacecraft_mass)



    print(f'Planet 0 position at t=0 [{tcomp.planet_position(0,0)[0]},{tcomp.planet_position(0,0)[1]}]')
    initial_position = np.array([0.3658614707234236,8.831289590080655e-05])
    initial_velocity = np.array([2.5593576572450947,8.92172415360237])#np.array([2.5593576572450947,0.09843921360236853])
    initial_time = 0
    steps = 1000
    days = 21
    time = 0.6 #days*c.day/c.yr
    final_time, position, velocity, acc = tcomp.simulate_trajectory(initial_time,initial_position,initial_velocity,time,time/steps)

    print(f'RESULTS with {steps} steps during {time} years')
    print(f'    Initial time {initial_time} year')
    print(f'      Final time {final_time[-1]} year')
    print(f'Initial position {position[0]} AU')
    print(f'  Final position {position[-1]} AU')
    print(f'Initial velocity [{velocity[0,0]} {velocity[0,1]}] AU/yr')
    print(f'  Final velocity [{velocity[-1,0]} {velocity[-1,1]}] AU/yr')

    tcomp.animate_orbits(0,0.6,360)
    #plt.plot(position[:,0], position[:,1])
    #plt.scatter(position[0,0],position[0,1],marker='x')
    #plt.scatter(position[-1,0],position[-1,1],marker='o')

    #plt.show()
