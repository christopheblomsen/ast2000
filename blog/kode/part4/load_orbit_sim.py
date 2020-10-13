# Egen kode
import os
import sys
try:
    import cPickle as pickle
except:
    import pickle
import orbit_sim

def orbit_sim_factory(filename,args):
    '''
    Constructs or loads the orbit_sim object from file or URL
    '''
    url = 'https://www.uio.no/studier/emner/matnat/astro/AST2000/h20/blogger/Flukten%20fra%20Hoth/data/simulated_orbits.pkl'
    if (os.path.exists(filename) == False or args.download==True):
        try:
            import requests
            r = requests.get(url, allow_redirects=True)

            open(filename, 'wb').write(r.content)
        except:
            print('You need to install requests to download file: pip install requests')

    if (os.path.exists(filename) == False or args.run_sim==True):
        orbit = orbit_sim.orbit_sim()
        orbit.sim()

        with open(filename, "wb") as output:
            pickle.dump(orbit, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as input:
            orbit = pickle.load(input)

    return orbit
