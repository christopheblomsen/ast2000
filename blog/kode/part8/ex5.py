from ast2000tools.relativity import RelativityExperiments
seed = 35511
experiments = RelativityExperiments(seed)

planet_idx = 6
experiments.laser_chase(planet_idx,
                        filename_1='laser_chase_frame_1.xml',
                        filename_2='laser_chase_frame_2.xml')
