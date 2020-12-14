from ast2000tools.relativity import RelativityExperiments
seed = 35511
experiments = RelativityExperiments(seed)

planet_idx = 6
experiments.spaceship_race(planet_idx,
                           filename_1='spaceship_race_frame_1.xml',
                           filename_2='spaceship_race_frame_2.xml',
