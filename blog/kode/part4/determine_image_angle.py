# Egen kode
from PIL import Image
import numpy as np
import argparse
import astrogation_computer

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='The filename to find the angle for')
args = parser.parse_args()

astro_comp = astrogation_computer.astrogation_computer()

angle = astro_comp.find_orientation_angle(args.filename)

print(f'Best angle match is {angle}')
