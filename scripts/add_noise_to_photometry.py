#! /usr/bin/python

import sys
import numpy as np
from numpy.random import normal

depths = np.array(sys.argv[1:]).astype(np.float32)

noise_magnitude = 0.2 * np.power(10, -0.4 * (depths - 23.9))

for line in sys.stdin:
    line_contents = line.split()
    object_id = line_contents[0]
    ra = line_contents[1]
    dec = line_contents[2]
    fluxes = np.array(line_contents[9:]).astype(np.float32)
    noise = noise_magnitude * normal(size=noise_magnitude.size)
    output_flux = fluxes + noise
    error = noise_magnitude
    # interleave the fluxes and flux errors
    final_result = np.empty((output_flux.size + error.size,), dtype=error.dtype)
    final_result[0::2] = output_flux
    final_result[1::2] = error
    print(" ".join([object_id, ra, dec]) + " " + " ".join(final_result.astype(str)))
