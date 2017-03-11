import numpy as np
import math

max_redshift = 4.45
dz = 0.00175
n_templates = 422 #188
red_max = 1.0
dred = 0.045
n_reddening_laws = 5

widths = [ 0, 14.9, 57.9, 100.6, 6.1, 51.8, 31.4, 176.9, 99.1, 0, 11.6, 0, 3.6, 16.8,
           16.4, 106, 0, 0, 29.6, 321.9, 360.5, 408.8, 20, 107.3, 46.3, 0, 2.6, 0,
           3.4, 11.1, 63.4, 0, 20.2, 26.5, 32.6, 22.5, 56, 17.1, 27.7, 34.5, 3.4,
           30.3, 42.8, 20.1, 0, 11.4, 19.3, 108.6, 7.8, 0, 8.1, 5.8, 0, 60.2, 39.6,
           14.6, 23.6, 0, 5.4, 0, 67.7, 27.5, 15.2, 0, 34.5, 0, 0.4, 0, 1.5, 0, 18.9,
           0, 0, 0, 31.8, 2.5, 1.7, 0, 0, 18.6, 33, 0, 90.9, 6.7, 1.4, 4.1, 11.1, 0,
           0, 10.1, 7.5, 9, 14.4, 1.6, 59, 55.9, 19.9, 38.9, 29, 0, 27.7, 47.5, 122.2,
           109.1, 28.8, 3.5, 0, 15.2, 54.2, 98.3, 32.8, 58.4, 112.4, 12.3, 13.7, 19.3,
           154.4, 265.5, 19.2, 19.2, 35.7, 19.9, 21.1, 29.9, 10.4, 0, 15.6, 136.7,
           218.1, 361.4 ]

id = 1

for redshift in np.arange(0, max_redshift, dz):
    for template in range(1, n_templates + 1):
        for reddening in np.arange(0, red_max, dred):
            for rlaw in range(n_reddening_laws):

                width = 0

                if template < 130:
                    width = widths[template]

                line_scale = 1.0

                if width > 0 and redshift > 2.75:
                    line_scale = 0.65 * redshift - 0.78  # approximate fit to the measured correction factor

                print id, 150.0, 2.0, 20, redshift, template, reddening, rlaw + 1, 1.0, line_scale

                id += 1
