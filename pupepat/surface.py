"""
pupepat.surface - Tools to fit a 2-D polynomial to a 2D surface

Author
    Curtis McCully (cmccully@lco.global)
    Rob Siverd

License
    GPL v3.0

February 2018
"""

import numpy as np


class SurfaceFitter(object):
    def __init__(self, degree):
        """
        param degree: maximum of the sum of x and y powers for the polynomial
        """
        self.exponents = list(self._get_exponent(degree))
        self.coefficients = None

    @staticmethod
    def _get_exponent(degree):
        """
        Make a list of exponents for the polynomial terms
        :param degree: Maximum of the sum of the x and y powers
        :return: num_exponents, list of x,y powers to use for the polynomial
        """
        x_powers, y_powers = np.meshgrid(range(degree + 1), range(degree + 1))
        powers_to_use = (x_powers + y_powers <= degree)
        return zip(x_powers[powers_to_use], y_powers[powers_to_use])

    def fit(self, x, y, z):
        """
        Fit a surface
        :param x: x coordinate grid for the surface
        :param y: y coordinate grid for the surface
        :param z: 2D surface data
        """
        design_matrix = np.zeros((x.size, len(self.exponents)), dtype=float)
        for k, (i, j) in enumerate(self.exponents):
            design_matrix[:, k] = x.flatten() ** i * y.flatten() ** j
        self.coefficients, _, _, _ = np.linalg.lstsq(design_matrix, z.flatten(), rcond=None)

    def eval(self, x, y):
        """
        Evaluate a surface given a grid of x and y coordinates
        :param x: x coordinate grid for the surface
        :param y: y coordinate grid for the surface
        :return: Surface evalueted at x, y
        """
        z = np.zeros(x.shape)
        for a, (i, j) in zip(self.coefficients, self.exponents):
            z += a * x ** i * y ** j
        return z
