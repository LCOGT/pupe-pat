"""
pupepat.ellipse - Tools for generating ellipses

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0

February 2018
"""
import numpy as np


def inside_ellipse(x, y, x0, y0, a, b, theta):
    """
    Calculate if x and y are in an ellipse.

    :param x: X coordinates to check if in ellipse
    :param y: Y coordinates to check if in ellipse
    :param x0: X center of the ellipse
    :param y0: Y center of the ellipse
    :param a: semimajor axis
    :param b: semiminor axis
    :param theta: orientation of the semimajor axis
    :return: boolean numpy array
    """
    x_recentered = x - x0
    y_recentered = y - y0

    semimajor_term = x_recentered * np.cos(theta) + y_recentered * np.sin(theta)
    semimajor_term /= a
    semimajor_term **= 2.0

    semiminor_term = y_recentered * np.cos(theta) - x_recentered * np.sin(theta)
    semiminor_term /= b
    semiminor_term **= 2.0

    return semimajor_term + semiminor_term <= 1.0


def generate_ellipse(x0, y0, a, b, theta):
    """
    Generate the x and y values of an ellipse to plot.

    :param x0: X coordinate of the center of the ellipse
    :param y0: Y coordinate of the center of the ellipse
    :param a: semi-major axis of the ellipse
    :param b: semi-minor axis of the ellipse
    :param theta: orientation of the semimajor axis (radians) measured counterclockwise from the X-axis.

    :return: x, y arrays of coordinates in the ellipse
    """
    # Angle around the ellipse
    phi = np.arange(0, 2.0 * np.pi, 0.05)
    x_unrotated = a * np.cos(phi)
    y_unrotated = b * np.sin(phi)
    # Rotate the coordinates

    x = x_unrotated * np.cos(theta) - y_unrotated * np.sin(theta)
    y = x_unrotated * np.sin(theta) + y_unrotated * np.cos(theta)

    x += x0
    y += y0

    # Close the ellipse back on itself
    np.append(x, x[0])
    np.append(y, y[0])
    return x, y


def calculate_donut_offset(inner_coordinate, outer_coordinate, sbig=False):
    """
    Calculate the offset of the inner coordinate from the outer coordinate.
    Will return the offset in pixels unless sbig is set to True, in which case it will return the offset in mm of M1
    offset based on calculations by Patrick Conway using zemax in Feb. 2018.
    (https://lcogt.slack.com/archives/C8V95PF9T/p1519653167000491)

    :param inner_coordinate: inner coordinate (x or y) in pixels
    :param outer_coordinate: outer coordinate (x or y) in pixels
    :param sbig: set to True if using a kb camera, and you want the offset to be in mm of M1 offset
    :return: offset
    """

    pixel_to_m1_mm = 9.679  # mm of M1 offset per pixel for an sbig on a 1m telescope

    offset = inner_coordinate - outer_coordinate  # pixels
    if sbig:
        offset /= pixel_to_m1_mm  # mm of M1 offset
    return offset
