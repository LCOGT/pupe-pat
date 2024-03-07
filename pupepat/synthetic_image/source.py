import abc
import numpy as np


class Source(abc.ABC):
    def __init__(self, x=None, y=None):
        """
        A Source need not know it's position on a plate until it's been added.
        :param x:
        :param y:
        """
        self.x, self.y = x, y
        self.set_position(x, y)

        # the half_width value is sub-class dependent, but it's reason
        # for being is to provide a consistent value for the size of the
        # decal. like this: decal.shape[0] = 2*half_width + 1
        self.half_width = 0
        super().__init__()

    def set_position(self, x, y):
        self.x = x
        self.y = y

    @property
    def position(self):
        return self.x, self.y

    @property
    @abc.abstractmethod
    def decal(self):
        pass


class FlatSource(Source):
    def __init__(self, radius, counts, x=None, y=None):
        self.radius = radius
        self.counts = counts
        super().__init__(x=x, y=y)

        self.half_width = radius  # same as radius (for FlatSource)

    def __str__(self):
        return 'FlatSource(radius={radius:0.2f}, counts={counts:0g}, ' \
               'x={x:}, y={y:})'.format(radius=self.radius, counts=self.counts,
                                        x=self.x, y=self.y)

    @property
    def decal(self):
        shape = (2*self.half_width + 1, 2*self.half_width + 1)
        _decal = np.zeros(shape)

        ci, cj = self.half_width, self.half_width  # center of decal
        ix, jx = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]))

        # distance points in the index arrays (ix, jx) to the decal center (ci, cj)
        dist = np.sqrt((ix-ci)**2 + (jx-cj)**2)
        _decal[np.where(dist < self.radius)] = self.counts  # fill the circle with counts
        return _decal


class CircularPupil(Source):
    def __init__(self, inner_radius, outer_radius, counts,
                 center_vector=(0, 0), x=None, y=None):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.counts = counts
        self.center_vector = center_vector
        super().__init__(x=x, y=y)

        # make sure annulus makes sense (validate parameters)
        # inner_radius must not intersect the outer radius.
        if inner_radius + (center_vector[0]**2 + center_vector[1]**2) > outer_radius:
            raise ValueError('combined inner_radius and length of center_vector must not exceed outer_radius')

        self.half_width = outer_radius

    def __str__(self):
        return 'CircularPupil(inner_radius={inner_radius:0.2f}, outer_radius={outer_radius: 0.2f}, ' \
               'counts={counts:0g}, x={x:}, y={y:})'.format(inner_radius=self.inner_radius,
                                                            outer_radius=self.outer_radius,
                                                            counts=self.counts, x=self.x, y=self.y)

    @property
    def decal(self):
        shape = (2*self.half_width + 1, 2*self.half_width + 1)
        _decal = np.zeros(shape, dtype=float)

        ci, cj = self.half_width, self.half_width  # center of decal
        ix, jx = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]))

        # calculate dist of all points (ix, jx) to the decal center (ci, cj)
        dist = np.sqrt((ix-ci)**2 + (jx-cj)**2)
        _decal[np.where(dist < self.outer_radius)] = self.counts  # fill the circle with counts

        # remake center and dist array for offset inner_radius
        ci -= self.center_vector[0]
        cj -= self.center_vector[1]
        dist = np.sqrt((ix-ci)**2 + (jx-cj)**2)
        _decal[np.where(dist < self.inner_radius)] = 0.0  # clear the counts within the inner_radius
        return _decal


class EllipticalPupil(Source):
    def __init__(self,
                 inner_a, inner_b, inner_phi,
                 outer_a, outer_b, outer_phi,
                 counts,
                 center_vector=(0, 0),
                 x=None, y=None):
        super().__init__(x=x, y=y)
        self.inner_a = inner_a
        self.inner_b = inner_b
        self.inner_phi = inner_phi
        self.outer_a = outer_a
        self.outer_b = outer_b
        self.outer_phi = outer_phi
        self.counts = counts
        # center_vector is the vector from inner to outer center
        self.center_vector = center_vector

        # TODO: possible validations:
        # A > B, so we can make that assumption (swap them, don't throw)
        # outer_b > inner_a
        # outer_b > inner_a + max(center_vector) (inner doesn't intersect outer)
        # outer_b > inner_a at worst when inner_phi and outer_phi orthogonal
        # etc

        self.half_width = max(outer_a, outer_b)  # major semi-axis

    def __str__(self):
        return 'EllipticalPupil(inner_a={inner_a:0.2f}, inner_b={inner_b:0.2f}, inner_phi={inner_phi:0.2f}, ' \
               'outer_a={outer_a:0.2f}, outer_b={outer_b:0.2f}, outer_phi={outer_phi:0.2f}, counts={counts:0g} ' \
               'x={x:}, y={y:})'.format(inner_a=self.inner_a, inner_b=self.inner_b, inner_phi=self.inner_phi,
                                        outer_a=self.outer_a, outer_b=self.outer_b, outer_phi=self.outer_phi,
                                        counts=self.counts, x=self.x, y=self.y)

    @property
    def decal(self):
        shape = (2*self.half_width + 1, 2*self.half_width + 1)
        _decal = np.zeros(shape, dtype=float)

        center_i, center_j = self.half_width, self.half_width  # center of decal
        X, Y = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]))

        # calculate dist of all points (ix, jx) to the decal center (ci, cj)
        dist = self.elliptical_distance_meshgrid(self.outer_a, self.outer_b, self.outer_phi,
                                                 center_i, center_j, X, Y)
        _decal[np.where(dist < 1.)] = self.counts  # fill the circle with counts

        # remake center and dist array for offset inner_a, -b
        center_i -= self.center_vector[0]
        center_j -= self.center_vector[1]
        dist = self.elliptical_distance_meshgrid(self.inner_a, self.inner_b, self.inner_phi,
                                                 center_i, center_j, X, Y)
        _decal[np.where(dist < 1.)] = 0  # clear the counts with inner_a, _b
        return _decal

    def elliptical_distance_meshgrid(self, a, b, phi, center_x, center_y, X, Y):
        return (((np.cos(phi) * (X - center_x) + np.sin(phi) * (Y - center_y)) ** 2 / a ** 2) +
                ((-np.sin(phi) * (X - center_x) + np.cos(phi) * (Y - center_y)) ** 2 / b ** 2))
