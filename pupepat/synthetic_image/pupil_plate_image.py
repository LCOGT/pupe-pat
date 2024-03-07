import numpy as np
from astropy.io import fits
from pupepat.synthetic_image.source import Source



class PupilPlateImage:
    def __init__(self, shape=(4096, 4096), seed=42):
        self.shape = shape
        self._random = np.random.RandomState(seed)  # mine

        self.sources = []
        self.cutout_margin_factor = 1.2

        # Default Noise parameters
        self.sky_level = 200.
        self.read_noise = 10.
        self.add_poisson_noise = True

    def render(self) -> np.array:
        """
        Given the properties that have been set (background, sources, noise, etc)
        generate the pupil plate image as an np.array

        :return: image
        """
        image = np.zeros(self.shape, dtype=float)

        # add the background sky level
        image += self.sky_level

        # add the sources
        for source in self.sources:
            # add (+=) source decal into the appropriate slice of the data
            slice_index_x = source.x - source.half_width
            slice_index_y = source.y - source.half_width

            decal = source.decal
            image[slice_index_x:slice_index_x + decal.shape[0],
                  slice_index_y:slice_index_y + decal.shape[1]] += decal

        # add Poisson noise (shot noise)
        if self.add_poisson_noise:
            image = np.float_(self._random.poisson(image))

        # add read noise
        read_noise = self._random.normal(scale=self.read_noise, size=image.size)
        read_noise.shape = image.shape
        image += read_noise

        return image

    def __str__(self) -> str:
        return super().__str__()

    def read_fits_from(self, filename):
        self.filename = filename

        # assume for now that we only read FITS file that we've written
        # and that those have known properties
        hdu_list = fits.open(filename)
        assert(len(hdu_list) == 1)
        self.hdu_primary = hdu_list[0]
        hdu_list.close()

    def write_fits_to(self, filename):
        """
        Basically following
        http://docs.astropy.org/en/stable/io/fits/#creating-a-new-fits-file
        http://docs.astropy.org/en/stable/io/fits/#using-astropy-io-fits
        :param filename:
        :return:
        """
        fits.writeto(filename, self.render(), fits.Header())

    def add_source(self, source, x, y):
        assert isinstance(source, Source)
        # set source position and add it to sources
        source.set_position(x, y)
        self.sources.append(source)

    def cutout_for_source(self, source):
        """The cutout is sliced from the PupilPlateImage.

        Any background or PSF noise is included in the cutout.

        :param source: the source to be extracted
        :return: np.array of cutout
        """
        cutout_half_width = int(self.cutout_margin_factor * source.half_width)
        cx_start = source.x - cutout_half_width
        cx_end = source.x + cutout_half_width + 1
        cy_start = source.y - cutout_half_width
        cy_end = source.y + cutout_half_width + 1

        image = self.render()
        cutout = image[cx_start:cx_end, cy_start:cy_end]
        return cutout

