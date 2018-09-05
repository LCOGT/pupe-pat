"""
pupepat.utils - Utility functions for pupil plate fits

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0

February 2018
"""
from astropy.stats import median_absolute_deviation
import numpy as np
from astropy.table import Table
from astropy.io import fits, ascii
import collections
import os
import sep
import logging
from PyPDF2 import PdfFileReader, PdfFileWriter

logger = logging.getLogger('pupepat')

#
# configuraton for everything: image filter, cutout, solver guesses, etc.
# this dictionary holds the default values. if --config-file is specified,
# that file is read and this dict is updated. It's written to the --output-dir
# as yaml file
#
# guess for the size of a defocused donut is 0.387" 4 FOCDMD
#
config = {
    'SEP': {
        'extract_SN_threshold': 5.0,
        'mask_threshold_scale_factor': 10.0,
        'min_area': 2000, # see RSilverd comment in story with eqn for this
        'box_size': 128,
    },
    'cutout': {
        'cutout_radius': 150,
        'cutout_edge_limit': 150,
    },
    'image_filters': {'bad_column_max': 400,   # not_plumbed
                      'bad_column_min': 200,  # not_plumbed
                      'edge_dist_px': 150,  # not_plumbed
                      'in_focus_scale_factor': 200,  # not_plumbed
    },
}

def update_mapping(dest, source):
    '''Update dest dictionary with key,value pairs from source dictionary.
    '''
    for key, value in source.items():
        if isinstance(value, collections.Mapping):
            dest[key] = update_mapping(dest.get(key, {}), value)
        else:
            dest[key] = value
    return dest


def offsets(x1, y1, x2, y2):
    return ((x2 - x1) ** 2.0 + (y2 - y1) ** 2.0) ** 0.5


def get_bias_corrected_data_in_electrons(hdu):
    """
    Estimate the bias level of image and substract it from the image.
    If bias estimate is negative, issue warning and scale data.

    Uses median_absolute_deviation (MAD) to compute standard deviation.
    Note: does not substract background from data to compute noise, as was done
    prior to Sept.2018 parameter tuning.

    :param hdu: fits hdu with the image in question

    :return: np.array data_e
    """
    gain = float(hdu[0].header['GAIN'])            # e-/count
    read_noise_e = float(hdu[0].header['RDNOISE']) # e-/pixel
    data_e = gain * hdu[0].data                    # counts to electrons

    # 1.48 here goes from median absolute deviation to standard deviation
    noise_e = 1.48 * median_absolute_deviation(data_e)
    #noise_e = 1.48 * median_absolute_deviation(hdu[0].data - sep.Background(hdu[0].data.astype(np.float)))

    estimated_bias_level_in_electrons = np.median(data_e) - noise_e * noise_e + read_noise_e * read_noise_e

    if estimated_bias_level_in_electrons < 0:
        noise_e = 1.48 * median_absolute_deviation(data_e)
        sqrt_median_e = np.sqrt(np.median(data_e))
        scale_factor = sqrt_median_e / noise_e
        warn = '''Negative bias {:0.2f}. Scaling data by (sqrt(median)/noise): ({:0.2f}/{:0.2f})={:0.2f}'''
        logger.warning(warn.format(estimated_bias_level_in_electrons, scale_factor, sqrt_median_e, noise_e ))
        data_e *= scale_factor
    else:
        # bias corrected data in electrons
        data_e = data_e - estimated_bias_level_in_electrons

    return data_e


def save_results(input_filename: str, best_fit_models: list,
                 output_table: Table, output_filename: str):
    """
    Save the best fit model to a text file and into the table

    :param input_filename: Filename of the image that was fit
    :param best_fit_models: list of astropy.modeling.Model for the best fit elliptical annulus
    :param output_table: astropy.tables.Table with the previous results
    :param output_filename: filename to save the results
    :return: output_table
    """
    header_keywords = ['M2ROLL', 'MOLTYPE', 'EXPTIME', 'FOCDMD', 'M2PITCH']
    hdu = fits.open(input_filename)

    if output_table is None:
        param_names = [param for param in best_fit_models[0].keys()]
        output_table = {param: [best_fit_model[param] for best_fit_model in best_fit_models]
                        for param in param_names}
        output_table['filename'] = [os.path.basename(input_filename)] * len(best_fit_models)
        for keyword in header_keywords:
            output_table[keyword] = [hdu[0].header[keyword]] * len(best_fit_models)
        output_table = Table(output_table)
    else:
        for i, best_fit_model in enumerate(best_fit_models):
            best_fit_model['filename'] = os.path.basename(input_filename)
            for keyword in header_keywords:
                best_fit_model[keyword] = hdu[0].header[keyword]
            logger.debug(best_fit_model)
            output_table.add_row(best_fit_model)

    output_table.write(output_filename, format='ascii', overwrite=True)
    return output_table


def make_cutout(data, x0, y0, r):
    '''Slice the data array centered on (x0,y0) from -r to r+1 in both dimensions.
    '''
    return data[int(y0 - r):int(y0 + r + 1), int(x0 - r):int(x0 + r + 1)]


def cutout_coordinates(cutout, x0, y0):
    x, y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))
    return np.sqrt((x - x0) ** 2.0 + (y - y0) ** 2.0)


def run_sep(data, mask_threshold=None):
    """Compute mask and background; then sep.extract sources.
    """
    extract_SN_threshold = config['SEP']['extract_SN_threshold']
    mask_threshold_scale_factor = config['SEP']['mask_threshold_scale_factor']
    box_size = config['SEP']['box_size']
    min_area = config['SEP']['min_area']

    if mask_threshold is None:
        mask_threshold = mask_threshold_scale_factor * np.sqrt(np.median(data)) + np.median(data)

    background = sep.Background(np.ascontiguousarray(data),
                                mask=np.ascontiguousarray(data > mask_threshold),
                                bw=box_size, bh=box_size)
    return sep.extract(data - background, extract_SN_threshold,
                       err=np.sqrt(data), minarea=min_area, deblend_cont=1.0, filter_kernel=None)


def prettify_focdmd(demanded_focus):
    return '{focdmd:+.0f}'.format(focdmd=demanded_focus).replace('-', 'm').replace('+', 'p')


def merge_pdfs(output_directory, output_table, output_pdf='pupe-pat'):
    table_path = os.path.join(output_directory, output_table)

    # Short circuit
    if not os.path.exists(table_path):
        return
    data = ascii.read(table_path)
    data.sort(['FOCDMD', 'filename'])
    pdf_files = np.array([os.path.join(output_directory,
                                       '{basename}_{sourceid}.pdf'.format(basename=os.path.splitext(row['filename'])[0],
                                                                          sourceid=row['sourceid']))
                          for row in data])
    for demanded_focus in np.unique(data['FOCDMD']):
        focus_set_indexes = data['FOCDMD'] == demanded_focus

        pdf_writer = PdfFileWriter()
        for pdf_file in pdf_files[focus_set_indexes]:
            pdf_writer.appendPagesFromReader(PdfFileReader(pdf_file))

        output_filename = os.path.join(output_directory,
                                       '{basename}-{focdmd}-donut.pdf'.format(basename=output_pdf,
                                                                              focdmd=prettify_focdmd(demanded_focus)))
        with open(output_filename, 'wb') as output_stream:
            pdf_writer.write(output_stream)


def estimate_scatter(a, axis=None):
    """
    Caclulate the inter-quartile range of *a* (scaled to normal).
    :param a: samples for the distribution
    :param axis: axis to calculate the scatter
    :return: scatter
    """
    # 2 Φ^−1(0.75)σ ≈ 1.349σ ≈ (27/20)σ for the interquartile range
    # 0.741301109 = 1.0 / (2**0.5 * special.erfinv(1.5 - 1.0) * 2.0)
    lower_quartile, upper_quartile = np.percentile(a, [25, 75], axis=axis)
    return 0.741301109 * (upper_quartile - lower_quartile)


def get_inliers(data, threshold):
    """
    Get the indexes of data that are within the threshold of the median
    :param data: Data to check for inliers
    :param threshold: n sigma threshold to consider a sample an inlier
    :return: Boolean array with True for inliers
    """
    scatter = estimate_scatter(data)
    scaled_absolute_deviation = np.abs(data - np.median(data)) / scatter
    return scaled_absolute_deviation <= threshold


def get_primary_extension(hdu: fits.HDUList):
    """
    Figure out which is the primary extension (1 if fpacked, 0 otherwise)
    :param hdu: astropy.io.fits.HDUList
    :return: int
    """
    _, file_extension = os.path.splitext(hdu.filename())
    if file_extension == '.fz':
        primary_extension = 1
    else:
        primary_extension = 0

    return primary_extension


def is_unstitched(hdu):
    """
    Check if a file has been stitched or not

    If the file has NAXIS3 in the header or there are multiple SCI extensions, we assume the file has not been
    stitched yet.
    :param hdu: astropy fits HDU list
    :rtype: bool
    """
    primary_extension = get_primary_extension(hdu)

    if 'NAXIS3' in hdu[primary_extension].header:
        return True

    sci_extension_counter = 0
    # Apparently passing false to this function makes it return tuples instead of printing to screen
    hdu_info = hdu.info(False)
    for hdu_row in hdu_info:
        if 'SCI' in hdu_row:
            sci_extension_counter += 1
    return sci_extension_counter > 1


def should_process_image(path, proposal_id='LCOEngineering'):
    """
    Test if we should try to process a given image.

    We assume that all images of interest will be defocused (abs(FOCDMD) > 0) and either an experimental or
    exposure obstype. We only consider images taken with the given proposal id. If the frame is from a Sinistro, we
    only consider frames that are already stitched.

    :param path: Full path to the image of interest
    :param proposal_id: Proposal ID for images to process
    :rtype: bool
    """
    hdu = fits.open(path)

    if is_unstitched(hdu):
        logger.error('File appears to be unstitched.', extra={'tags': {'filename': path}})
        return False

    primary_extension = get_primary_extension(hdu)

    header = hdu[primary_extension].header
    should_process = np.abs(header.get('FOCDMD', 0.0)) > 0.0
    should_process &= header.get('PROPID') == proposal_id
    should_process &= header.get('OBSTYPE') == 'EXPOSE' or header.get('OBSTYPE') == 'EXPERIMENTAL'
    return should_process
