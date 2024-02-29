"""
pupepat.fitting - Tools to fit models to pupil plate images

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0

February 2018
"""
import numpy as np
from astropy.io import fits
import astroscrappy
from pupepat.ellipse import inside_ellipse
from pupepat.utils import make_cutout, cutout_coordinates, run_sep, config
from pupepat.utils import get_bias_corrected_data_in_electrons, source_is_valid
from pupepat.plot import plot_best_fit_ellipse
from types import SimpleNamespace
import os
import emcee
import corner
from scipy.optimize import minimize

import logging

logger = logging.getLogger('pupepat')


def elliptical_annulus(x, y, x0_inner=0.0, y0_inner=0.0, a_inner=1.0, b_inner=1.0, theta_inner=0.0, amplitude_inner=1.0,
                       x0_outer=0.0, y0_outer=0.0, a_outer=1.0, b_outer=1.0, theta_outer=0.0, amplitude_outer=1.0,
                       x_slope=0.0, y_slope=0.0, background=0.0):
    """
    2D Elliptical Annulus

    :param x: X positions to evaluate the model
    :param y: Y positions to evaluate the model
    :param x0_inner: X center of the inner boundary of the annulus
    :param y0_inner: Y center of the inner boundary of the annulus
    :param a_inner: Semi-major axis length of the inner boundary of the annulus
    :param b_inner: Semi-minor axis length of the inner boundary of the annulus
    :param theta_inner: The rotation angle of the semi-major axis in radians of the inner boundary of the annulus
    :param amplitude_inner: Amplitude of the region inside the inner boundary of the annulus
    :param x0_outer: X center of the outer boundary of the annulus
    :param y0_outer: Y center of the outer boundary of the annulus
    :param a_outer:  Semi-major axis length of the outer boundary of the annulus
    :param b_outer:  Semi-minor axis length of the outer boundary of the annulus
    :param theta_outer: The rotation angle of the semi-major axis in radians of the outer boundary of the annulus
    :param amplitude_outer: Amplitude of the region inside the annulus
    :param x_slope: Slope in the X-direction for the flux in the annulus
    :param y_slope: Slope in the Y-direction for the flux in the annulus
    :param background: Background level

    :return: numpy array with amplitude_inner inside the inner boundary of the annulus and amplitude_outer in the
             annulus
    """

    result = np.zeros(x.shape)
    # Include a gradient in the flux
    inside_outer = inside_ellipse(x, y, x0_outer, y0_outer, a_outer, b_outer, theta_outer)
    result[inside_outer] = amplitude_outer + x_slope * x[inside_outer] + y_slope * y[inside_outer]
    result[inside_ellipse(x, y, x0_inner, y0_inner, a_inner, b_inner, theta_inner)] = amplitude_inner
    result[~inside_outer] = background
    return result


def fit_defocused_image(filename, plot_basename):
    logger.info('Extracting sources', extra={'tags': {'filename': os.path.basename(filename)}})
    hdu = fits.open(filename)
    if not hdu[0].header.get('L1STATBI', None):
        data = get_bias_corrected_data_in_electrons(hdu)
    else:
        data = hdu[0].data

    # repair any hot pixels and/or cosmic rays (gain=1.0 b/c we're already in units of electrons)
    hpcr_mask, cleaned_data = astroscrappy.detect_cosmics(data, readnoise=hdu[0].header['RDNOISE'], gain=1.0,
                                                          sigclip=20.0, sigfrac=0.5, objlim=20.0, cleantype='idw')

    data = cleaned_data
    sources = run_sep(data, hdu[0].header)
    logger.info('Found {num_s} sources'.format(num_s=len(sources)),
                extra={'tags': {'filename': os.path.basename(filename)}})

    best_fit_models = [fit_cutout(data, source, plot_basename, os.path.basename(filename),
                                  hdu[0].header, i)
                       for i, source in enumerate(sources)]
    return best_fit_models


def cutout_has_multiple_sources(data, cutout, header, background):
    # Run sep again and make sure there is only one source in the cutout
    cutout_sources = run_sep(cutout, header, background_mask_threshold=25.0 * np.sqrt(background) + background)

    if len(cutout_sources) > 1:
        logger.error('Too many sources detected in cutout. Likely source crowding.')

    return len(cutout_sources) > 1


def fit_cutout(data, source, plot_filename, image_filename, header, id, fit_circle=True, fit_gradient=False):
    logger.debug('Validating source', extra={'tags': {'filename': image_filename, 'x': source['x'], 'y': source['y']}})
    if not source_is_valid(data, source):
        return
    logger.info('Fitting source', extra={'tags': {'filename': image_filename, 'x': source['x'], 'y': source['y']}})

    cutout_radius = config['cutout']['cutout_radius']
    # backed out this change because it needs more testing
    #cutout_padding = config['cutout']['cutout_padding']
    #pad = cutout_padding * max((source['ymax']-source['ymin']), (source['xmax']-source['xmin']))
    #cutout = data[int(source['ymin']-pad):int(source['ymax']+pad),
    #              int(source['xmin']-pad):int(source['xmax']+pad)]
    cutout = make_cutout(data, source['x'], source['y'], cutout_radius)
    r = cutout_coordinates(cutout, cutout_radius + 1, cutout_radius + 1)

    # backed out this change because it needs more testing
    ##inner_radius_guess = 3.4 * np.abs(header['FOCDMD']) # Rob Siverd, personal communiation
    ##outer_radius_guess = 3.0 * inner_radius_guess \
    ##                     if header['FOCDMD'] > 0 else 2.25 * inner_radius_guess

    inner_brightness_guess = np.median(cutout[r < 5])
    inside_donut_guess = cutout > (20.0 * np.sqrt(np.median(inner_brightness_guess)) + inner_brightness_guess)
    inner_radius_guess = np.percentile(r[inside_donut_guess].flatten(), 5)

    # Remake the cutout with +- 4 * inner_radius_guess
    cutout_radius = 4.0 * inner_radius_guess
    cutout = make_cutout(data, source['x'], source['y'], cutout_radius)
    x0 = y0 = cutout_radius + 1.0
    r = cutout_coordinates(cutout, x0, y0)

    background = np.median(data)
    if cutout_has_multiple_sources(data, cutout, header, background):
        return

    inside_donut_guess = cutout > (20.0 * np.sqrt(np.median(inner_brightness_guess)) + background)
    outer_radius_guess = np.percentile(r[inside_donut_guess].flatten(), 90)

    brightness_guess = np.median(cutout[np.logical_and(r < outer_radius_guess, r > inner_radius_guess)])
    inner_brightness_guess = np.median(cutout[r < inner_radius_guess])

    x, y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))
    error = np.sqrt(np.abs(cutout))
    if fit_circle:
        best_fit_parameter_names = ['x0_inner', 'y0_inner', 'r_inner', 'amplitude_inner',
                                    'x0_outer', 'y0_outer', 'r_outer', 'amplitude_outer', 'background']
        p0 = [x0, y0, inner_radius_guess, inner_brightness_guess, x0, y0, outer_radius_guess,
              brightness_guess, background]
        ln_likelihood_function = ln_likelihood_circle

    else:
        best_fit_parameter_names = ['x0_inner', 'y0_inner', 'a_inner', 'b_inner', 'theta_inner', 'amplitude_inner',
                                    'x0_outer', 'y0_outer', 'a_outer', 'b_outer', 'theta_outer', 'amplitude_outer',
                                    'x_slope', 'y_slope', 'background']
        p0 = [x0, y0, inner_radius_guess, inner_radius_guess, inner_brightness_guess, x0, y0,
              outer_radius_guess, outer_radius_guess, brightness_guess, 0.0, 0.0, background]
        ln_likelihood_function = ln_likelihood_ellipse

    best_fit_solution = minimize(lambda *args: -ln_likelihood_function(*args), p0, method='Powell',
                                 args=(x, y, cutout, error), options={'maxiter': 20000})
    logger.debug(best_fit_solution)

    best_fit_parameters = {best_fit_parameter_name: best_fit_parameter
                           for best_fit_parameter_name, best_fit_parameter in zip(best_fit_parameter_names,
                                                                                  best_fit_solution['x'])}
    best_fit_model = SimpleNamespace(**best_fit_parameters)
    best_fit_model.param_names = best_fit_parameter_names
    if fit_circle:
        best_fit_model.a_inner = best_fit_model.r_inner
        best_fit_model.b_inner = best_fit_model.r_inner
        best_fit_model.a_outer = best_fit_model.r_outer
        best_fit_model.b_outer = best_fit_model.r_outer
        best_fit_model.theta_inner = 0.0
        best_fit_model.theta_outer = 0.0
        best_fit_model.x_slope = 0.0
        best_fit_model.y_slope = 0.0

    plot_best_fit_ellipse(plot_filename + '_{id}.pdf'.format(id=id), cutout, best_fit_model, header)

    logging_tags = {parameter: getattr(best_fit_model, parameter) for parameter in best_fit_model.param_names}
    logging_tags['filename'] = image_filename
    for keyword in ['xmin', 'xmax', 'ymin', 'ymax']:
        logging_tags[keyword] = float(source[keyword])

    logger.debug('Best fit parameters for PUPE-PAT model',
                 extra={'tags': logging_tags})

    fit_results = {param: getattr(best_fit_model, param) for param in best_fit_model.param_names}
    fit_results['x0_centroid'], fit_results['y0_centroid'] = x0, y0
    fit_results['sourceid'] = id
    return fit_results


def ln_likelihood_ellipse(theta, x, y, data, data_error):
    x0_inner, y0_inner, a_inner, b_inner, theta_inner, amplitude_inner, \
    x0_outer, y0_outer, a_outer, b_outer, theta_outer, amplitude_outer, \
    x_slope, y_slope, background = theta
    model = elliptical_annulus(x, y, x0_inner=x0_inner, y0_inner=y0_inner, a_inner=a_inner, b_inner=b_inner,
                               theta_inner=theta_inner,
                               amplitude_inner=amplitude_inner, x0_outer=x0_outer, y0_outer=y0_outer, a_outer=a_outer,
                               b_outer=b_outer, theta_outer=theta_outer, amplitude_outer=amplitude_outer,
                               x_slope=x_slope,
                               y_slope=y_slope, background=background)
    inv_sigma2 = 1.0 / (data_error ** 2.0)
    return -0.5 * (np.sum((data - model) ** 2 * inv_sigma2 - np.log(inv_sigma2)))


def run_mcmc(x, y, data_error, theta0, data, likelihood_function, labels, output_plot='pupe-pat-mcmc.pdf', nwalkers=50):
    n_dim = len(theta0)
    starting_positions = [np.array(theta0) * (1.0 + np.random.uniform(-0.1, 0.1, len(theta0))) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, n_dim, likelihood_function, args=(x, y, data, data_error), threads=4)
    sampler.run_mcmc(starting_positions, 1000)
    samples = sampler.chain[:, 200:, :].reshape((-1, n_dim))
    fig = corner.corner(samples, labels=labels)
    fig.savefig(output_plot)


def ln_likelihood_circle(theta, x, y, data, data_error):
    x0_inner, y0_inner, r_inner, amplitude_inner, \
    x0_outer, y0_outer, r_outer, amplitude_outer, background = theta
    model = elliptical_annulus(x, y, x0_inner=x0_inner, y0_inner=y0_inner, a_inner=r_inner, b_inner=r_inner,
                               theta_inner=0.0, amplitude_inner=amplitude_inner, x0_outer=x0_outer,
                               y0_outer=y0_outer, a_outer=r_outer, b_outer=r_outer, theta_outer=0.0,
                               amplitude_outer=amplitude_outer, x_slope=0.0, y_slope=0.0, background=background)
    inv_sigma2 = 1.0 / (data_error ** 2.0)
    return -0.5 * (np.sum((data - model) ** 2 * inv_sigma2 - np.log(inv_sigma2)))


ellipse_labels = ['x0_inner', 'y0_inner', 'a_inner', 'b_inner', 'theta_inner', 'amplitude_inner',
                  'x0_outer', 'y0_outer', 'a_outer', 'b_outer', 'theta_outer', 'amplitude_outer', 'x_slope', 'y_slope',
                  'background', 'ln_f']

circle_labels = ['x0_inner', 'y0_inner', 'r_inner', 'amplitude_inner', 'x0_outer', 'y0_outer', 'r_outer',
                 'amplitude_outer', 'background', 'ln_f']
