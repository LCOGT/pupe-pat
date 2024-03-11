"""
pupepat.plot - Plotting utils for visualizing the results of the pupil plate fits

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0

February 2018
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from pupepat.ellipse import generate_ellipse, calculate_donut_offset
import numpy as np
import os


def plot_best_fit_ellipse(output_filename, data_cutout, best_fit_model, header):
    """
    Plot the best fit elliptical annulus on top of the data.

    :param output_filename: Filename to save the plot
    :param data_cutout: numpy array with the cutout of the data used for the fit
    :param best_fit_model: astropy.modeling.Model corresponding to the best fit ellipse.
    :param header: dict-like with M2PITCH, M2ROLL, and FOCDMD keywords
    """
    pyplot.clf()
    pyplot.imshow(data_cutout, origin='lower',
                  vmin=np.percentile(data_cutout, 40),  # median is 50
                  vmax=np.percentile(data_cutout, 99))  # filter hottest pixels
    pyplot.colorbar()
    x_inner, y_inner = generate_ellipse(best_fit_model.x0_inner, best_fit_model.y0_inner, best_fit_model.a_inner,
                                        best_fit_model.b_inner, best_fit_model.theta_inner)
    pyplot.plot(x_inner, y_inner, color='cyan')
    x_outer, y_outer = generate_ellipse(best_fit_model.x0_outer, best_fit_model.y0_outer, best_fit_model.a_outer,
                                        best_fit_model.b_outer, best_fit_model.theta_outer)
    pyplot.plot(x_outer, y_outer, color='cyan')
    pyplot.title(os.path.basename(os.path.splitext(output_filename)[0]))
    pyplot.text(data_cutout.shape[1] * 0.01, data_cutout.shape[0] * 0.95,
                'Inner: {x: 0.3f}, {y: 0.3f}'.format(x=best_fit_model.x0_inner, y=best_fit_model.y0_inner),
                color='white')
    pyplot.text(data_cutout.shape[1] * 0.01, data_cutout.shape[0] * 0.915,
                'Outer: {x: 0.3f}, {y: 0.3f}'.format(x=best_fit_model.x0_outer, y=best_fit_model.y0_outer),
                color='white')

    # If using an sbig camera on a 1m, calculate the M1 offset, otherwise give offset in pixels
    if 'kb' in header['INSTRUME'].lower() and '1m0' in header['TELID'].lower():
        pyplot.text(data_cutout.shape[1] * 0.01, data_cutout.shape[0] * 0.88,
                    'M1 Offset (mm): {x: 0.3f}, {y: 0.3f}'.format(x=calculate_donut_offset(best_fit_model.x0_inner,
                                                                                           best_fit_model.x0_outer,
                                                                                           sbig=True),
                                                                  y=calculate_donut_offset(best_fit_model.y0_inner,
                                                                                           best_fit_model.y0_outer,
                                                                                           sbig=True)),
                    color='white')
    else:
        pyplot.text(data_cutout.shape[1] * 0.01, data_cutout.shape[0] * 0.88,
                    'I-O Offset (pix): {x: 0.3f}, {y: 0.3f}'.format(x=calculate_donut_offset(best_fit_model.x0_inner,
                                                                                             best_fit_model.x0_outer),
                                                                    y=calculate_donut_offset(best_fit_model.y0_inner,
                                                                                             best_fit_model.y0_outer)),
                    color='white')
    pyplot.text(data_cutout.shape[1] * 0.6, data_cutout.shape[0] * 0.95,
                'M2 Pitch: {pitch: 0.3f}"'.format(pitch=header['M2PITCH']), color='white')
    pyplot.text(data_cutout.shape[1] * 0.6, data_cutout.shape[0] * 0.915,
                'M2 Roll: {roll: 0.3f}"'.format(roll=header['M2ROLL']), color='white')
    pyplot.text(data_cutout.shape[1] * 0.6, data_cutout.shape[0] * 0.88,
                'FOCDMD: {focdmd: 0.3f} mm'.format(focdmd=header['FOCDMD']), color='white')
    pyplot.tight_layout()
    pyplot.savefig(output_filename)


def plot_quiver(data, focus_surface, output_plot, dpi=200):
    pyplot.clf()
    pyplot.grid(True)
    pyplot.gcf().set_size_inches(11, 10, forward=True)
    # Make a plotting range a little larger than the M2 pitch and roll grid
    x_grid, y_grid = np.meshgrid(1.1 * np.linspace(min(data['M2ROLL']), max(data['M2ROLL']), 10 * dpi),
                                 1.1 * np.linspace(min(data['M2PITCH']), max(data['M2PITCH']), 10 * dpi))
    center_offset_surface = focus_surface[0].eval(x_grid, y_grid) ** 2.0
    center_offset_surface += focus_surface[1].eval(x_grid, y_grid) ** 2.0
    center_offset_surface **= 0.5
    pyplot.contourf(x_grid, y_grid, center_offset_surface, 15)
    pyplot.quiver(data['M2ROLL'], data['M2PITCH'], data['x0_inner'] - data['x0_outer'],
                  data['y0_inner'] - data['y0_outer'], headlength=3, headaxislength=3.0)
    best_roll = x_grid.ravel()[np.argmin(center_offset_surface)]
    best_pitch = y_grid.ravel()[np.argmin(center_offset_surface)]
    pyplot.scatter(best_roll, best_pitch, marker='X', s=200, c='r', lw=0,
                   label='ROLL={roll:+0.0f}, PITCH={pitch:+0.0f}'.format(roll=best_roll, pitch=best_pitch))
    pyplot.tick_params(axis='both', which='major', labelsize=18)
    pyplot.title("Donut Decenter (Inner - Outer)", fontsize=22)
    pyplot.xlabel("M2 Roll tilt (arcsec)", fontsize=20)
    pyplot.ylabel("M2 Pitch tilt (arcsec)", fontsize=20)
    pyplot.legend(loc='upper right', framealpha=0.9, fontsize=18)
    pyplot.tight_layout()
    pyplot.savefig(output_plot)
