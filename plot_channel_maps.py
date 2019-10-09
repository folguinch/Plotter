import argparse
from builtins import range, zip

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from loaders import load_cube, load_fits
from parsers import common_lines
from src.map_plotter import MapsPlotter
from src.utils import auto_vminmax, auto_levels
from utils import all_mapfig_setup, plot_single_map, get_shape, plot_contours

def plot_channel_maps_parser():
    # Parser help
    h = "Plot channel maps"
    parser = argparse.ArgumentParser(add_help=False, parents=[common_lines()])
    parser.add_argument('--section', nargs=1, default=['channel_maps'],
            help="Section of the config file")
    parser.add_argument('--continuum', dest='imagenames', nargs=1,
            default=[], metavar='CONTINUUM', 
            help="Background image file for color plot")
    parser.add_argument('--every', metavar='STEP', nargs=1, default=[1],
            type=int,
            help="Step between channels for range")
    parser.add_argument('cubename', nargs=1,
            help="Cube file name")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--chanran', metavar=('INITIAL','FINAL'), nargs=2, 
            default=None, type=int,
            help="Channel range to plot (final not included)")
    group.add_argument('--chans', metavar='CHAN', nargs='*', 
            default=None, type=int,
            help="List of channels to plot")
    parser.set_defaults(func=plot_channel_maps, loaders=[load_cube, load_fits])

    return {'chanmap': (parser, h)}

def plot_channel_maps(args):
    # Number of channels
    if args.chanran:
        ichans = range(*(args.chanran+args.every))
    elif args.chans:
        ichans = args.chans
    nchans = len(ichans)
    args.logger.info('Number of channels = %i', nchans)
    assert nchans>0

    # Keyword arguments for tile plotter
    opts = {}
    opts['nrows'], opts['ncols'] = get_shape(args, nchans)
    assert opts['nrows']*opts['ncols'] >= nchans
    args.logger.info('Rows, columns = %i, %i', opts['nrows'], opts['ncols'])

    # Velocity axis
    args.cube = args.cube.with_spectral_unit(u.km/u.s,
            velocity_convention='radio')
    vel = args.cube.spectral_axis
    ind = np.argsort(vel[ichans])
    ichans = list(np.array(ichans)[ind])
    if args.auto_velshift:
        ind = len(ichans)//2
        vel = vel - vel[ichans][ind]
    elif args.vlsr is not None:
        vel = vel - args.vlsr * u.km/u.s
    else:
        pass

    # Setup tile plotter
    cubewcs = args.cube.wcs.sub(('longitude','latitude'))
    args.logger.debug('Initializing figure')
    fig = MapsPlotter(config=args.config[0], section=args.section[0],
            projection=cubewcs, **opts)
    cen, radius, markers, orientation = all_mapfig_setup(fig)

    # Get global vmin and vmax
    if len(args.images) == 1:
        args.logger.debug('Using input background image')
        img = args.images[0]
        data = np.squeeze(img.data)
    else:
        data = args.cube.unmasked_data[ichans,:,:].value
    vmin, vmax = auto_vminmax(data)
    vmin = fig.config.getfloat('vmin', fallback=vmin)
    vmax = fig.config.getfloat('vmax', fallback=vmax)
    args.logger.info('Cube vmin, vmax = %.3e, %.3e', vmin, vmax)

    # Levels
    if len(args.images) == 1:
        vminline = fig.config.getfloat('vminline', fallback=None)
        vmaxline = fig.config.getfloat('vmaxline', fallback=None)
    else:
        vminline = vmin
        vmaxline = vmax
    levels = auto_levels(args.cube.unmasked_data[ichans,:,:].value,
            n=fig.config.getint('nlevels', fallback=10), 
            stretch=fig.config.get('stretch', fallback='log'),
            vmin=vminline, vmax=vmaxline)
    args.logger.info('Global levels: %r', levels)

    # Plot
    for loc,i in zip(fig.axes, ichans):
        cbar = ichans.index(i)==0
        #ax = fig.get_mapper(ax, include_cbar=cbar, vmin=vmin, vmax=vmax, a=a)
        bmaj = args.cube.beams[i].major.to(u.deg).value
        bmin = args.cube.beams[i].minor.to(u.deg).value
        bpa = args.cube.beams[i].pa.to(u.deg).value

        if len(args.images) == 0:
            self_contours = True
            img = fits.PrimaryHDU(args.cube.unmasked_data[i,:,:].value,
                    header=cubewcs.to_header())
            img.header['BMAJ'] = bmaj
            img.header['BMIN'] = bmin
            img.header['BPA'] = bpa
            img.header['BUNIT'] = args.cube.header['BUNIT']
        else:
            self_contours = False
            contours = fits.PrimaryHDU(args.cube.unmasked_data[i,:,:].value,
                    header=cubewcs.to_header())
            contours.header['BMAJ'] = bmaj
            contours.header['BMIN'] = bmin
            contours.header['BPA'] = bpa
            optscont = {'colors': 'w'}
            contours = [(contours, optscont)]

        label = '%.2f %s' % (vel[i].value, vel.unit.to_string('latex_inline'))

        ax = plot_single_map(loc, fig, img, args.logger, cen=cen, 
                radius=radius, self_contours=self_contours, levels=levels, 
                cbar_orientation=orientation if cbar else None,
                markers=markers, axlabel=label, vmin=vmin, vmax=vmax, 
                include_cbar=cbar)
        if not self_contours:
            plot_contours(ax, contours, levels, zorder=2)

    fig.auto_config()

    return fig

