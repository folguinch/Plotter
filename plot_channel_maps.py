import os
import argparse
from builtins import range, zip
from configparser import ConfigParser

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

import utils as ut
from loaders import load_cube, load_fits
from parsers import common_lines
from src.map_plotter import MapsPlotter
from src.maths import quick_rms
from src.utils import auto_vminmax, auto_levels

def plot_channel_maps_parser():
    # Parser help
    h = "Plot channel maps"
    parser = argparse.ArgumentParser(add_help=False, parents=[common_lines()])
    parser.add_argument('--section', nargs=1, default=['channel_maps'],
            help="Section of the config file")
    parser.add_argument('--continuum', dest='imagenames', nargs=1,
            default=[], metavar='CONTINUUM', 
            help="Background image file for color plot")
    parser.add_argument('--every', metavar='STEP', nargs=1, default=[None],
            type=int,
            help="Step between channels for range")
    parser.add_argument('--freq', nargs=1, default=[None],
            help="Line rest frequency in GHz")
    parser.add_argument('cubename', nargs=1,
            help="Cube file name")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--chanran', metavar=('INITIAL','FINAL'), nargs=2, 
            default=None, type=int,
            help="Channel range to plot (final not included)")
    group.add_argument('--chans', metavar='CHAN', nargs='*', 
            default=None, type=int,
            help="List of channels to plot")
    parser.set_defaults(func=plot_channel_maps, loaders=[load_cube, load_fits])

    return {'chanmap': (parser, h)}

def get_channel_indices(args, nrms=5.):
    # Initial moments
    if args.chanran:
        aux1, aux2 = args.chanran
        ind = range(aux1, aux2+1)
    elif args.chans:
        ichans = args.chans
        nrows, ncols = ut.get_shape(args, len(ichans))
        args.logger.info('Number of channels = %i', len(ichans))
        return ichans, nrows, ncols
    else:
        raise ValueError('No channels were selected')

    # Get rms
    rms = quick_rms(args.cube.unmasked_data[ind,:,:].value) * args.cube.unit
    args.logger.info('Preliminary rms: %s', rms)
    mask = args.cube >= nrms*rms
    mask = np.squeeze(mask.include())
    nvalid = np.sum(mask, axis=(-1,-2))

    # Final channel limits
    for i, nval in zip(ind, nvalid[ind]):
        if nval==0:
            args.logger.info('Channel %i below threshold', i)
            aux1 = i+1
        else:
            break
    for i, nval in zip(ind[::-1],nvalid[ind][::-1]):
        if nval==0:
            args.logger.info('Channel %i below threshold', i)
            aux2 = i-1
        else:
            break
    args.logger.info('Final channel range: %i, %i, %i', aux1, aux2,
            args.every[0])

    # Final results
    ichans = range(aux1, aux2+1, args.every[0])
    nrows, ncols = ut.get_shape(args, len(ichans), minimize=True)
    args.logger.info('Number of channels = %i', len(ichans))
    args.logger.info('Figure size = %i', len(ichans))
    while nrows*ncols<len(ichans):
        i0 = ichans[0]
        i1 = ichans[-1]
        if nvalid[i0]<=nvalid[i1]:
            args.logger.info('Dropping channel %i to match shape', i0)
            ichans = ichans[1:]
        else:
            args.logger.info('Dropping channel %i to match shape', i1)
            ichans = ichans[:-1]

    return ichans, nrows, ncols

def plot_channel_maps(args):
    # Line config
    if args.lineconfig:
        linecfg = ut.read_config(args.lineconfig[0])
    else:
        linecfg = None

    # Selected lines
    if args.lines:
        if linecfg is None:
            lines = None
        else:
            lines = args.lines
    elif linecfg is not None:
        lines = linecfg.sections()
    else:
        lines = None

    # Iterate over lines
    if linecfg is None:
        return _plot_single_channel_maps(args)
    else:
        plotbase = os.path.expanduser(args.plotname[0])
        for line in lines:
            # Assign values for this iteration
            if args.every[0] is None:
                args.every = [linecfg.getint(line, 'every', fallback=1)]
            args.chanran = map(int, linecfg.get(line, 'chanran').split())
            freq = linecfg.get(line, 'freq').split()
            freq = float(freq[0]) * u.Unit(freq[1])
            args.freq = [freq.to(u.GHz).value]

            # Plot
            fig = _plot_single_channel_maps(args)

            # Figure title
            fig.set_title(linecfg.get(line, 'title', fallback=''), ha='left')
            
            # Save fig
            plotname = os.path.splitext(plotbase)
            plotname = ('.%s.chanmap' % line).join(plotname)
            args.plotname = [plotname]
            args.logger.info('Saving figure: %s', args.plotname[0])
            ut.postprocess(fig, args)

    return None

def _plot_single_channel_maps(args):
    # Every step
    if args.every[0] is None:
        args.every = [1]
    
    # Keyword arguments for tile plotter
    opts = {}
    ichans, opts['nrows'], opts['ncols'] = get_channel_indices(args)
    assert opts['nrows']*opts['ncols'] >= len(ichans)
    args.logger.info('Rows, columns = %i, %i', opts['nrows'], opts['ncols'])

    # Velocity axis
    if args.freq[0] is not None:
        args.logger.info('Using rest frequency: %.3f GHz', args.freq[0])
        args.cube = args.cube.with_spectral_unit(u.GHz,
                velocity_convention='radio')
        args.cube = args.cube.with_spectral_unit(u.km/u.s,
                velocity_convention='radio', rest_value=args.freq[0]*u.GHz)
    else:
        args.cube = args.cube.with_spectral_unit(u.km/u.s,
                velocity_convention='radio')
    vel = args.cube.spectral_axis
    ind = np.argsort(vel[ichans])
    ichans = list(np.array(ichans)[ind])
    if args.auto_velshift:
        ind = len(ichans)//2
        vel_shift = vel[ichans][ind]
    elif args.vlsr is not None:
        vel_shift = args.vlsr * u.km/u.s
    elif args.sources is not None:
        args.logger.warn('Using first source for velocity shift')
        vel_shift = args.sources[0].get_quantity('vlsr')
        vel_shift = vel_shift.to(u.km/u.s)
    elif args.source is not None:
        vel_shift = args.source.get_quantity('vlsr')
        vel_shift = vel_shift.to(u.km/u.s)
    elif args.atsources is not None:
        wcs = WCS(img.header, naxis=['longitude','latitude'])
        pos = args.atsources[i].position
        pix = wcs.all_world2pix([[pos.ra.deg,pos.dec.deg]],0)[0]
        pix = map(int, pix[::-1])
        spec = np.squeeze(args.cube.unmasked_data)[:,pix[0],pix[1]][ichans]
        maxind = np.nanargmax(spec)
        vel_shift = vel[ichans][maxind]
    else:
        args.logger.warn('Could not find velocity shift')
        vel_shift = 0.
    args.logger.info('Shifting velocity axis by: %s', vel_shift)
    vel = vel - vel_shift

    # Setup tile plotter
    cubewcs = args.cube.wcs.sub(('longitude','latitude'))
    args.logger.debug('Initializing figure')
    fig = MapsPlotter(config=args.config[0], section=args.section[0],
            projection=cubewcs, **opts)
    cen, radius, markers, orientation = ut.all_mapfig_setup(fig)

    # Data
    if args.images is None:
        data = args.cube.unmasked_data[ichans,:,:].value
    elif len(args.images) == 1:
        args.logger.debug('Using input background image')
        img = args.images[0]
        data = np.squeeze(img.data)
    else:
        raise NotImplementedError

    # rms, nsigma
    try:
        rms = float(fig.get_value('rms'))
    except TypeError:
        rms = None
    nsigma = float(fig.get_value('nsigma', default=5.))

    # Get global vmin and vmax
    vmin, vmax = auto_vminmax(data)
    vmin = fig.config.getfloat('vmin', fallback=vmin)
    vmax = fig.config.getfloat('vmax', fallback=vmax)
    args.logger.info('Color map vmin, vmax = %.3e, %.3e', vmin, vmax)

    # Levels
    levels = auto_levels(args.cube.unmasked_data[ichans,:,:].value,
            rms=rms, nsigma=nsigma)
    args.logger.info('Global line levels: %r', levels)

    # Plot
    for loc,i in zip(fig.axes, ichans):
        cbar = ichans.index(i)==0
        #ax = fig.get_mapper(ax, include_cbar=cbar, vmin=vmin, vmax=vmax, a=a)
        bmaj = args.cube.beams[i].major.to(u.deg).value
        bmin = args.cube.beams[i].minor.to(u.deg).value
        bpa = args.cube.beams[i].pa.to(u.deg).value
        overplots = ut.get_overplots(args, i)

        if args.images is None or len(args.images) == 0:
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
            contcolor = fig.get_value('contour_color', default='w', ax=loc)
            optscont = {'colors': contcolor}
            contours = [(contours, optscont)]

        label = '%.2f %s' % (vel[i].value, vel.unit.to_string('latex_inline'))

        ax = ut.plot_single_map(loc, fig, img, args.logger, cen=cen[0], 
                radius=radius[0], self_contours=self_contours, levels=levels, 
                cbar_orientation=orientation if cbar else None,
                markers=markers, axlabel=label, vmin=vmin, vmax=vmax, 
                include_cbar=cbar, contours=overplots)
        if not self_contours:
            ut.plot_contours(ax, contours, levels, zorder=2)

    fig.auto_config()

    return fig

