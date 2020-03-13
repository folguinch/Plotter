import os
import argparse
from configparser import ConfigParser
from builtins import zip, range

import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np

import loaders
from parsers import common_lines
from src.map_plotter import MapsPlotter
import utils as ut 
import utils_image as utimg

def plot_moments_parser():
    # Parser help
    h = "Plot moment maps"

    # Parser
    parser = argparse.ArgumentParser(add_help=False, parents=[common_lines()])
    parser.add_argument('--section', nargs=1, default=['moments'],
            help="Section of the config file")
    parser.add_argument('--moments', nargs='*', type=int,
            help="Number of the moment for each axis")
    parser.add_argument('--momentbase', nargs=1, 
            help="Base name for saving moments")
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument('--imagenames', nargs='*',
            help="List of images to plot")
    group1.add_argument('--cubename', nargs=1,
            help="Cube file name")
    parser.set_defaults(func=plot_moments, linecfg=ConfigParser(),
            loaders=[loaders.load_fits,loaders.load_cube])

    return {'plotmoments': (parser, h)}

def get_axlabel(mom):
    n = mom if (mom <= 20) else (mom % 10)
    suffix = { 1: "st", 2: "nd", 3: "rd" }.get(n, 'th')
    return '%i%s moment' % (mom, suffix)

def plot_moments(args):
    if args.moments is None and len(args.images)>0:
        args.moments = list(range(len(args.images)))
    elif args.moments is None:
        args.moments = [0]

    # Read configuration
    args.logger.info('Reading line configuration file: %s', args.lineconfig[0])
    args.linecfg.read(args.lineconfig[0])

    # Selected lines
    if args.lines:
        lines = args.lines
    else:
        lines = args.linecfg.sections()

    plotbase = os.path.expanduser(args.plotname[0])
    for line in lines:
        # Plot
        if line not in args.linecfg.sections():
            args.logger.warn('%s not in configuration file, skipping..', line)
            continue
        fig = _plot_moments(args, line)

        # Figure title
        if args.linecfg:
            fig.set_title(args.linecfg.get(line, 'title'))
        
        # Save fig
        plotname = os.path.splitext(plotbase)
        plotname = ('.%s.moments' % line).join(plotname)
        args.plotname = [plotname]
        args.logger.info('Saving figure: %s', args.plotname[0])
        ut.postprocess(fig, args)

def _plot_moments(args, line):
    # Keyword arguments for tile plotter
    opts = {}
    opts['nrows'], opts['ncols'] = ut.get_shape(args, len(args.moments),
            default_cols=len(args.moments))

    # Setup tile plotter
    fig = MapsPlotter(config=args.config[0], section=args.section[0],
            **opts)
    cen, radius, markers, orientation = ut.all_mapfig_setup(fig)

    # Iterate over images
    peak = None
    for i,(loc,mom) in enumerate(zip(fig.axes, args.moments)):
        # Get image
        if args.images is not None:
            img = args.images[i]
        else:
            if args.momentbase:
                filename = args.momentbase[0]+'.%s.moment%i.fits' % (line,mom)
            else:
                filename = None
            img = utimg.moment(args.cube, mom, args.linecfg[line],
                    filename=filename)

        # Overplots
        overplots = ut.get_overplots(args, i)

        # Labels
        label = fig.get_value('axlabel', get_axlabel(mom), loc, sep=',')
        label = ut.get_axis_label(args, i, label)

        # Recenter
        ceni = cen[0] if len(cen)==1 else cen[i]
        radiusi = radius[0] if len(radius)==1 else radius[i]
        if mom != 1:
            ax = ut.plot_single_map(loc, fig, img, args.logger, cen=ceni,
                    radius=radiusi, self_contours=mom==0 and overplots is None, 
                    cbar_orientation=orientation, markers=markers,
                    contours=overplots,
                    axlabel=label)
            if mom==0:
                peak = np.unravel_index(np.nanargmax(np.squeeze(img.data)),
                        np.squeeze(img.data).shape)

        else:
            if args.auto_velshift and peak is not None:
                vel_shift = np.squeeze(img.data)[peak]
            elif args.vlsr is not None:
                if len(args.vlsr)==1:
                    vel_shift = args.vlsr[0]
                else:
                    vel_shift = args.vlsr[i]
            elif args.sources is not None:
                vel_shift = args.sources[i].get_quantity('vlsr')
                vel_shift = vel_shift.to(u.km/u.s).value
            elif args.source is not None:
                vel_shift = args.source.get_quantity('vlsr')
                vel_shift = vel_shift.to(u.km/u.s).value
            elif args.atsources is not None:
                wcs = WCS(img.header, naxis=['longitude','latitude'])
                pos = args.atsources[i].position
                pix = wcs.all_world2pix([[pos.ra.deg,pos.dec.deg]],0)[0]
                pix = map(int, pix[::-1])
                vel_shift = np.squeeze(img.data)[pix[0],pix[1]]
            else:
                args.logger.warn('Could not find vlsr')
                vel_shift = 0.
            args.logger.info('Velocity shift: %f', vel_shift)
            img.data = img.data - vel_shift

            with plt.style.context('bwr'):
                ax = ut.plot_single_map(loc, fig, img, args.logger, cen=ceni, 
                        radius=radiusi,
                        self_levels=True, cbar_orientation=orientation,
                        markers=markers, dtype='velocity',
                        skip_marker_label=args.moments.index(mom)>0,
                        axlabel=label)

    fig.auto_plot()
    fig.auto_config()

    return fig
