import os, argparse
from configparser import ConfigParser

import astropy.units as u

from loaders import load_cube, load_fits
from utils import get_quantity

def plot_from_cube_parser():
    # Parser help
    h = "Plot data from cube (moments 0, 1, and channel maps)"

    # Parser
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--continuum', dest='imagenames', nargs=1,
            default=[], metavar='CONTINUUM', 
            help="Background image file for color plot")
    parser.add_argument('--every', metavar='STEP', nargs=1, default=[1],
            type=int,
            help="Step between channels for range")
    parser.add_argument('--save_moments', nargs=1, default=[None],
            help="Moments base file name")
    parser.add_argument('--line_list', nargs='*',
            help="Lines to plot (default all in config)")
    parser.add_argument('lineconfig', nargs=1,
            help="Configuration file of the desired lines")
    parser.add_argument('cubename', nargs=1,
            help="Cube file name")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--chanran', metavar=('INITIAL','FINAL'), nargs=2, 
            default=None, type=int,
            help="Channel range to plot (final not included)")
    group.add_argument('--chans', metavar='CHAN', nargs='*', 
            default=None, type=int,
            help="List of channels to plot")
    parser.set_defaults(func=plot_from_cube, lines=ConfigParser(),
            loaders=[load_cube, load_fits])

    return {'fromcube': (parser, h)}

def plot_from_cube(args):
    # Read configuration
    args.logger.info('Reading line configuration file: %s', args.lineconfig)
    args.lines.read(args.lineconfig)

    # Base arguments
    bargs = get_base_args(args)

    # Base plotname
    plotname = os.path.splitext(args.plotname[0])

    for line in args.lines.sections():
        if args.line_list is not None and line not in args.line_list:
            continue
        args.logger.info('Plotting line: %s', line)
        config = args.lines[line]

        # Get arguments
        subcube, momargs, chanargs = process_cube(args.cube, config, 
                args.logger, args.save_moments[0], args.every) 

        # Plot
        momargs = bargs + ['_%s_moments'.join(plotname) % line, 'plotmoments'] + momargs
        chanargs = bargs + [ '_%s_chanmap'.join(plotname) % line, 'chanmap'] \
                + chanargs + [ 'dummy_cube'] 
        plot_all(args, subcube, momargs, chanargs)

    return None

def process_cube(cube, config, logger, filename=None, every=[1]):
    # New arguments
    momargs = ['--moments', '0', '1']
    chanargs = []

    # Data needed from config
    logger.info('Extracting sub-cube:')
    vlsr = get_quantity(config['vlsr'])
    restfreq = get_quantity(config['freq'])
    chmin, chmax = map(int, config['chanran'].split())
    logger.info('Source v_LSR = %s', vlsr)
    logger.info('Rest frequency = %s', restfreq)
    logger.info('Channel range = %i, %i', chmin, chmax)
    momargs += ['--vlsr', str(vlsr.to(u.km/u.s).value)]
    chanargs += ['--vlsr', str(vlsr.to(u.km/u.s).value)]

    # Convert velocity
    cube = cube.with_spectral_unit(u.km/u.s,
            velocity_convention='radio', rest_value=restfreq)
    vel = cube.spectral_axis
    logger.info('Velocity range = %s, %s',vel[chmin], vel[chmax])
    subcube = cube[chmin:chmax+1,:,:]
    chanargs += ['--chanran', '0', '%i' % (len(subcube.spectral_axis),)]
    chanargs += ['--every'] + map(str, every)

    # Moment 0
    mom0 = subcube.moment(order=0)

    # Mask data for moment 1
    if 'low' in config:
        low = get_quantity(config['low'])
        mask = subcube >= low
    else:
        mask = subcube.mask
    if 'up' in config:
        up = get_quantity(config['up'])
        mask = mask & (subcube<=up)
    subcube = subcube.with_mask(mask)
    mom1 = subcube.moment(order=1)

    # Save moments
    if filename:
        m0name = filename+'.%s.moment0.fits' % config.name
        m1name = filename+'.%s.moment1.fits' % config.name
        mom0.write(m0name, overwrite=True)
        mom1.write(m1name, overwrite=True)
        momargs += ['--imagenames', m0name, m1name]
    else:
        # Create a temp file
        raise NotImplementedError

    return subcube, momargs, chanargs

def plot_all(args, cube, momargs, chanargs):
    # Get argparsers
    momparser, chanparser = get_parsers()
    margs = momparser.parse_args(momargs)
    cargs = chanparser.parse_args(chanargs)
    cargs.cube = cube

    # Plot moments
    margs.logger = args.logger
    for loader in margs.loaders:
        margs = loader(margs)
    margs.shape = [1, 2]
    momfig = margs.func(margs)
    margs.post(momfig, margs)
    
    # Plot chanmaps
    cargs.logger = args.logger
    cargs.images = args.images
    chanfig = cargs.func(cargs)
    cargs.post(chanfig, cargs)

    return momfig, chanfig

def get_base_args(args):
    # Triggers
    bargs = []
    if args.png:
        bargs += ['--png']
    if args.axlabel:
        bargs += ['--axlabel']
    if args.legend:
        bargs += ['--legend']
    if args.detailaxlabel:
        bargs += ['--detailaxlabel']

    # Mutually exclusive
    if args.shape:
        bargs += ['--shape'] + map(str, args.shape)

    bargs += ['--loglevel'] + args.loglevel
    bargs += args.config

    return bargs

def get_parsers():
    from parsers import global_parser
    from plot_moments import plot_moments_parser
    from plot_channel_maps import plot_channel_maps_parser
    
    momparser = argparse.ArgumentParser(add_help=False,
            parents=[global_parser()])
    subparsers1 = momparser.add_subparsers()
    for key, val in plot_moments_parser().items():
        subparser1 = subparsers1.add_parser(key, parents=[val[0]],
                help=val[1])

    chanparser = argparse.ArgumentParser(add_help=False,
            parents=[global_parser()])
    subparsers2 = chanparser.add_subparsers()
    for key, val in plot_channel_maps_parser().items():
        subparser2 = subparsers2.add_parser(key, parents=[val[0]],
                help=val[1])

    return momparser, chanparser
