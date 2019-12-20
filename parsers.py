import argparse

from utils import postprocess
from logger import get_logger

def global_parser():
    # Command line options
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--loglevel', default=['info'], nargs=1,
            help='Console logging level')
    parser.add_argument('--png', action='store_true',
            help='Save also as png')
    parser.add_argument('--axlabel', action='store_true',
            help='Automatic select axis labels')
    parser.add_argument('--legend', action='store_true',
            help='Plot legend')
    parser.add_argument('--overplot', nargs='*', default=None,
            help="Overplot files (extension dependent)")
    parser.add_argument('--oplabel', metavar='LABEL', nargs='*', default=None,
            help="Overplot label for legend")
    parser.add_argument('--opcolor', metavar='COLOR', nargs='*', default=None,
            help="Overplot colors")
    parser.add_argument('config', nargs=1,
            help='Plot configuration file')
    parser.add_argument('plotname', nargs=1,
            help='Plot file name')
    group1 = parser.add_mutually_exclusive_group(required=False)
    group1.add_argument('--shape', metavar=('ROWS','COLS'), nargs=2, 
            default=None, type=int,
            help="Shape of figure tile")
    group1.add_argument('--rows', nargs=1, default=None, type=int,
            help="Number of rows")
    group1.add_argument('--cols', nargs=1, default=None, type=int,
            help="Number of columns")
    group2 = parser.add_mutually_exclusive_group(required=False)
    group2.add_argument('--detailaxlabel', action='store_true',
            help='Automatic detailed axis labels')
    group2.add_argument('--axlabels', nargs='*', default=None,
            help='Axis labels (one per each image)')
    group3 = parser.add_mutually_exclusive_group(required=False)
    group3.add_argument('--overall', action='store_true',
            help='Plot all overplots over each image')
    group3.add_argument('--opmapping', nargs='*', type=int,
            help='Indices of axes where to plot the overplots')
    parser.set_defaults(post=postprocess, cube=None, images=[], overplots=[],
            data=[], logger=get_logger)

    return parser

def common_lines():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--lines', nargs='*', 
            help="Lines to plot (default all in config)")
    parser.add_argument('--lineconfig', nargs=1,
            help="Configuration file of the desired lines")
    group1 = parser.add_mutually_exclusive_group(required=False)
    group1.add_argument('--auto_velshift', action='store_true',
            help="Determine the velocity shift automatically")
    group1.add_argument('--vlsr', nargs='*', type=float,
            help="LSR velocity to shift spectral data in km/s")
    try:
        import astroSource.source as src
        group1.add_argument('--source', action=src.LoadSourcefromConfig,
                help="Source for the vlsr velocity")
        group1.add_argument('--sources', nargs='*',
                action=src.LoadSourcesfromConfig,
                help="Sources for the vlsr velocity")
        group1.add_argument('--atsources', nargs='*',
                action=src.LoadSourcesfromConfig,
                help="Shift velocity using the source position")
    except ImportError:
        pass

    return parser

def common_maps():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--global_levels', action='store_true',
            help='Use same levels for all plots')

    return parser

