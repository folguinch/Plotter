#!/usr/bin/python
import os, argparse

from loaders import *
from utils import postprocess
# Here import all plotting functions in the directory
from plot_channel_maps import plot_channel_maps

def main():
    # Command line options
    parser = argparse.ArgumentParser(description='Data plotting tools.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    parser.add_argument('--png', action='store_true',
            help='Save also as png')
    parser.add_argument('config', nargs=1,
            help='Plot configuration file')
    parser.add_argument('plotname', nargs=1,
            help='Plot file name')
    parser.set_defaults(post=postprocess)
    group1 = parser.add_mutually_exclusive_group(required=False)
    group1.add_argument('--shape', metavar=('ROWS','COLS'), nargs=2, 
            default=None, type=int,
            help="Shape of figure tile")
    group1.add_argument('--rows', nargs=1, default=None, type=int,
            help="Number of rows")
    group1.add_argument('--cols', nargs=1, default=None, type=int,
            help="Number of columns")
    # Subparsers
    chanmaps = subparsers.add_parser('chanmap', 
            help="Plot channel maps")
    chanmaps.add_argument('--section', nargs=1, default=['channel_maps'],
            help="Section of the config file")
    chanmaps.add_argument('--imagename', nargs=1, default=None,
            help="Background image file for color plot")
    chanmaps.add_argument('cubename', nargs=1,
            help="Cube file name")
    group2 = chanmaps.add_mutually_exclusive_group(required=True)
    group2.add_argument('--chanran', metavar=('INITIAL','FINAL'), nargs=2, 
            default=None, type=int,
            help="Channel range to plot (final not included)")
    group2.add_argument('--chans', metavar='CHAN', nargs='*', 
            default=None, type=int,
            help="List of channels to plot")
    chanmaps.set_defaults(func=plot_channel_maps, cube=None, image=None,
            loaders=[load_cube, load_fits])
    args = parser.parse_args()
    for loader in args.loaders:
        args = loader(args)
    figure = args.func(args)
    args.post(figure, args)

if __name__=='__main__':
    main()
