#!/usr/bin/python
import argparse

from loaders import *
from parsers import global_parser
# Here import all plotting argparsers in the directory
from plot_channel_maps import plot_channel_maps_parser
from plot_maps import plot_maps_parser
from plot_moments import plot_moments_parser
from plot_from_cube import plot_from_cube_parser
from plot_data import plot_data_parser
from plot_pvmap import plot_pvmaps_parser

def main():
    # Evaluate parser functions
    subpar = {}
    subpar.update(plot_channel_maps_parser())
    subpar.update(plot_maps_parser())
    subpar.update(plot_moments_parser())
    subpar.update(plot_from_cube_parser())
    subpar.update(plot_data_parser())
    subpar.update(plot_pvmaps_parser())

    # Command line options
    parser = argparse.ArgumentParser(description='Data plotting tools.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            parents=[global_parser()], conflict_handler='resolve')
    # Subparsers
    subparsers = parser.add_subparsers()
    # Add subparsers
    for key,(p,h) in subpar.items():
        subparser = subparsers.add_parser(key, parents=[p], help=h)
    
    # Parse arguments
    args = parser.parse_args()
    args.logger(__name__, args)
    for loader in args.loaders:
        args = loader(args)
    args = load_overplot(args)
    figure = args.func(args)
    if figure is not None:
        args.post(figure, args)

if __name__=='__main__':
    main()
