import argparse

from loaders import load_fits
from parsers import common_maps
import plot_maps as pltmaps

def plot_pvmaps_parser():
    # Parser help
    h = "Plot pv maps"

    # Parser
    parser = argparse.ArgumentParser(add_help=False, parents=[common_maps()])
    parser.add_argument('--section', nargs=1, default=['maps_plots'],
            help="Section of the config file")
    parser.add_argument('imagenames', nargs='*',
            help="List of images to plot")
    parser.set_defaults(func=plot_pvmaps, loaders=[load_fits])

    return {'plotpvmaps': (parser, h)}

def plot_pvmaps(args):
    return pltmaps.plot_maps(args, dtype='pvmap')
