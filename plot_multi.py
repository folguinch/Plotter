import argparse
from builtins import zip

import loaders
from parsers import common_maps
from src.multi_plotter import MultiPlotter
from utils import all_mapfig_setup, plot_single_map, get_shape, get_overplots, get_axis_label

def plot_multi_parser():
    # Parser help
    h = "Plot different type of plots"

    # Parser
    parser = argparse.ArgumentParser(add_help=False, parents=[common_maps()])
    parser.add_argument('--section', nargs=1, default=['maps_plots'],
            help="Section of the config file")
    parser.set_defaults(func=plot_multi, loaders=[])

    return {'multiplot': (parser, h)}

def plot_multi(args):
    # Create plot
    fig = MultiPlotter(config=args.config[0])

    # Define data loaders
    load = get_loaders(fig._config, args)

    # Plot
    fig.auto_plot(load)

    return fig

def get_loaders(cfg, args):
    aux = {}
    for section in cfg.sections():
        ctype = cfg[section]['type']
        if ctype in ['map', 'contour_map']:
            aux[section] = loaders.multi_map
        elif ctype == 'spectrum':
            aux[section] = loaders.multi_spectrum
        elif ctype == 'moment':
            aux[section] = loaders.multi_moment_map
        else:
            raise NotImplementedError('No type: %s' % ctype)

    return aux

