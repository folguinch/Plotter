import argparse
from builtins import zip

from loaders import load_fits
from parsers import common_maps
from src.map_plotter import MapsPlotter
from utils import all_mapfig_setup, plot_single_map, get_shape, get_overplots, get_axis_label

def plot_maps_parser():
    # Parser help
    h = "Plot FITS maps"

    # Parser
    parser = argparse.ArgumentParser(add_help=False, parents=[common_maps()])
    parser.add_argument('--section', nargs=1, default=['maps_plots'],
            help="Section of the config file")
    parser.add_argument('imagenames', nargs='*',
            help="List of images to plot")
    parser.set_defaults(func=plot_maps, loaders=[load_fits])

    return {'plotmaps': (parser, h)}

def plot_maps(args, dtype='intensity'):
    # Keyword arguments for tile plotter
    opts = {}
    opts['nrows'], opts['ncols'] = get_shape(args, len(args.images),
            default_cols=1)
    assert opts['nrows']*opts['ncols'] >= len(args.images)
    args.logger.info('Rows, columns = %i, %i', opts['nrows'], opts['ncols'])

    # Setup tile plotter
    args.logger.debug('Initializing figure')
    fig = MapsPlotter(config=args.config[0], section=args.section[0],
            **opts)
    cen, radius, markers, artists, orientation = all_mapfig_setup(fig)

    # Iterate over images
    for i,(loc, img) in enumerate(zip(fig.axes, args.images)):
        label = fig.get_value('axlabel', '', loc, sep=',')
        label = get_axis_label(args, i, label)
        overplots = get_overplots(args, i)

        # Global color bar
        if args.global_color and i>0:
            orientation = None
            cbar = False
        else:
            cbar = True

        # Recenter
        ceni = cen[0] if len(cen)==1 else cen[i]
        radiusi = radius[0] if len(radius)==1 else radius[i]

        # Plot
        ax = plot_single_map(loc, fig, img, args.logger, contours=overplots, 
                cen=ceni, radius=radiusi, self_contours=args.selflevels, 
                cbar_orientation=orientation, markers=markers, axlabel=label,
                dtype=dtype, nsigmalevel=args.nsigmalevel[0],
                include_cbar=cbar, artists=artists)

    fig.auto_config(legend=args.legend, dtype=dtype)

    return fig
