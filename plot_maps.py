from builtins import zip

from src.map_plotter import MapsPlotter
from utils import all_mapfig_setup, plot_single_map, get_shape, get_overplots, get_axis_label

def plot_maps(args):
    # Keyword arguments for tile plotter
    opts = {}
    opts['nrows'], opts['ncols'] = get_shape(args, len(args.images),
            default_cols=1)
    assert opts['nrows']*opts['ncols'] >= len(args.images)

    # Setup tile plotter
    fig = MapsPlotter(config=args.config[0], section=args.section[0],
            **opts)
    cen, radius, markers, orientation = all_mapfig_setup(fig)

    # Iterate over images
    for i,(loc, img) in enumerate(zip(fig.axes, args.images)):
        label = fig.get_value('axlabel', None, loc, sep=',')
        label = get_axis_label(args, i, label)
        overplots = get_overplots(args, i)
        ax = plot_single_map(loc, fig, img, contours=overplots, cen=cen, 
                radius=radius, self_contours=args.selflevels, cbar_orientation=orientation,
                markers=markers, axlabel=label)

    fig.auto_config(legend=args.legend)

    return fig
