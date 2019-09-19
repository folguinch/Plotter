from builtins import zip

from src.map_plotter import MapsPlotter
from utils import all_mapfig_setup, plot_single_map, get_shape

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
    for ax, img in zip(fig.axes, args.images):
        ax = plot_single_map(ax, fig, img, cen=cen, radius=radius,
                self_levels=True, cbar_orientation=orientation,
                markers=markers)

    fig.auto_config()

    return fig
