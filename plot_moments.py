from builtins import zip, range

import numpy as np
import matplotlib.pyplot as plt

from src.map_plotter import MapsPlotter
from utils import all_mapfig_setup, plot_single_map, get_shape, get_axis_label

def get_axlabel(mom):
    n = mom if (mom <= 20) else (mom % 10)
    suffix = { 1: "st", 2: "nd", 3: "rd" }.get(n, 'th')
    return '%i%s moment' % (mom, suffix)

def plot_moments(args):
    if args.moments is None:
        args.moments = list(range(len(args.images)))
    assert len(args.images)==len(args.moments)

    # Keyword arguments for tile plotter
    opts = {}
    opts['nrows'], opts['ncols'] = get_shape(args, len(args.images),
            default_cols=len(args.images))
    assert opts['nrows']*opts['ncols'] >= len(args.images)

    # Setup tile plotter
    fig = MapsPlotter(config=args.config[0], section=args.section[0],
            **opts)
    cen, radius, markers, orientation = all_mapfig_setup(fig)

    # Iterate over images
    peak = None
    for i,(loc,mom,img) in enumerate(zip(fig.axes, args.moments, args.images)):
        label = fig.get_value('axlabel', get_axlabel(mom), loc, sep=',')
        label = get_axis_label(args, i, label)
        if mom != 1:
            ax = plot_single_map(loc, fig, img, cen=cen, radius=radius,
                    self_contours=True, cbar_orientation=orientation,
                    markers=markers, axlabel=label)
            if mom==0:
                peak = np.unravel_index(np.nanargmax(np.squeeze(img.data)),
                        np.squeeze(img.data).shape)

        else:
            if args.auto_velshift and peak is not None:
                vel_shift = np.squeeze(img.data)[peak]
            elif 'vel_shift' in fig.config:
                vel_shift = fig.config.getfloat('vel_shift')
            else:
                vel_shift = 0.
            img.data = img.data - vel_shift

            with plt.style.context('bwr'):
                ax = plot_single_map(loc, fig, img, cen=cen, radius=radius,
                        self_levels=True, cbar_orientation=orientation,
                        markers=markers, dtype='velocity',
                        skip_marker_label=args.moments.index(mom)>0,
                        axlabel=label)

    fig.auto_config()

    return fig
