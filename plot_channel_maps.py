from builtins import range, zip

import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from src.map_plotter import MapsPlotter
from src.utils import auto_vminmax, auto_levels
from utils import get_shape

def plot_channel_maps(args):
    # Number of channels
    if args.chanran:
        ichans = range(*args.chanran)
    elif args.chans:
        ichans = args.chans
    nchans = len(ichans)
    assert nchans>0

    # Keyword arguments for tile plotter
    opts = {}
    opts['nrows'], opts['ncols'] = get_shape(args, nchans)
    assert opts['nrows']*opts['ncols'] >= nchans

    # Velocity axis
    args.cube = args.cube.with_spectral_unit(u.km/u.s,
            velocity_convention='radio')
    vel = args.cube.spectral_axis
    ind = np.argsort(vel[ichans])
    ichans = list(np.array(ichans)[ind])
    ind = len(ichans)//2
    vel = vel - vel[ichans][ind]

    # Setup tile plotter
    cubewcs = args.cube.wcs.sub(('longitude','latitude'))
    fig = MapsPlotter(config=args.config[0], section=args.section[0],
            projection=cubewcs, **opts)

    # Get vmin and vmax
    if len(args.images) == 1:
        data = np.squeeze(args.images[0].data)
        unit = u.Unit(args.images[0].header['BUNIT'])
        imagewcs = WCS(args.images[0].header).sub(('longitude','latitude'))
    else:
        data = args.cube.unmasked_data[ichans,:,:].value
    vmin, vmax = auto_vminmax(data)
    vmin = fig.config.get('vmin', fallback=vmin)
    vmax = fig.config.get('vmax', fallback=vmax)
    a = fig.config.get('a', fallback=100)

    # Levels
    levels = auto_levels(args.cube.unmasked_data[ichans,:,:].value,
            n=fig.config.getint('nlevels', fallback=10), 
            stretch=fig.config.get('stretch', fallback='log'))
    print(levels)

    # Center position
    if 'center' in fig.config and 'radius' in fig.config:
        cen = SkyCoord(fig.config['center'], unit=(u.hourangle, u.deg),
                frame='fk5')
        radius = fig.config['radius'].split()
        radius = float(radius[0]) * u.Unit(radius[1])
    else:
        cen = radius = None

    # Plot
    for ax,i in zip(fig.axes, ichans):
        cbar = ichans.index(i)==0
        ax = fig.get_mapper(ax, include_cbar=cbar, vmin=vmin, vmax=vmax, a=a)

        if len(args.images) == 0:
            data = args.cube.unmasked_data[i,:,:].value
            unit = args.cube.unmasked_data[i,:,:].unit
            wcs = cubewcs
            contours = contours_wcs = None
        else:
            contours = args.cube.unmasked_data[i,:,:].value
            contours_wcs = cubewcs
            wcs = imagewcs

        lev_color = 'w'

        ax.plot_map(data, wcs=wcs, r=radius, position=cen, contours=contours,
                contours_wcs=contours_wcs, levels=levels, colors=lev_color)
        
        # Color bar
        if cbar:
            orientation = 'horizontal' if fig.config.getboolean('hcbar') else \
                    'vertical'
            ax.plot_cbar(fig.fig, 
                    label='Intensity (%s)' % unit.to_string('latex_inline'),
                    orientation=orientation)

        # Velocity labels
        label = '%.2f %s' % (vel[i].value, vel.unit.to_string('latex_inline'))
        ax.annotate(label, xy=(0.1,0.9), xytext=(0.1,0.9),
                xycoords='axes fraction', zorder=3, color='k', 
                backgroundcolor='w')

    fig.auto_config()

    return fig

