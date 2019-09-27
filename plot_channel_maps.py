from builtins import range, zip

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from src.map_plotter import MapsPlotter
from src.utils import auto_vminmax, auto_levels
from utils import all_mapfig_setup, plot_single_map, get_shape

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
    cen, radius, markers, orientation = all_mapfig_setup(fig)

    # Get global vmin and vmax
    if len(args.images) == 1:
        img = args.images[0]
        data = np.squeeze(img.data)
    else:
        data = args.cube.unmasked_data[ichans,:,:].value
    vmin, vmax = auto_vminmax(data)
    vmin = fig.config.getfloat('vmin', fallback=vmin)
    vmax = fig.config.getfloat('vmax', fallback=vmax)
    a = fig.config.getfloat('a', fallback=100)

    # Levels
    levels = auto_levels(args.cube.unmasked_data[ichans,:,:].value,
            n=fig.config.getint('nlevels', fallback=10), 
            stretch=fig.config.get('stretch', fallback='log'))
    print(levels)

    # Plot
    for ax,i in zip(fig.axes, ichans):
        cbar = ichans.index(i)==0
        #ax = fig.get_mapper(ax, include_cbar=cbar, vmin=vmin, vmax=vmax, a=a)
        bmaj = args.cube.beams[i].major.to(u.deg).value
        bmin = args.cube.beams[i].minor.to(u.deg).value
        bpa = args.cube.beams[i].pa.to(u.deg).value

        if len(args.images) == 0:
            img = fits.PrimaryHDU(args.cube.unmasked_data[i,:,:].value,
                    header=cubewcs.to_header())
            img.header['BMAJ'] = bmaj
            img.header['BMIN'] = bmin
            img.header['BPA'] = bpa
            img.header['BUNIT'] = args.cube.header['BUNIT']
            contours = None
        else:
            contours = fits.PrimaryHDU(args.cube.unmasked_data[i,:,:].value,
                    header=cubewcs.to_header())
            contours.header['BMAJ'] = bmaj
            contours.header['BMIN'] = bmin
            contours.header['BPA'] = bpa

        label = '%.2f %s' % (vel[i].value, vel.unit.to_string('latex_inline'))

        ax = plot_single_map(ax, fig, img, contours=contours, cen=cen, 
                radius=radius, levels=levels, 
                cbar_orientation=orientation if cbar else None,
                markers=markers, axlabel=label, vmin=vmin, vmax=vmax, a=a,
                include_cbar=cbar)

    fig.auto_config()

    return fig

