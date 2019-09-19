import os, math

import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from src.utils import auto_vminmax, auto_levels

def postprocess(fig, args):
    # Save figure
    plotname = os.path.expanduser(args.plotname[0])
    fig.savefig(plotname)
    if args.png:
        plotname = os.path.splitext(plotname)[1] + '.png'
        fig.savefig(plotname)

def get_shape(args, total, default_cols=3):
    if args.shape:
        rows, cols = args.shape
    elif args.rows:
        rows = args.rows
        cols = math.ceil(float(total)/rows)
    elif args.cols:
        rows = math.ceil(float(total)/args.cols)
        cols = args.cols
    else:
        cols = default_cols
        rows = math.ceil(float(total)/cols)
    
    return int(rows), int(cols)

def get_cbar_orientation(fig):
    if fig.config.getboolean('vcbar'):
        orientation = 'vertical'
    elif fig.config.getboolean('hcbar'):
        orientation = 'horizontal'
    else:
        orientation = None
        print('WARNING: Map plot does not have colorbar')

    return orientation

def get_center(fig):
    if 'center' in fig.config and 'radius' in fig.config:
        cen = SkyCoord(fig.config['center'], unit=(u.hourangle, u.deg),
                frame='fk5')
        radius = fig.config['radius'].split()
        radius = float(radius[0]) * u.Unit(radius[1])
    else:
        cen = radius = None

    return cen, radius

def get_markers(config):
    markers = []
    for i, m in enumerate(config.get('markers', fallback='').split(',')):
        if m == '':
            break
        marker = {'loc':SkyCoord(m, unit=(u.hourangle, u.deg),frame='fk5')}
        try:
            marker['style'] = config.get('mstyles', 
                    fallback='x').replace(',',' ').split()[i]
        except IndexError:
            marker['style'] = 'x'
        try:
            marker['color'] = config.get('mcolors', 
                    fallback='r').replace(',',' ').split()[i]
        except IndexError:
            marker['color'] = 'r'
        try:
            marker['label'] = config.get('mlabels',
                    fallback=None).split(',')[i]
            try:
                loclabel = config.get('mlabloc', 
                        fallback=None).split(',')
                marker['labloc'] = SkyCoord(loclabel[i], unit=(u.hourangle,
                    u.deg),frame='fk5')
            except (AttributeError, IndexError):
                marker['labloc'] = marker['loc']
        except AttributeError:
            pass
        except IndexError:
            marker['label'] = ''
        markers += [marker]

    return markers

def all_mapfig_setup(fig):
    # Center position
    cen, radius = get_center(fig)

    # Color bar
    orientation = get_cbar_orientation(fig)

    # Markers
    markers = get_markers(fig.config)

    return cen, radius, markers, orientation 

def get_axis_label(args, n, detail):
    if args.axlabel:
        label = '(' + chr(ord('a')+n) + ')'
    else:
        label = ''

    if args.detailaxlabel:
        label = ('%s %s' % (label, detail)).strip()
    elif args.axlabels:
        try:
            label = ('%s %s' % (label, args.axlabels[n])).strip()
        except IndexError:
            label = None
    elif args.axlabel:
        pass
    else:
        label = None

    return label

def plot_single_map(ax, fig, img, cen=None, radius=None, 
        dtype='intensity', self_levels=True, cbar_orientation=None, 
        markers=None, skip_marker_label=False, axlabel=None):
    # Data
    data = np.squeeze(img.data)
    unit = u.Unit(img.header['BUNIT'])
    wcs = WCS(img.header).sub(('longitude','latitude'))

    # Get vmin and vmax
    vmin, vmax = auto_vminmax(data, dtype=dtype)
    vmin = float(fig.get_value('vmin', vmin, ax))
    vmax = float(fig.get_value('vmax', vmax, ax))
    a = fig.config.getfloat('a', fallback=100)

    # Levels
    if not self_levels or dtype == 'velocity':
        levels = None
    else:
        levels = fig.get_value('levels', None, ax, sep=',')
        if levels is not None:
            levels = map(float, levels.split())
        else:
            levels = auto_levels(data,
                    n=fig.config.getint('nlevels', fallback=10), 
                    stretch=fig.config.get('stretch', fallback='log'))
    print(levels)
    lev_color = fig.config.get('level_color', fallback='w')

    # Create axis, auto determine if cbar is needed
    ax = fig.get_mapper(ax, vmin=vmin, vmax=vmax, a=a, projection=wcs)

    # Plot data
    ax.plot_map(data, wcs=wcs, r=radius, position=cen, contours=data,
            contours_wcs=wcs, levels=levels, colors=lev_color)

    # Plot markers
    if markers is not None:
        ax.plot_markers(markers, skip_label=skip_marker_label)

    # Color bar
    if cbar_orientation is not None:
        ax.plot_cbar(fig.fig, 
                label='Intensity (%s)' % unit.to_string('latex_inline'),
                orientation=cbar_orientation)

    # Axis label
    if axlabel:
        ax.annotate(axlabel, xy=(0.1,0.9), xytext=(0.1,0.9), 
            xycoords='axes fraction', color='k', backgroundcolor='w')

    # Beam
    if 'BMAJ' in img.header and 'BMIN' in img.header and 'BPA' in img.header:
        ax.plot_beam(img.header)

    return ax
