import os, math

import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from src.utils import auto_vminmax, auto_levels
from src.logger import get_logger

def postprocess(fig, args):
    # Save figure
    plotname = os.path.expanduser(args.plotname[0])
    fig.savefig(plotname)
    if args.png:
        plotname = os.path.splitext(plotname)[0] + '.png'
        fig.savefig(plotname)

def get_quantity(x):
    return _get_quantity(*tuple(x.split()))

def _get_quantity(x,y):
    return float(x)*u.Unit(y)

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

def get_markers(fig):
    config = fig.config
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
        marker['color'] = fig.get_value('mcolors', 'r', n=i, sep=',')
        marker['size'] = float(fig.get_value('msizes', 50, n=i))
        try:
            marker['label'] = config.get('mlabels',
                    fallback=None).split(',')[i].strip()
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
        if marker['label']=='-':
            marker['label'] = ''
        try:
            marker['legend'] = config.get('mlegend', 
                    fallback=None).split(',')[i].strip()
        except AttributeError:
            marker['legend'] = None
        except IndexError:
            marker['legend'] = None
        if marker['legend']=='-' or marker['legend']=='':
            marker['legend'] = None
        markers += [marker]

    return markers

def get_overplots(args, n):
    if args.overplots:
        if args.overall:
            over = args.overplots
        elif n in args.opmapping:
            assert len(args.opmapping)==len(args.overplots[0])
            over = [args.overplots[0][i] \
                    for i, j in enumerate(args.opmapping) if j==n]
            over = [over, args.overplots[1], args.overplots[2]]
        else:
            try:
                over = args.overplots[0][n]
                over = [[over], args.overplots[1], args.overplots[2]]
            except IndexError:
                over = None
    else:
        over = None
    return over

def all_mapfig_setup(fig):
    # Center position
    cen, radius = get_center(fig)

    # Color bar
    orientation = get_cbar_orientation(fig)

    # Markers
    markers = get_markers(fig)

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

def plot_single_map(loc, fig, img, logger, contours=None, cen=None, radius=None, 
        dtype='intensity', levels=None, self_contours=True, cbar_orientation=None,
        markers=None, skip_marker_label=False, axlabel=None,
        **kwargs):
    # Data
    data = np.squeeze(img.data)
    unit = u.Unit(img.header['BUNIT'])
    wcs = WCS(img.header).sub(('longitude','latitude'))
    cbarlabel = '%s (%s)' % (dtype.capitalize(), unit.to_string('latex_inline'))

    # Get vmin and vmax
    if kwargs.get('vmin') is not None and kwargs.get('vmax') is not None:
        vmin = kwargs['vmin']
        vmax = kwargs['vmax']
    else:
        vmin, vmax = auto_vminmax(data, dtype=dtype)
        vmin = float(fig.get_value('vmin', vmin, loc))
        vmax = float(fig.get_value('vmax', vmax, loc))
    logger.info('Plotting data with vmin, vmax: %.3e, %.3e', vmin, vmax)

    # Create axis, auto determine if cbar is needed
    logger.info('Creating axis with%s color bar',
            ['out',''][int(kwargs.get('include_cbar', False))])
    ax = fig.get_mapper(loc, vmin=vmin, vmax=vmax, a=kwargs.get('a'), 
            projection=wcs, include_cbar=kwargs.get('include_cbar'))

    # Levels
    nlevels = fig.config.getint('nlevels', fallback=10)
    if not self_contours or dtype=='velocity':
        levels = None
    elif levels is not None:
        pass
    else:
        levels = fig.get_value('levels', levels, loc, sep=',')
        if levels is not None:
            levels = map(float, levels.split())
        else:
            pass
    print(levels)

    # Plot data
    lev_color = kwargs.get('contour_color', 
            fig.config.get('contour_color', fallback='w'))
    self_levs = ax.plot_map(data, wcs=wcs, r=radius, position=cen,
            self_contours=self_contours and dtype!='velocity', 
            levels=levels, colors=lev_color,
            label=cbarlabel, nlevels=nlevels)
    if self_contours and levels is None:
        levels = self_levs
        print(levels)

    # Plot contours
    if contours is not None and dtype!='velocity':
        if not self_contours and levels is None:
            levels = auto_levels(None, vmin=contours[1], vmax=contours[2],
                    stretch=ax.stretch, n=nlevels)
        plot_contours(ax, contours[0], levels)

    # Plot markers
    if markers is not None:
        ax.plot_markers(markers, skip_label=skip_marker_label)

    # Color bar
    if fig.has_cbar(loc) and cbar_orientation is not None:
        ax.plot_cbar(fig.fig, orientation=cbar_orientation,
                labelpad=fig.config.getfloat('labelpad', fallback=10))

    # Axis label
    if axlabel:
        ax.annotate(axlabel, xy=(0.1,0.9), xytext=(0.1,0.9), 
            xycoords='axes fraction', color='k', backgroundcolor='w', zorder=3)

    # Beam
    if 'BMAJ' in img.header and 'BMIN' in img.header and 'BPA' in img.header:
        color = fig.get_value('beamcolor', 'w', loc)
        ax.plot_beam(img.header, color=color)

    # Legend
    #if legend:
    #    ax.legend(auto=True, loc=4, match_colors=True,
    #            fancybox=fig.config.getboolean('fancybox', fallback=False),
    #            framealpha=fig.config.getfloat('framealpha', fallback=None),
    #            facecolor=fig.config.get('facecolor', fallback=None))

    return ax

def plot_contours(ax, contours, levels, zorder=2):
    for j,(img,opt) in enumerate(contours):
        data = np.squeeze(img.data)
        wcs = WCS(img.header).sub(('longitude','latitude'))
        ax.plot_contours(data, levels, wcs=wcs, zorder=zorder+j, **opt)

