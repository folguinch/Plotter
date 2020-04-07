import argparse
import math
import os
from configparser import ConfigParser

import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import myutils.image_utils as myimg_utils

from src.utils import auto_vminmax, auto_levels
from src.logger import get_logger

def postprocess(fig, args):
    # Save figure
    plotname = os.path.expanduser(args.plotname[0])
    fig.savefig(plotname)
    if args.png:
        plotname = os.path.splitext(plotname)[0] + '.png'
        fig.savefig(plotname)

def read_config(filename, parser=None):
    # Check that file exists
    if not os.path.isfile(filename):
        raise IOError('%s does not exist' % filename)

    # Create parser
    if parser is None:
        parser = ConfigParser()

    # Read parser
    parser.read(filename)

    return parser

def get_quantity(x):
    return _get_quantity(*tuple(x.split()))

def _get_quantity(x,y):
    return float(x)*u.Unit(y)

def get_shape(args, total, default_cols=3, minimize=False):
    if args.shape:
        rows, cols = args.shape
    elif args.rows:
        rows = args.rows[0]
        cols = math.ceil(float(total)/rows)
    elif args.cols:
        rows = math.ceil(float(total)/args.cols[0])
        cols = args.cols[0]
    else:
        cols = default_cols
        rows = math.ceil(float(total)/cols)

    if minimize and int(cols*rows)>total:
        newarg = argparse.ArgumentParser()
        newarg.set_defaults(shape=None, rows=None)
        newarg.add_argument('--cols', nargs=1, type=int, default=[default_cols])
        rows, cols = get_shape(newarg.parse_args(['--cols', '%i' % cols]), total)
    
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
        cens = fig.config['center'].split(',')
        cen = []
        for c in cens:
            cen += [SkyCoord(c, unit=(u.hourangle, u.deg), frame='fk5')]
        radii = fig.config['radius'].split(',')
        radius = []
        for r in radii:
            rr = r.split()
            radius += [float(rr[0]) * u.Unit(rr[1])]
    else:
        cen = radius = [None]

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
        if marker.get('label')=='-':
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
            assert len(args.opmapping)==len(args.overplots)
            over = [args.overplots[i] \
                    for i, j in enumerate(args.opmapping) if j==n]
            #over = [over, args.overplots[1], args.overplots[2]]
        else:
            try:
                over = args.overplots[0][n]
                over = [[over], args.overplots[1], args.overplots[2]]
            except IndexError:
                over = None
    else:
        over = None
    return over

def get_levels(data, loc, fig, levels=None, self_contours=True, 
        dtype='intensity', logger=None, base=2):
    if not self_contours or dtype=='velocity':
        levels = None
    elif levels is not None:
        pass
    elif 'rms' in fig.config:
        # Levels as multiples of the rms
        rms = float(fig.get_value('rms', ax=loc))
        nsigma = int(fig.get_value('nsigma', default=5, ax=loc))
    else:
        levels = fig.get_value('levels', levels, loc, sep=',')
        if levels is not None:
            try:
                levels = map(float, levels.split())
            except:
                if logger is not None:
                    logger.warn('Could not determine levels')
                self_contours = False
        else:
            pass
    return levels, self_contours

def get_extent(img, loc, fig, logger):
    # Get coordinate axes
    xaxis, yaxis = myimg_utils.get_coord_axes(img)
    xaxis = xaxis.to(u.arcsec).value
    yaxis = yaxis.to(u.km/u.s).value

    # Coodinate offset
    try:
        xoffset = float(fig.get_value('xoffset', default=None, ax=loc))
        logger.info('Shifting x axis by: %f', xoffset)
        xaxis = xaxis - offset
    except TypeError:
        xaxis = xaxis - xaxis[len(xaxis)//2]
    try:
        yoffset = float(fig.get_value('yoffset', default=None, ax=loc))
        logger.info('Shifting x axis by: %f', yoffset)
        yaxis = yaxis - offset
    except TypeError:
        yaxis = yaxis - yaxis[len(yaxis)//2]

    # Extent
    extent = [xaxis[0], xaxis[-1], yaxis[0], yaxis[-1]]
    return extent

def all_mapfig_setup(fig):
    # Center position
    cen, radius = get_center(fig)

    # Color bar
    orientation = get_cbar_orientation(fig)

    # Markers
    markers = get_markers(fig)

    return cen, radius, markers, orientation 

def get_axis_label(args, n, detail):
    label = '(' + chr(ord('a')+n) + ')'

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
        dtype='intensity', levels=None, self_contours=True, nsigmalevel=None,
        global_levels=False, cbar_orientation=None, markers=None, 
        skip_marker_label=False, axlabel=None, **kwargs):
    # Change default self_contours
    if dtype=='velocity':
        self_contours = False
    if levels is not None:
        global_levels = True
    
    # Data
    data = np.squeeze(img.data)
    try:
        unit = u.Unit(img.header['BUNIT'])
    except KeyError:
        logger.info('Unit not found')
        unit = ''
    if dtype == 'pvmap':
        wcs = None
        extent = get_extent(img, loc, fig, logger)
    else:
        wcs = WCS(img.header).sub(('longitude','latitude'))
        extent = None

    # Colorbar label
    realdtype = dtype if dtype!='pvmap' else 'intensity'
    try:
        cbarlabel = '%s (%s)' % (realdtype.capitalize(), unit.to_string('latex_inline'))
    except AttributeError:
        cbarlabel = '%s' % (realdtype.capitalize(),)
    cbarlabel = fig.get_value('cbarlabel', cbarlabel, loc)

    # rms, nsigma
    try:
        rms = float(fig.get_value('rms', ax=loc))
    except TypeError:
        rms = None
    nsigma = float(fig.get_value('nsigma', default=5., ax=loc))
    if nsigmalevel is None:
        try:
            nsigmalevel = int(fig.get_value('nsigmalevel', default=None,
                ax=loc))
        except TypeError:
            pass

    # Get vmin and vmax
    if kwargs.get('vmin') is not None and kwargs.get('vmax') is not None:
        vmin = kwargs['vmin']
        vmax = kwargs['vmax']
    else:
        vmin, vmax = auto_vminmax(data, dtype=dtype, rms=rms)
        vmin = float(fig.get_value('vmin', vmin, loc))
        vmax = float(fig.get_value('vmax', vmax, loc))
    logger.info('Plotting data with vmin, vmax: %.3e, %.3e', vmin, vmax)

    # Create axis, auto determine if cbar is needed
    logger.info('Creating axis with%s color bar',
            ['out',''][int(kwargs.get('include_cbar', False))])
    ax = fig.get_mapper(loc, vmin=vmin, vmax=vmax, a=kwargs.get('a'), 
            projection=wcs, include_cbar=kwargs.get('include_cbar'))

    # Levels
    if self_contours and global_levels and levels is None:
        levels = auto_levels(data, rms=rms, nsigma=nsigma)
    elif not self_contours and global_levels and levels is None:
        nlevels = int(fig.get_value('nlevels', default=None, ax=loc))
        if rms is not None and nlevels is not None:
            levels = auto_levels(rms=rms, nsigma=nsigma, nlevels=nlevels)
        else:
            raise Exception('Could not determine global levels for contours')
    lev_color = kwargs.get('contour_color', 
            fig.config.get('contour_color', fallback='w'))

    # Contour line width
    try:
        linewidths = map(int, fig.get_value('linewidths', ax=loc,
            sep=',').split())
        if len(linewidths)==1:
            linewidths = linewidths[0]
    except (TypeError,AttributeError):
        linewidths = None

    # Plot data
    ax.plot_map(data, wcs=wcs, r=radius, position=cen, 
            self_contours=self_contours, levels=levels, rms=rms, nsigma=nsigma,
            nsigmalevel=nsigmalevel, colors=lev_color, label=cbarlabel, 
            extent=extent, linewidths=linewidths)

    # Plot contours
    if contours is not None:
        if not self_contours and 'levels' in fig.config:
            levels = map(float, 
                    fig.get_value('levels', ax=loc, sep=',').split())
        else:
            levels = None

        try:
            contour_rms = float(fig.get_value('rms_contour', default=None, 
                ax=loc))
        except TypeError:
            logger.warn('rms_contour not in config file')
            contour_rms = None
        
        # Plot
        logger.info('Overplot contours: %s', levels)
        plot_contours(ax, contours, levels=levels, nsigma=nsigma,
                rms=contour_rms, linewidths=linewidths)

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

def plot_contours(ax, contours, levels=None, nsigma=5., zorder=2, rms=None,
        linewidths=None):
    for j,(img,opt) in enumerate(contours):
        data = np.squeeze(img.data)
        wcs = WCS(img.header).sub(('longitude','latitude'))
        ax.plot_contours(data, levels=levels, wcs=wcs, rms=rms, nsigma=nsigma, 
                zorder=zorder+j, linewidths=linewidths, **opt)

def splot(loc, fig, data, logger, overplots=[], cols=[0,1], errorcol=None):
    # Axis
    ax = fig.get_axis(loc)

    # Iterate over data
    xunit = None
    yunit = None
    if overplots is None:
        overplots = []
    for i, d in enumerate([data]+overplots):
        # Error column
        if errorcol is not None:
            try:
                try:
                    err = d[:,errcol]
                except TypeError:
                    err = d[0][errcol]
                    if yunit is None:
                        yunit = d[1][errcol]
                    else:
                        err = err.to(yunit)
                fmt = fig.get_value('pfmt', n=i)
            except:
                logger.warn('No error column')
                err = None
                fmt = fig.get_value('lfmt', n=i)
        else:
            err = None
            fmt = fig.get_value('lfmt', n=i)

        # Data
        try:
            # Normal array
            x = d[:,cols[0]]
            y = d[:,cols[1]]
        except TypeError:
            # Array with units
            x = d[0][cols[0]]
            if xunit is None:
                xunit = d[1][cols[0]]
            x = x * d[1][cols[0]]
            x = x.to(xunit).value

            y = d[0][cols[1]]
            if yunit is None:
                yunit = d[1][cols[1]]
            y = y * d[1][cols[1]]
            y = y.to(yunit).value
        
        # plot
        if err is None:
            ax.plot(x, y, fmt)
        else:
            raise NotImplementedError
            #ax.errorbar(d[cols[0]], d[cols[1]])
    return ax, xunit, yunit

