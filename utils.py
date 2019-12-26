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
        rows = args.rows[0]
        cols = math.ceil(float(total)/rows)
    elif args.cols:
        rows = math.ceil(float(total)/args.cols[0])
        cols = args.cols[0]
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
        dtype='intensity', levels=None, self_contours=True,
        global_levels=False, cbar_orientation=None, markers=None, 
        skip_marker_label=False, axlabel=None, **kwargs):
    # Change default self_contours
    if dtype=='velocity':
        self_contours = False
    if levels is not None:
        global_levels = True
    
    # Data
    data = np.squeeze(img.data)
    unit = u.Unit(img.header['BUNIT'])
    wcs = WCS(img.header).sub(('longitude','latitude'))
    cbarlabel = '%s (%s)' % (dtype.capitalize(), unit.to_string('latex_inline'))
    cbarlabel = fig.get_value('cbarlabel', cbarlabel, loc)

    # rms, nsigma
    try:
        rms = float(fig.get_value('rms', ax=loc))
    except TypeError:
        rms = None
    nsigma = float(fig.get_value('nsigma', default=5., ax=loc))

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
    #levels, self_contours = get_levels(data, loc, fig, levels=levels,
    #        self_contours=self_contours, dtype=dtype, logger=logger)
    #logger.info('Levels: %r', levels)

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

    # Plot data
    ax.plot_map(data, wcs=wcs, r=radius, position=cen, 
            self_contours=self_contours, levels=levels, rms=rms, nsigma=nsigma,
            colors=lev_color, label=cbarlabel)

    # Plot contours
    if contours is not None:
        if not self_contours and 'levels' in fig.config:
            levels = map(float, 
                    fig.get_value('levels', ax=loc, sep=',').split())
        try:
            contour_rms = float(fig.get_value('rms_contour', default=None, 
                ax=loc))
        except TypeError:
            logger.warn('rms_contour not in config file')
            contour_rms = None

        plot_contours(ax, contours, levels=levels, nsigma=nsigma,
                rms=contour_rms)

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

def plot_contours(ax, contours, levels=None, nsigma=5., zorder=2, rms=None):
    for j,(img,opt) in enumerate(contours):
        data = np.squeeze(img.data)
        wcs = WCS(img.header).sub(('longitude','latitude'))
        ax.plot_contours(data, levels=levels, wcs=wcs, rms=rms, nsigma=nsigma, 
                zorder=zorder+j, **opt)

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

