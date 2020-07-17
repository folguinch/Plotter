import os

import astropy.units as u
from astropy.wcs import WCS
import numpy as np
#try:
from myutils.array_utils import load_mixed_struct_array as data_loader
#except ImportError:
#    data_loader = np.loadtxt

from utils_image import moment as get_moment
from src.utils import auto_vminmax

def load_cube(args):
    if args.cubename:
        args.cube = get_cube(args.cubename[0])
    else:
        args.cube = None

    return args

def load_fits(args):
    if args.imagenames:
        for img in args.imagenames:
            args.images += [get_fits(img)]
    else:
        args.images = None

    return args

def load_overplot(args):
    if args.overplot:
        if args.oplabel:
            if len(args.oplabel)==1:
                labels = args.oplabel*len(args.overplot)
            else:
                assert len(args.oplabel)==len(args.overplot)
                labels = args.oplabel
        else:
            labels = [None]*len(args.overplot)

        if args.opcolor:
            if len(args.opcolor)==1:
                colors = args.opcolor*len(args.overplot)
            else:
                assert len(args.opcolor)==len(args.overplot)
                colors = args.opcolor
        else:
            colors = ['w']*len(args.overplot)
                
        over = []
        for l,o,c in zip(labels, args.overplot, colors):
            ext = os.path.splitext(o)[1]
            label = None if l=='-' else l
            d = {'label': label, 'colors': c}
            if ext=='.fits':
                img = get_fits(o)
                over += [(img, d)]
            elif ext=='.dat':
                over += [(get_dat(o), d)]
            else:
                raise NotImplementedError
        args.overplots = over
    return args

def load_data(args):
    for filename in args.filenames:
        args.logger.info('Loading: %s', filename)
        args.data += [data_loader(os.path.expanduser(filename))]
    return args

def get_fits(filename):
    from astropy.io import fits
    return fits.open(os.path.expanduser(filename))[0]

def get_cube(filename):
    from spectral_cube import SpectralCube
    return SpectralCube.read(os.path.expanduser(filename))

def get_dat(filename):
    return np.loadtxt(os.path.expanduser(filename))

def get_spec_at_pos(cfg, cube, filename=None):
    # Options
    xunit = u.Unit(cfg.get('xunit', fallback=1))
    if 'vlsr' in cfg:
        vlsr = cfg.getquantity('vlsr')
    else:
        vlsr = None
    if 'restfreq' in cfg:
        restfreq = cfg.getquantity('restfreq', fallback=None)
    else:
        restfreq = None

    # Spectral axis
    if xunit.is_equivalent(u.km/u.s):
        cube0 = cube.with_spectral_unit(xunit, velocity_convention='radio',
                rest_value=restfreq)
    elif xunit.is_equivalent(u.Hz) and vlsr is not None and \
            restfreq is not None:
        # vlsr to freq
        freq_to_vel = u.doppler_radio(restfreq)
        flsr = vlsr.to(xunit, equivalencies=freq_to_vel)

        # Convert to velocity
        cube0 = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                rest_value=flsr)

        # Convert back
        cube0 = cube0.with_spectral_unit(xunit, velocity_convention='radio',
                rest_value=restfreq)
    elif xunit.is_equivalent(u.Hz):
        # Convert to requested unit
        cube0 = cube.with_spectral_unit(xunit, velocity_convention='radio',
                rest_value=restfreq)
    else:
        cube0 = cube

    # Spectra position
    coord = cfg.getskycoord('position')
    wcs = cube.wcs.sub(['longitude', 'latitude'])
    x, y = wcs.all_world2pix([[coord.ra.degree,coord.dec.degree]], 0)[0]
    x, y = int(x), int(y)

    # Spectral range
    if 'freqrange' in cfg or 'velrange' in cfg:
        if 'freqrange' in cfg:
            xlim = cfg.getquantity('freqrange')
        elif 'velrange' in cfg:
            xlim = cfg.getquantity('velrange')
        aux = cube0.spectral_slab(xlim[0], xlim[1])
        xaxis = aux.spectral_axis
        spec = aux[:, y, x]

        # Shift velocity
        if xaxis.unit.is_equivalent(u.km/u.s) and vlsr is not None:
            xaxis = xaxis - vlsr
    else:
        spec = cube0[:, y, x]
        xaxis = cube0.spectral_axis

        # Shift velocity
        if xaxis.unit.is_equivalent(u.km/u.s) and vlsr is not None:
            xaxis = xaxis - vlsr

    if filename:
        with open(filename, 'w') as out:
            out.write('#v\tF\n')
            out.write('#{0.unit}\t{1.unit}\n'.format(xaxis, spec))
            for dt in zip(xaxis, spec[:]):
                out.write('{0.value:f}\t{1.value:f}\n'.format(*dt))

    return xaxis, spec

# For multi plot
def multi_map(cfg):
    data = get_fits(cfg['filename'])
    projection = WCS(data.header).sub(['longitude','latitude'])
    
    return [[data, projection]], projection

def multi_spectrum(cfg):
    data = get_cube(cfg['filename'])
    data = get_spec_at_pos(cfg, data, filename=cfg.get('outspec',
        fallback=None))
    projection = 'rectilinear'

    return [data], projection

def multi_moment_map(cfg):
    mom = cfg.getint('moment')
    if 'filename' in cfg and os.path.isfile(cfg['filename']):
        data, projection = multi_map(cfg)
        data = data[0][0]
    else:
        cube = get_cube(cfg['cubename'])
        data = get_moment(cube, mom, cfg, filename=cfg.get('filename',
            fallback=None))
        projection = WCS(data.header).sub(['longitude','latitude'])
    if mom==1 and 'vlsr' in cfg:
        bunit = u.Unit(data.header['BUNIT'])
        vlsr = cfg.getquantity('vlsr').to(bunit).value
        data.data = data.data - vlsr
    elif mom==2 and cfg.get('stretch', fallback='linear').lower()=='log':
        mask = data.data<=0.
        data.data[mask] = float('nan')

    return [[data, projection]], projection

