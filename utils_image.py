import astropy.units as u
import numpy as np

import logger
import utils

import src.maths

LOG = logger.get_logger(__name__)

def get_restfreq(cube):
    try:
        restfreq = cube.header['RESTFRQ'] * u.Hz
    except KeyError:
        restfreq = cube.header['RESTFREQ'] * u.Hz
    return restfreq

def get_spectral_limits(cube, cfg, vlsr=None):
    # Turn everything into frequency or velocity
    spec_unit = cube.spectral_axis.unit
    rangefmt = '{0.value[0]:.4f} {0.value[1]:.4f} {0.unit}'
    if 'chan_width' in cfg and 'freq' in cfg and vlsr is not None:
        # Get spectral axis
        spaxis = cube.spectral_axis
        restfreq = get_restfreq(cube)
        linefreq = utils.get_quantity(cfg['freq'])

        # Convert to velocity from line
        vel_to_freq = u.doppler_radio(restfreq)
        spaxis = spaxis.to(linefreq.unit, equivalencies=vel_to_freq)
        freq_to_vel = u.doppler_radio(linefreq)
        spaxis = spaxis.to(vlsr.unit, equivalencies=freq_to_vel)
        spaxis = spaxis - vlsr

        # Closest value to 0.0
        ind = np.nanargmin(np.abs(spaxis.value))
        chanwidth = int(cfg['chan_width'])
        if ind-chanwidth < 0:
            chmin = 0
        else:
            chmin = ind - chanwidth
        if ind+chanwidth >= len(spaxis):
            chmax = len(spaxis)-1
        else:
            chmax = ind + chanwidth
        LOG.info('Using channel range = %i %i', chmin, chmax)
        return chmin, chmax
    elif 'freq_range' in cfg:
        freq_range = utils.get_quantity(cfg['freq_range'])
        LOG.info('Using frequency range = %s', 
                rangefmt.fmt(freq_range.to(u.GHz)))
        if spec_unit.is_equivalent(freq_range.unit):
            return freq_range[0], freq_range[1]
        else:
            # Convert to spectral units
            restfreq = get_restfreq(cube)
            freq_to_vel = u.doppler_radio(restfreq)
            lim = freq_range.to(spec_unit, equivalencies=freq_to_vel)
            return np.min(lim), np.max(lim)
    elif 'vel_range' in cfg:
        vel_range = utils.get_quantity(cfg['vel_range'])
        LOG.info('Using velocity range = %s', 
                rangefmt.fmt(vel_range.to(u.km/u.s)))
        if spec_unit.is_equivalent(vel_range.unit):
            return args.vel_range[0], args.vel_range[1]
        else:
            # Convert to spectral units
            restfreq = get_restfreq(cube)
            freq_to_vel = u.doppler_radio(restfreq)
            lim = args.vel_range.to(spec_unit, equivalencies=freq_to_vel)
            return np.min(lim), np.max(lim)
    elif 'chan_range' in cfg:
        chmin, chmax = sorted(map(int, cfg['chan_range'].split()))
        LOG.info('Channel range = %i %i', chmin, chmax)
        # Get middle channel
        if 'midchan' in cfg:
            midchan = int(linecfg['midchan'])
            LOG.info('Using channel %i as line center', midchan)
            mindiff = min(abs(midchan-chmin), abs(midchan-chmax))
            chmin = midchan - mindiff
            chmax = midchan + mindiff
            LOG.info('Re-centred channel range = %i, %i', chmin, chmax)
        return chmin, chmax
    else:
        LOG.warn('No spectral limit')
        return None, None

def get_subcube(cube, cfg, vlsr=None):
    # Extract the spectral slab
    low, high = get_spectral_limits(cube, cfg, vlsr=vlsr)
    if hasattr(low, 'unit'):
        subcube = cube.spectral_slab(low, high)
    elif low is None:
        subcube = cube
    else:
        subcube = cube[low:high+1,:,:]

    if 'RMS' in cube.header and 'RMS' not in subcube.header:
        subcube.header['RMS'] = cube.header['RMS']

    return subcube

def moment(cube, mom, linecfg, vlsr=None, filename=None):
    # Info from linecfg
    #try:
    #    restfreq = utils.get_quantity(linecfg['freq'])
    #except:
    #    restfreq = utils.get_quantity(linecfg['restfreq'])
    #chmin, chmax = map(int, linecfg['chanran'].split())
    #LOG.info('Rest frequency = %s', restfreq)
    #LOG.info('Channel range = %i, %i', chmin, chmax)

    ## Get middle channel
    #try:
    #    midchan = int(linecfg['midchan'])
    #    LOG.info('Using channel %i as line center')
    #    mindiff = min(abs(midchan-chmin), abs(midchan-chmax))
    #    chmin = midchan - mindiff
    #    chmax = midchan + mindiff
    #    LOG.info('Re-centred channel range = %i, %i', chmin, chmax)
    #except:
    #    pass

    ## Convert velocity
    #cube = cube.with_spectral_unit(u.km/u.s,
    #        velocity_convention='radio', rest_value=restfreq)
    #vel = cube.spectral_axis
    #LOG.info('Velocity range = %s, %s',vel[chmin], vel[chmax])
    #subcube = cube[chmin:chmax+1,:,:]

    # Get subcube
    LOG.info('Extracting sub-cube:')
    subcube = get_subcube(cube, linecfg, vlsr=vlsr)

    # Convert to velocity
    try:
        restfreq = utils.get_quantity(linecfg['freq'])
    except:
        restfreq = utils.get_quantity(linecfg['restfreq'])
    LOG.info('Line rest frequency = %s', restfreq)
    subcube = subcube.with_spectral_unit(u.km/u.s,
            velocity_convention='radio', rest_value=restfreq)
    vel = subcube.spectral_axis
    LOG.info('Velocity range = %s, %s',vel[0], vel[-1])

    # Moment
    if mom>0:
        if 'low' in linecfg:
            low = utils.get_quantity(linecfg['low'])
            LOG.info('Using lower flux limit: %s', 
                    '{0.value:.3f} {0.unit}'.format(low))
            mask = subcube >= low
        elif 'nsigma' in linecfg:
            nsigma = float(linecfg['nsigma'])
            if 'rms' in linecfg:
                rms = utils.get_quantity(linecfg['rms'])
                LOG.info('Using config rms: %s',
                        '{0.value:.3e} {0.unit}'.format(rms))
            elif 'RMS' in cube.header:
                rms = cube.header['RMS'] * cube.unit
                LOG.info('Using header rms: %s',
                        '{0.value:.3e} {0.unit}'.format(rms))
            else:
                try:
                    rms = src.maths.quick_rms(cube.unmasked_data)
                except TypeError:
                    rms = src.maths.quick_rms(cube.unmasked_data[:])
                LOG.info('Using cube rms: %s',
                        '{0.value:.3e} {0.unit}'.format(rms))
            low = nsigma*rms
            LOG.info('%isigma lower flux limit: %s', nsigma, 
                    '{0.value:.3f} {0.unit}'.format(low))
            mask = subcube >= low
        else:
            mask = subcube.mask
        if 'up' in linecfg:
            up = utils.get_quantity(linecfg['up'])
            mask = mask & (subcube<=up)
        subcube = subcube.with_mask(mask)
    if mom==2:
        mmnt = subcube.linewidth_fwhm()
    else:
        mmnt = subcube.moment(order=mom)

    if filename:
        LOG.info('Saving moment: %s', filename)
        mmnt.write(filename, overwrite=True)

    return mmnt.hdu

def line_center(cube, linecfg, nsigma=5, rms=None):
    # Line limits
    restfreq = utils.get_quantity(linecfg['freq'])
    chmin, chmax = map(int, linecfg['chanran'].split())
    LOG.info('Channel range = %i, %i', chmin, chmax)

    # rms
    if rms is None:
        rms = maths.quick_rms(cube.unmasked_data.value)

    # Get spectral axis
    cube = cube.with_spectral_unit(u.km/u.s,
            velocity_convention='radio', rest_value=restfreq)
    vel = cube.spectral_axis
    subcube = cube[chmin:chmax+1,:,:] 
    
    # Mask
    mask = subcube > nsigma*rms*cube.unit
    subcube = subcube.with_mask(mask)
    mask = np.any(mask, axis=0)

    # Iterate over spectra
    I, J = np.indices(mask.shape)
    #for i, j in 

