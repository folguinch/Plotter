import astropy.units as u

import logger
import utils

LOG = logger.get_logger(__name__)

def moment(cube, mom, linecfg, filename=None):
    # Info from linecfg
    restfreq = utils.get_quantity(linecfg['freq'])
    chmin, chmax = map(int, linecfg['chanran'].split())
    LOG.info('Extracting sub-cube:')
    LOG.info('Rest frequency = %s', restfreq)
    LOG.info('Channel range = %i, %i', chmin, chmax)

    # Convert velocity
    cube = cube.with_spectral_unit(u.km/u.s,
            velocity_convention='radio', rest_value=restfreq)
    vel = cube.spectral_axis
    LOG.info('Velocity range = %s, %s',vel[chmin], vel[chmax])
    subcube = cube[chmin:chmax+1,:,:]

    # Moment
    if mom>0:
        if 'low' in linecfg:
            low = utils.get_quantity(linecfg['low'])
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

