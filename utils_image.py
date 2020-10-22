import astropy.units as u
import numpy as np

import myutils.astro_tools.cube_utils as cubeutils

import logger
import utils
import src.maths

LOG = logger.get_logger(__name__)

def moment(cube, mom, linecfg, vlsr=None, filename=None, filenamebase=None):
    # Get moment
    mmnt = cubeutils.moment_from_config(cube, mom, linecfg, vlsr=vlsr,
            filenamebase=filenamebase, filename=filename)

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

