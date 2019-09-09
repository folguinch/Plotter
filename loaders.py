import os

def load_cube(args):
    from spectral_cube import SpectralCube
    if args.cubename:
        args.cube = SpectralCube.read(os.path.expanduser(args.cubename[0]))
    else:
        args.cube = None

    return args

def load_fits(args):
    from astropy.io import fits
    if args.imagename:
        args.image = fits.open(os.path.expanduser(args.imagename[0]))[0]
    else:
        args.image = None

    return args

