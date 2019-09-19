import os

def load_cube(args):
    from spectral_cube import SpectralCube
    if args.cubename:
        args.cube = SpectralCube.read(os.path.expanduser(args.cubename[0]))
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

def get_fits(filename):
    from astropy.io import fits
    return fits.open(os.path.expanduser(filename))[0]
