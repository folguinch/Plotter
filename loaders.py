import os

from src.utils import auto_vminmax

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
                
        rmin, rmax = float('inf'), -float('inf')
        over = []
        for l,o,c in zip(labels, args.overplot, colors):
            ext = os.path.splitext(o)[1]
            label = None if l=='-' else l
            d = {'label': label, 'colors': c}
            if ext=='.fits':
                img = get_fits(o)
                over += [(img, d)]
                vmin, vmax = get_img_ranges(img)
                rmin = min(rmin, vmin)
                rmax = max(rmin, vmax)
            elif ext=='.dat':
                over += [(get_dat(o), d)]
            else:
                raise NotImplementedError
        args.overplots = [over, rmin, rmax]
    return args

def get_fits(filename):
    from astropy.io import fits
    return fits.open(os.path.expanduser(filename))[0]

def get_dat(filename):
    import numpy as np
    return np.loadtxt(os.path.expanduser(filename))

def get_img_ranges(img):
    return auto_vminmax(img.data)

