import os

def postprocess(fig, args):
    # Save figure
    plotname = os.path.expanduser(args.plotname[0])
    fig.savefig(plotname)
    if args.png:
        plotname = os.path.splitext(plotname)[1] + '.png'
        fig.savefig(plotname)

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
    
    return rows, cols
