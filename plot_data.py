import argparse

import loaders
import src.plotter as plt
import utils 

def plot_data_parser():
    # Parser help
    h = "Plot data files"

    # Parser
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--section', nargs=1, default=['data_plots'],
            help="Section of the config file")
    parser.add_argument('--columns', metavar=['COL1','COL2'], nargs=2, default=[0,1],
            help="Data columns to plot")
    parser.add_argument('--errcols', metavar=['COL'], nargs=1, default=[None],
            help="Error column")
    parser.add_argument('filenames', nargs='*',
            help='Data files (1 set per plot, else use overplot)')
    parser.set_defaults(func=plot_data, loaders=[loaders.load_data])

    return {'plotdata': (parser, h)}

def plot_data(args):
    # Keyword arguments for tile plotter
    opts = {}
    opts['nrows'], opts['ncols'] = utils.get_shape(args, len(args.data),
            default_cols=1)
    assert opts['nrows']*opts['ncols'] >= len(args.data)
    args.logger.info('Rows, columns = %i, %i', opts['nrows'], opts['ncols'])

    # Setup tile plotter
    args.logger.debug('Initializing figure')
    fig = plt.NPlotter(config=args.config[0], section=args.section[0],
            **opts)

    # Iterate over data
    xunits = []
    yunits = []
    for i,(loc, data) in enumerate(zip(fig.axes, args.data)):
        label = fig.get_value('axlabel', None, loc, sep=',')
        label = utils.get_axis_label(args, i, label)
        overplots = utils.get_overplots(args, i)
        
        ax, xunit, yunit = utils.splot(loc, fig, data, args.logger, 
                overplots=overplots, cols=args.columns, 
                errorcol=args.errcols[0])

        xunits += [xunit]
        yunits += [yunit]

    fig.auto_config(legend=args.legend, units=(xunits, yunits))
        
    return fig
    
