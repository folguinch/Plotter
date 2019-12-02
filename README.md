## Usage

- Global command line options:

```bash
usage: plotter.py [-h] [--loglevel LOGLEVEL] [--png] [--axlabel] [--legend]
                  [--overplot [OVERPLOT [OVERPLOT ...]]]
                  [--oplabel [LABEL [LABEL ...]]]
                  [--opcolor [COLOR [COLOR ...]]]
                  [--shape ROWS COLS | --rows ROWS | --cols COLS]
                  [--detailaxlabel | --axlabels [AXLABELS [AXLABELS ...]]]
                  [--overall | --opmapping [OPMAPPING [OPMAPPING ...]]]
                  config plotname {fromcube,plotmaps,chanmap,plotmoments} ...

Data plotting tools.

positional arguments:
  config                Plot configuration file
  plotname              Plot file name
  {fromcube,plotmaps,chanmap,plotmoments}
    fromcube            Plot data from cube (moments 0, 1, and channel maps)
    plotmaps            Plot FITS maps
    chanmap             Plot channel maps
    plotmoments         Plot moment maps

optional arguments:
  -h, --help            show this help message and exit
  --loglevel LOGLEVEL   Console logging level (default: ['info'])
  --png                 Save also as png (default: False)
  --axlabel             Automatic select axis labels (default: False)
  --legend              Plot legend (default: False)
  --overplot [OVERPLOT [OVERPLOT ...]]
                        Overplot files (extension dependent) (default: None)
  --oplabel [LABEL [LABEL ...]]
                        Overplot label for legend (default: None)
  --opcolor [COLOR [COLOR ...]]
                        Overplot colors (default: None)
  --shape ROWS COLS     Shape of figure tile (default: None)
  --rows ROWS           Number of rows (default: None)
  --cols COLS           Number of columns (default: None)
  --detailaxlabel       Automatic detailed axis labels (default: False)
  --axlabels [AXLABELS [AXLABELS ...]]
                        Axis labels (one per each image) (default: None)
  --overall             Plot all overplots over each image (default: False)
  --opmapping [OPMAPPING [OPMAPPING ...]]
                        Indices of axes where to plot the overplots (default:
                        None)
```

- `fromcube` options:

``` bash
usage: plotter.py config plotname fromcube [-h] [--continuum CONTINUUM]
                                           [--every STEP]
                                           [--save_moments SAVE_MOMENTS]
                                           [--line_list [LINE_LIST [LINE_LIST ...]]]
                                           [--chanran INITIAL FINAL | --chans [CHAN [CHAN ...]]]
                                           lineconfig cubename

positional arguments:
  lineconfig            Configuration file of the desired lines
  cubename              Cube file name

optional arguments:
  -h, --help            show this help message and exit
  --continuum CONTINUUM
                        Background image file for color plot
  --every STEP          Step between channels for range
  --save_moments SAVE_MOMENTS
                        Moments base file name
  --line_list [LINE_LIST [LINE_LIST ...]]
                        Lines to plot (default all in config)
  --chanran INITIAL FINAL
                        Channel range to plot (final not included)
  --chans [CHAN [CHAN ...]]
                        List of channels to plot
```

- `plotmaps` options:

```bash
usage: plotter.py config plotname plotmaps [-h] [--section SECTION]
                                           [--selflevels]
                                           [imagenames [imagenames ...]]

positional arguments:
  imagenames         List of images to plot

optional arguments:
  -h, --help         show this help message and exit
  --section SECTION  Section of the config file
  --selflevels       Plot contours from input image
```

- `chanmap` options:

```bash
usage: plotter.py config plotname chanmap [-h] [--auto_velshift | --vlsr VLSR]
                                          [--section SECTION]
                                          [--continuum CONTINUUM]
                                          [--every STEP]
                                          (--chanran INITIAL FINAL | --chans [CHAN [CHAN ...]])
                                          cubename

positional arguments:
  cubename              Cube file name

optional arguments:
  -h, --help            show this help message and exit
  --auto_velshift       Determine the velocity shift automatically
  --vlsr VLSR           LSR velocity to shift spectral data in km/s
  --section SECTION     Section of the config file
  --continuum CONTINUUM
                        Background image file for color plot
  --every STEP          Step between channels for range
  --chanran INITIAL FINAL
                        Channel range to plot (final not included)
  --chans [CHAN [CHAN ...]]
                        List of channels to plot
```

- `plotmoments` options:

```bash
usage: plotter.py config plotname plotmoments [-h]
                                              [--auto_velshift | --vlsr VLSR]
                                              [--section SECTION]
                                              [--moments [MOMENTS [MOMENTS ...]]]
                                              --imagenames IMAGENAMES
                                              [IMAGENAMES ...]

optional arguments:
  -h, --help            show this help message and exit
  --auto_velshift       Determine the velocity shift automatically
  --vlsr VLSR           LSR velocity to shift spectral data in km/s
  --section SECTION     Section of the config file
  --moments [MOMENTS [MOMENTS ...]]
                        Number of the moment of the corresponding image
  --imagenames IMAGENAMES [IMAGENAMES ...]
                        List of images to plot
```
