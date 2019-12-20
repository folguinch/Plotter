# Plotter

Set of tools to plot different types of data.

## Download

```bash
git clone --recurse-submodules git@github.com:folguinch/plotter.git
```

## Requirements

- Numpy
- Matplotlib
- astropy
- spectral_cube

Matplotlib styles are available in `src/styles`. These are useful to define 
color schemes for e.g. map plots. To use these, copy or link them into
`~/.config/matplotlib/stylelib/`

## Description

At the moment 4 types of plots have been implemented:
- `plotmaps`: plot (2-D) images from FITS files.
- `chanmaps`: channel maps from a data cube.
- `plotmoments`: plot cube moments.
- `fromcube`: shortcut for plotting moments 0, 1 and chanmaps from the same line.

Each plot is configured from a configuration file, although some options can be
overwritten from the command line. This file configures the geometry (including 
color bars) of the plot as well as labels and markers. All the options available 
for configuring the geometry are available in `src/examples/default.cfg`. 
Additional arguments in the configuration file are used to confire labels, 
markers, legends, etc. See the directory `examples`.

## Plotting maps

### Configuration file parameters

When making a tile plot, several properties may be specific for certain plots. 
Global parameter for all the plots in the tile should be given by writing:
```INI
parameter_name = value1 value2 ...
```
Note that some parameters, like coordinates, are coma separated, e.g:
```INI
markers = ra1 dec1, ra2 dec2, ...
```
In general parameters that involve text are coma separated.

Parameter for a specfic `row_number, column_number` pair can be given by:
```INI
parameter_nameIJ = value1 value2 ...
```
where `I = row_number` and `J = column_number`.

*IMPORTANT:* Note that the behaviour may depend on the parameter, e.g. if 
there is a list of `rms` values each one will correspond to one plot in the 
tile. In this case, the indexing is from left to right and the from top to 
bottom, and starting at 0, e.g. the index of `row, column = 1, 0` of a 
3-by-3 tile is `3`.

### Levels for contours

With the exception of channel maps, contours can be from the background image,
using the command line option `self_contours`, and/or from other images, using 
the `overplot` option.

There are 2 ways of determining the contour levels:
- directly from the `levels` parameter in the configuration file,
- or rms based contours. 
The `rms` can be given in the configuration file together with the number of 
rms values for the first level (`nsigma`). The default value for `nsigma` is 5.
The parameter `nlevels` is used when the `global_levels` option is used and the 
`self_contours` is not, but this requires an input `rms` value to determine the 
levels.  

Hint: when the `overplot`s come from different frequencies the `global_levels` 
should not be used, while when plotting models fitting the data `global_levels`
(and probably `self_contours`) are useful to compare like with like. In the 
former, each overplot will have its own levels. In the latter, if 
`self_contours` is not used then `rms` and `nlevels` values should be given.

### For plotting maps

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
