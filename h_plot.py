import sys
import os.path
from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
from matplotlib import rc
from matplotlib import colors
import matplotlib.patches as patch
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LogFormatterExponent
import matplotlib.colorbar as colorbar
import colormaps as cmaps
from matplotlib.ticker import ScalarFormatter


def density2d_to_points(x, y, ax=None, bins=150, masklim=5, xlim=None,
                        ylim=None, cmap='viridis', norm_func=colors.LogNorm,
                        interpolation='nearest', label=None,
                        zorder=1, point_color='k', point_size=2):
    '''Plot a 2D histogram over points when the density of points is above
    a given threshold.

    Parameters
    ----------
    x : array_like shape (N,)
      An array containing the x values to be plotted
    y : array_like shape (N,)
      An array containing the y values to be plotted
    ax : matplotlib figure axes
      The matplotlib figure axes to plot on.
    bins : int or array_like or [int, int] or [array, array], optional
      The bin specification for the 2d histogram, following the
      parameterization set by numpy.histogram2d().  Defaults to 150.
    masklim : int
      The number threshold of points per bin, below which 2d histogram pixels
      are masked to only show individual data points
    xlim : array_like, shape(2,), optional
      The x range to be included in the 2d histogram
    ylim : array_like, shape(2,), optional
      The y range to be included in the 2d histogram
    cmap : str
      The matplotlib color map to use for the 2d histogram
    norm : matplotlib.colors normalization
      The color normalization routine from matplotlib.  Defaults to log
      normalization.
    interpolation : str
      The interpolation scheme used to plot the 2d histogram, see
      matplotlib.imshow() for more details.  Defaults to nearest.
    label : str
      The label for the data points.
    zorder : int
      The plotting order for the plotted distribution.  Defaults to 1

    Returns
    -------
    points : matplotlib line2D collection
      Matplotlib collection of the points that are plotted
    cset : matplotlib patches collection
      Matplotlib collection describing the plotted 2d histogram

    See also
    --------
    numpy.histogram2d : numpy 2d histogram function
    matplotlib.pyplot.plot : matplotlib line and point plotting function
    matplotlib.pyplot.imshow : matplotlib image plotting function
    '''

    check_axes(ax)
    check_data(x)
    check_data(y)
    range = [set_range(x, xlim), set_range(y, ylim)]

    points = ax.plot(x, y, color=point_color, marker='.', ms=point_size,
                     ls='none', zorder=zorder, rasterized=True, label=label)

    H, xedges, yedges = np.histogram2d(x, y, bins=bins, range=range)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    masked_H = np.ma.masked_less(H, masklim).transpose((1, 0))
    if np.isfinite(np.nanmin(masked_H)) and np.isfinite(np.nanmax(masked_H)):
        norm = norm_func(np.nanmin(masked_H), np.nanmax(masked_H))
    else:
        norm = None

    cset = ax.imshow(masked_H, extent=extent, cmap=cmap,
                     origin='lower', interpolation=interpolation, norm=norm,
                     zorder=zorder, aspect='auto')

    return points, cset


def log_colorbar(cset, fig=None, ax=None, location="right", size="5%",
                 pad=0.05, format="%d", labelsize=13, fontsize=18,
                 label='Num'):
    ''' Add a log scale colorbar to the input axes, with properties set according
    to the matplotlib colorbar framework.

    Parameters
    ----------
    cset : patch set image
      This is the colored set that the colorbar will be illustrating.
    ax : matplotlib axes
      Axes to attach the colorbar to.  Will return error if no axes are given
    location : str
      String defining the location of the colorbar with respect to the axes (
      "top", "bottom", "right", "left").  Defaults to "right".
    size : str
      Percent width of the colorbar compared to the axes size.  Defaults to 5%
    pad : float
      Fraction of image to use as the axes padding for the colorbar.  Defaults
      to 0.05.
    format : str
      Format for the colorbar tick labels.  Defaults to "%d".
    labelsize : float
      Size of the tick labels in pts.  Defaults to 13 pt.
    fontsize : float
      Size of the colorbar label in pts.  Defaults to 18 pt.
    label : str
      Title for the colorbar.  Defaults to "Num".

    Returns
    -------
    cbar : matplotlib colorbar
      Matplotlib collection of the points that are plotted
    divider : mpl location of axes "ax"
      Matplotlib collection describing the plotted 2d histogram
    cax : matplotlib axes appended to axes "ax"
      Matplotlib axes where the colorbar "cbar" is attached to the figure

    See also
    --------
    mpl_toolkits.axes_grid1.make_axes_locatable : matplotlib function
      to locate axes and attach new ones
    matplotlib.pyplot.colorbar : matplotlib colorbar function
    matplotlib.ticker.LogFormatterExponent : matplotlib log formatting function
    '''

    check_axes(ax)
    check_figure(fig)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes(location, size=size, pad=pad)
    cbar = fig.colorbar(cset, cax=cax, format=format)
    cbar.formatter = LogFormatterExponent(base=10)
    cbar.ax.tick_params(labelsize=labelsize)
    cbar.set_label(label, fontsize=fontsize)
    cbar.ax.tick_params(which="both", axis="y", direction="in")

    return cbar, divider, cax


def h_colorbar(cset, fig=None, ax=None, location="right",
               orientation='vertical', size="5%", pad=0.05, format="%d",
               labelsize=13, fontsize=18, label='Num',
               log=True):
    ''' Add a log scale colorbar to the input axes, with properties set according
    to the matplotlib colorbar framework.

    Parameters
    ----------
    cset : patch set image
      This is the colored set that the colorbar will be illustrating.
    ax : matplotlib axes
      Axes to attach the colorbar to.  Will return error if no axes are given
    location : str
      String defining the location of the colorbar with respect to the axes (
      "top", "bottom", "right", "left").  Defaults to "right".
    orientation : str
      String defining the orientation of the scale in the colorbar ("vertical",
      "horizontal").  Defaults to "vertical".
    size : str
      Percent width of the colorbar compared to the axes size.  Defaults to 5%
    pad : float
      Fraction of image to use as the axes padding for the colorbar.  Defaults
      to 0.05.
    format : str
      Format for the colorbar tick labels.  Defaults to "%d".
    labelsize : float
      Size of the tick labels in pts.  Defaults to 13 pt.
    fontsize : float
      Size of the colorbar label in pts.  Defaults to 18 pt.
    label : str
      Title for the colorbar.  Defaults to "Num".
    log : boolean
      If true will make the colorbar logarithmic

    Returns
    -------
    cbar : matplotlib colorbar
      Matplotlib collection of the points that are plotted
    divider : mpl location of axes "ax"
      Matplotlib collection describing the plotted 2d histogram
    cax : matplotlib axes appended to axes "ax"
      Matplotlib axes where the colorbar "cbar" is attached to the figure

    See also
    --------
    mpl_toolkits.axes_grid1.make_axes_locatable : matplotlib function
      to locate axes and attach new ones
    matplotlib.pyplot.colorbar : matplotlib colorbar function
    matplotlib.ticker.LogFormatterExponent : matplotlib log formatting function
    '''

    check_axes(ax)
    check_figure(fig)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes(location, size=size, pad=pad)
    cbar = fig.colorbar(cset, cax=cax, format=format, orientation=orientation)
    if log:
        formatter = ScalarFormatter()
        # formatter.set_powerlimits((0,2))
        formatter.set_scientific(False)
        cbar = fig.colorbar(cset, cax=cax, format=formatter,
                            orientation=orientation)
    else:
        cbar = fig.colorbar(cset, cax=cax, format=format,
                            orientation=orientation)
    cbar.set_label(label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=labelsize)
    if location == 'right':
        cax.yaxis.set_ticks_position(location)
        cbar.ax.tick_params(which="both", axis="y", direction="in")
    elif location == 'left':
        cax.yaxis.set_ticks_position(location)
        cax.yaxis.set_label_position('left')
        cbar.ax.tick_params(which="both", axis="y", direction="in")
    elif location == 'top':
        cax.xaxis.set_ticks_position(location)
        cax.xaxis.set_label_position('top')
        cbar.ax.tick_params(which="both", axis="x", direction="in")
    elif location == 'bottom':
        cax.xaxis.set_ticks_position(location)
        cax.xaxis.set_label_position('bottom')
        cbar.ax.tick_params(which="both", axis="x", direction="in")

    return cbar, divider, cax


def set_range(a, range):
    '''Checks the given range and if none is given, sets range to a.min(),
    a.max()

    Parameters
    ----------
    a : array_like shape (N,)
      Data to be plotted.
    range : array_like, shape(2,), optional
      range to set

    Returns
    -------
    range : array_like, shape(2,2)
      returns the manually or automatically assigned (x,y) range for the input
      data

    '''

    if range is None:
        lower_limit, upper_limit = a.min(), a.max()
        if not (np.isfinite(lower_limit) and np.isfinite(upper_limit)):
            raise ValueError(
                "autodetected range of [{}, {}] is not finite".format(
                    lower_limit, upper_limit))
    elif len(range) != 2:
        raise ValueError(
            'supplied range must be in the form [lower_limit, upper_limit]')
    elif len(range) == 2:
        lower_limit, upper_limit = range
        if not (np.isfinite(lower_limit) and np.isfinite(upper_limit)):
            raise ValueError(
                "supplied range of [{}, {}] is not finite".format(
                    lower_limit, upper_limit))
        if lower_limit > upper_limit:
            raise ValueError(
                'max must be larger than min in range parameter.')
    return [lower_limit, upper_limit]


def check_axes(ax):
    '''Checks to verify that figure axes are given and throws an error otherwise

    Parameters
    ----------
    ax : matplotlib figure axes
      The matplotlib figure axes to plot on.

    Returns
    -------
    pass : None
      passes upon successful checking of the presence of axes

    '''
    if ax is None:
        raise TypeError('must provide axes to plot on.')
    pass


def check_figure(fig):
    '''Checks to verify that figure axes are given and throws an error otherwise

    Parameters
    ----------
    fig : matplotlib figure
      The matplotlib figure to plot with.

    Returns
    -------
    pass : None
      passes upon successful checking of the presence of a figure

    '''

    if fig is None:
        raise TypeError('must provide figure to plot on.')
    pass


def check_data(a):
    '''Checks to verify that figure axes are given and throws an error otherwise

    Parameters
    ----------
    a : array_like shape (N,)
      Data to be plotted.

    Returns
    -------
    pass : None
      passes upon successful checking of the presence of data

    '''
    if len(a) < 1:
        raise ValueError('data array must have a minimum length of 1.')
    pass
