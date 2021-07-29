import numpy as np
import matplotlib.transforms as transforms
import matplotlib.patches as mpatches
from spec_tools import air_conversion


def plot_window(elem, ax, color='lightgray', air=False, path='./windows/'):

    # Calculate the vacuum wavelengths of each pixel in the APOGEE chips
    chip_start = [4.180476, 4.200510, 4.217064]
    n_pixels = [3028, 2495, 1991]

    wave_a = chip_start[0] + np.arange(3028) * 6e-6
    wave_b = chip_start[1] + np.arange(2495) * 6e-6
    wave_c = chip_start[2] + np.arange(1991) * 6e-6
    wave_filt = np.concatenate((np.concatenate((wave_a, wave_b)), wave_c))
    wave_filt = 10.**wave_filt

    if air:
        wave_filt = air_conversion(wave_filt)

    # Read in the window strengths from the mask files
    window_strength = np.genfromtxt(f'{path}{elem}.mask')

    # Only select the pixels where the window strength is greater than 0
    x_window = wave_filt[window_strength > 0.]

    # Set up an axis transformation to plot the windows with x values on the data scale (wavelength)
    # but have the y axis tied to the axes (so we can plot across the vertical range of the figure)
    trans = transforms.blended_transform_factory(
        ax.transData, ax.transAxes)

    # Set initial values to make the window finding work and set up the shapes list
    # (each window will end up a shape in this list)
    end_point = -1
    width = 0.3

    # Loop through the pixels with window strength > 0 and if the pixels are adjacent expand the current window
    # if the next pixel is not adjacent to the last, save the last window as an axis patch and add it to
    # the plot
    for i in range(len(x_window)):
        if i < end_point:
            continue
        for j in range(len(x_window) - 1 - i):
            if x_window[j + 1 + i] - x_window[i + j] < 0.3:
                width = x_window[j + 1 + i] - x_window[i]
                if len(x_window) - 1 == j + 1 + i:
                    rect = mpatches.Rectangle((x_window[i] - 0.14, 0), width=width + 0.28, height=1,
                                              transform=trans, color=color, ec='None',
                                              alpha=0.5, zorder=1)
                    ax.add_patch(rect)
                    width = 0.3
            else:
                rect = mpatches.Rectangle((x_window[i] - 0.14, 0), width=width + 0.28, height=1,
                                          transform=trans, color=color, ec='None',
                                          alpha=0.5, zorder=1)
                ax.add_patch(rect)
                width = 0.3
                break

        end_point = i + j + 1

        if end_point > len(x_window[:]) - 1:
            break


def plotly_window(elem, fig, shapes=None, color='lightgray', air=False, path='./windows/'):

    # Calculate the vacuum wavelengths of each pixel in the APOGEE chips
    chip_start = [4.180476, 4.200510, 4.217064]
    n_pixels = [3028, 2495, 1991]

    wave_a = chip_start[0] + np.arange(3028) * 6e-6
    wave_b = chip_start[1] + np.arange(2495) * 6e-6
    wave_c = chip_start[2] + np.arange(1991) * 6e-6
    wave_filt = np.concatenate((np.concatenate((wave_a, wave_b)), wave_c))
    wave_filt = 10.**wave_filt

    if air:
        wave_filt = air_conversion(wave_filt)

    if shapes is None:
        shapes = []

    # Read in the window strengths from the mask files
    window_strength = np.genfromtxt(f'{path}{elem}.mask')

    # Only select the pixels where the window strength is greater than 0
    x_window = wave_filt[window_strength > 0.]

    # Set initial values to make the window finding work and set up the shapes list
    # (each window will end up a shape in this list)
    end_point = -1
    width = 0.3

    # Loop through the pixels with window strength > 0 and if the pixels are adjacent expand the current window
    # if the next pixel is not adjacent to the last, save the last window as an axis patch and add it to
    # the plot
    for i in range(len(x_window)):
        if i < end_point:
            continue
        for j in range(len(x_window) - 1 - i):
            if x_window[j + 1 + i] - x_window[i + j] < 0.3:
                width = x_window[j + 1 + i] - x_window[i]
                if len(x_window) - 1 == j + 1 + i:
                    shapes.append(dict(
                        type="rect",
                        xref="x",
                        yref="paper",
                        x0=x_window[i] - 0.14,
                        y0=0,
                        x1=x_window[i] + 0.14 + width,
                        y1=1,
                        fillcolor=color,
                        opacity=0.5,
                        line_width=0,
                        layer='below'
                    ))
                    width = 0.3
            else:
                shapes.append(dict(
                    type="rect",
                    xref="x",
                    yref="paper",
                    x0=x_window[i] - 0.14,
                    y0=0,
                    x1=x_window[i] + 0.14 + width,
                    y1=1,
                    fillcolor=color,
                    opacity=0.5,
                    line_width=0,
                    layer='below'
                ))
                width = 0.3
                break

        end_point = i + j + 1

        if end_point > len(x_window[:]) - 1:
            break

    rect = fig.update_layout(shapes=shapes)
    return shapes
