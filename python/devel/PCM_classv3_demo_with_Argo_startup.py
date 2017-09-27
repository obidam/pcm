# -*coding: UTF-8 -*-
__author__ = 'gmaze@ifremer.fr'

# In[IMPORT]:
import os
import pprint
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy
import cartopy.crs as ccrs
from cartopy.examples.arrows import sample_data
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import seaborn as sns

import sys

sys.dont_write_bytecode = True  # This prevent to create .pyc files, useful for debug and reload


def discrete_colorbar(N, cmap, ticklabels='default'):
    """Add a colorbar with a discrete colormap.

        N: number of colors.
        cmap: colormap instance, eg. cm.jet.
    """
    # source: https://stackoverflow.com/questions/18704353/correcting-matplotlib-colorbar-ticks
    mappable = plt.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, N + 0.5)
    colorbar = plt.colorbar(mappable, shrink=0.5)
    colorbar.set_ticks(np.linspace(0, N, N))
    if 'default' in ticklabels:
        ticklabels = range(N)
    colorbar.set_ticklabels(ticklabels)
    return colorbar


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.
    """
    # source: https://stackoverflow.com/questions/18704353/correcting-matplotlib-colorbar-ticks
    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in xrange(N + 1)]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, N)


def init_map(region=[-180, 180, -90, 90], dpi=80, figsize=(7, 7), dxtick=60, dytick=30):
    """Create a cartopy map

        Return fig,ax,proj

        Default parameters:
            region=[-180,180,-90,90]
            dpi=80
            figsize=(7,7)
            dxtick=60
            dytick=30

        Example:
            fig,ax,proj = init_map()
            plt.show()

            fig,ax,proj = init_map(region=[-80,0,0,70],dxtick=10,dytick=5)
            plt.show()

        gmaze@ifremer.fr
    """

    # Init figures and cartopy projection
    fig = plt.subplots(nrows=1, ncols=1, figsize=figsize, dpi=dpi, facecolor='w', edgecolor='k')
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.set_extent(region)

    # Add the decorum:
    ax.coastlines(linewidth=0.1)
    # ax.stock_img()
    ax.add_feature(cartopy.feature.LAND)
    # ax.add_feature(cartopy.feature.OCEAN)
    #     ax.add_feature(cartopy.feature.LAND,facecolor=bgcolor,edgecolor='k')
    #     ax.add_feature(cartopy.feature.OCEAN,zorder=0)#,zorder=0,facecolor='k')

    # Add latitude/longitude lines:
    gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.xlabels_bottom = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabel_style = {'size': 6, 'color': 'black'}
    gl.ylabel_style = {'size': 6, 'color': 'black'}

    # Add latitude/longitude labels:
    xticks = np.arange(region[0], region[1] + np.finfo(np.float32).eps, dxtick)
    yticks = np.arange(region[2], region[3] + np.finfo(np.float32).eps, dytick)
    ax.set_xticks(xticks, crs=proj)
    ax.set_yticks(yticks, crs=proj)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # Output
    return fig, ax, proj


# In[Load pcmpack]:

import pcmpack
from pcmpack import v3 as pcm
dir(pcm)


# In[Load Argo data]
#%time
# Define where to find the Argo file on STD:
ncroot = '.'
ncroot = '/Users/gmaze/data/ARGO/copoda_db/setup_H/db_thd_config6_last/gmm'

# Load Argo Training dataset:
ncfile = 'NATL_HOMOGENEOUS_variables_7subset_1.nc'
dtrain = xr.open_mfdataset(os.path.join(ncroot, ncfile))

# Load all Argo dataset:
# ncfile = 'NATL_HOMOGENEOUS_variables_7subset_*.nc'
# dset = xarray.open_mfdataset(os.path.join(ncroot, ncfile))

# And Create the array X(Nz,Np): The field to classify with a GMM,
# Np profiles with Nz depth levels.
lon, lat = dtrain['LONGITUDE'], dtrain['LATITUDE']
X, Xunit = dtrain['TEMP'], dtrain['TEMP'].units
Z = dtrain['DEPTH']

# Size of the training set X:
[Np, Nz] = X.shape
print "Number of features (Depth Levels): ", Nz
print "Number of samples (N profiles): ", Np

