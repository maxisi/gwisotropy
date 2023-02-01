from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns

# Convert pts to inches
fig_width_pt = 2*246.0
inches_per_pt = 1.0/72.27               
# Golden ratio
fig_ratio = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height = fig_width*fig_ratio

figsize_column = (fig_width, fig_height)

def sky_scatter(x, ax=None, **kws):
    fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot(111, projection='mollweide')
    cnorm = matplotlib.colors.Normalize(vmin=0, vmax=1)

    vecs = x.reshape(np.prod(x.shape[:2]), 3)

    # idxs = random.choice(arange(len(nhats)), 512, replace=False)
    c = np.linalg.norm(vecs, axis=1)
    idxs = ...
    lat = np.arcsin(vecs[idxs,2]/c)
    lon = np.arctan2(vecs[idxs,1], vecs[idxs,0])
    c = np.linalg.norm(vecs, axis=1)
    def_kws = dict(marker='.', alpha=0.5, s=5, cmap='magma')
    def_kws.update(**kws)
    ax.scatter(lon, lat, c=c, norm=cnorm, **def_kws)
    return ax

def sky_hex(x, ax=None, **kws):
    fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot(111, projection='mollweide')
    cnorm = matplotlib.colors.Normalize(vmin=0, vmax=1)

    vecs = x.reshape(np.prod(x.shape[:2]), 3)

    # idxs = random.choice(arange(len(nhats)), 512, replace=False)
    c = np.linalg.norm(vecs, axis=1)
    idxs = ...
    lat = np.arcsin(vecs[idxs,2]/c)
    lon = np.arctan2(vecs[idxs,1], vecs[idxs,0])
    c = np.linalg.norm(vecs, axis=1)
    def_kws = dict(cmap='magma')
    def_kws.update(**kws)
    ax.hexbin(lon, lat, C=c, norm=cnorm, **def_kws)
    return ax

def sky_density(x, cmap='viridis', bin_number=15):
    import cartopy.crs as ccrs
    vecs = x.reshape(np.prod(x.shape[:-1]), 3)
    
    # idxs = random.choice(arange(len(nhats)), 512, replace=False)
    c = np.linalg.norm(vecs, axis=1)
    lat = np.arcsin(vecs[:,2]/c)
    lon = np.arctan2(vecs[:,1], vecs[:,0])
    
    lon_edges = np.linspace(-np.pi, np.pi, bin_number + 1)
    lat_edges = np.arcsin(np.linspace(-1., 1., bin_number + 1))

    # calculate 2D histogram, the shape of hist is (bin_number, bin_number)
    data, lon_edges, lat_edges = np.histogram2d(lon, lat, 
                                                bins=[lon_edges, lat_edges],
                                                density=True)
    # compute bin centers
    x = 0.5*(lon_edges[1:] + lon_edges[:-1])
    y = 0.5*(lat_edges[1:] + lat_edges[:-1])

    X, Y = np.degrees(np.meshgrid(x, y))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())

    ax.contourf(X, Y, data.T, levels=5,
                transform=ccrs.PlateCarree(),
                cmap=cmap)
    ax.set_global()
    return ax

def add_colorbar(ax, key='N', cmap='viridis'):
    cmap = 'viridis'
    color_ticks = [-1, 0, 1]
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    cm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cm.set_array([])

    cbaxes = ax.inset_axes([0.7, 1, 0.3, 0.05])

    cb = plt.colorbar(cm, orientation='horizontal', cax=cbaxes, ticks=color_ticks)
    # cb.ax.set_xticklabels([str(v) for v in color_ticks], fontsize=20)
    cbaxes.set_xlabel(r"$|v_%s|$" % key, fontdict={'size': 12}, labelpad=-6)
    cbaxes.xaxis.set_ticks_position('top')
    cbaxes.xaxis.set_label_position('top')
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params('x', length=0, width=0.5, which='major', labelsize=10)
    return cb
