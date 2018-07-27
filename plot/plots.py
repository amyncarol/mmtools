def scatter_hist(x, y, z, x_label, y_label, z_label):
    """
    Plot scatter plot x vs y and also histograms for x and y, respectively. 
    Also shows the values of z as colors of the scatter points.

    Args:
    x, y, z: data
    x_label, y_label, z_label: labels

    Returns:
    the figure

    """
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # the scatter plot:
    axScatter.scatter(x, y, c=z/max(z), s = max(z)*50, alpha=0.5, cmap='gist_rainbow')
    axScatter.set_xlim((min(x), max(x)))
    axScatter.set_ylim((min(y), max(y)))
    axScatter.set_xlabel(x_label)
    axScatter.set_ylabel(y_label)
    
    # create colorbar for z dimension
    norm = mpl.colors.Normalize(vmin=0, vmax=z.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap='gist_rainbow')
    cmap.set_array([])
    cbaxes = inset_axes(axScatter, width="30%", height="3%", loc=1)
    cbar = plt.colorbar(cmap, cax=cbaxes, ticks=np.arange(0, max(z), 0.5), orientation='horizontal')
    cbar.set_label(z_label)
    
    # histograms
    binwidth = 0.5
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth
    bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(x, bins=bins, alpha=0.5)
    axHisty.hist(y, bins=bins, orientation='horizontal', alpha=0.5)
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    axHistx.set_ylabel('Count')
    axHisty.set_xlabel('Count')
    axHistx.set_xticks([])
    axHisty.set_yticks([])
    
    return fig