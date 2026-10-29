from matplotlib import pyplot as plt
import numpy as np
from .spectrum import Spectrum
from .spectral_lines import get_wavelength, format_line

from astropy.nddata import block_reduce
from astropy.nddata import CCDData



def _get_spectra_xy(spec):
    if isinstance(spec, tuple):
        if len(spec) == 2:
            x, y = spec
            w = None
        else:
            x, y, w = spec
    elif isinstance(spec, Spectrum):
        x, y = spec.wavelength, spec.flux
        if spec.uncertainty is not None:
            w = spec.uncertainty
    else:
        raise TypeError(f"Unknown spectral type: {type(spec)}")

    return x, y, w


def plot_wrapped_spectrum(inputs, styles=None, *, 
        xlim = None,
        N_chunks = 8,
        redshift = None,
        label_lines = None,
        xlabel = "wavelength", 
        ylabel = "flux", 
        show_uncertainty = False,
        **kwargs):

    if styles is None:
        styles = {}
    if label_lines is None:
        label_lines = []

    if xlim is None:
        first_spectrum = next(iter(inputs.values()))
        x = _get_spectra_xy(first_spectrum)[0]
        xmin = np.min(x) 
        xmax = np.max(x)
        dx_tot = xmax - xmin
        pad = 0.001
        xmin -= pad * dx_tot
        xmax += pad*dx_tot
        xlim = (xmin, xmax)


    dx = np.diff(xlim)[0] / N_chunks

    fig, axs = plt.subplots(N_chunks, 1, figsize=(16, 2*N_chunks))

    for i in range(N_chunks):
        ax = axs[i]

        for label, spec in inputs.items():
            x, y, w = _get_spectra_xy(spec)
            style = styles.get(label, {})
            xlim_frame = xlim[0] + i*dx, xlim[0] + (1+i) * dx
            filt = x >= xlim_frame[0]
            filt &= x < xlim_frame[1]
            p = ax.plot(x[filt], y[filt], label=label, **style)

            if show_uncertainty and w is not None:
                ax.fill_between(x[filt], y[filt] - w[filt], y[filt] + w[filt],
                        alpha=0.2, color=p[0].get_color())
        
        if i == N_chunks - 1:
            xlabel_this = xlabel
        else:
            xlabel_this = ""

        plot_line_labels(spec, label_lines, redshift, ax=ax)
        ax.set(
            xlim = xlim_frame,
            xlabel = xlabel_this, 
            ylabel = ylabel, 
            **kwargs
        )

        if i == N_chunks - 1:
            ax.legend()
        

    fig.tight_layout()


def plot_line_labels(spec, lines, redshift, ax=None):
    if redshift is None:
        return
    if ax is None:
        ax = plt.gca()

    specwave, flux, _ = _get_spectra_xy(spec)
        
    for line in lines:
        wave = get_wavelength(line)
        wave *= 1 + redshift
            
        idx = np.argmin(np.abs(specwave - wave))
        y0 = np.max(flux[max(0, idx-5):min(idx+5, len(specwave))])

        ax.annotate(format_line(line), (wave, y0), (0, 2), 
                     xycoords="data", textcoords="offset fontsize", 
                     fontsize=6, rotation=90, ha="center",
                     arrowprops=dict(arrowstyle="-", linewidth=1, color="black"),
                    )



def show_colour_image(image,
               figsize=(10, 10),
               dpi = None,
               clabel=None,
               clim=None,
               show_colorbar=True, show_ticks=True,
               fig=None, ax=None, input_ratio=None,
              subplot_kw=dict()):
    """
    Show an image in matplotlib with some basic astronomically-appropriate stretching.

    Parameters
    ----------
    image
        The image to show
    percl : number
        The percentile for the lower edge of the stretch (or both edges if ``percu`` is None)
    percu : number or None
        The percentile for the upper edge of the stretch (or None to use ``percl`` for both)
    figsize : 2-tuple
        The size of the matplotlib figure in inches
    """

    if dpi is None:
        dpi = image.shape[0] / figsize[0]

    if (fig is None and ax is not None) or (fig is not None and ax is None):
        raise ValueError('Must provide both "fig" and "ax" '
                         'if you provide one of them')
    elif fig is None and ax is None:
        if figsize is not None:
            # Rescale the fig size to match the image dimensions, roughly
            image_aspect_ratio = image.shape[0] / image.shape[1]
            figsize = (max(figsize) * image_aspect_ratio, max(figsize))

        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi, subplot_kw=subplot_kw)


    # To preserve details we should *really* downsample correctly and
    # not rely on matplotlib to do it correctly for us (it won't).

    # So, calculate the size of the figure in pixels, block_reduce to
    # roughly that,and display the block reduced image.

    # Thanks, https://stackoverflow.com/questions/29702424/how-to-get-matplotlib-figure-size

    fig_size_pix = fig.get_size_inches() * fig.dpi

    ratio = (image.shape[0:2] // fig_size_pix).max()

    if ratio <= 1:
        ratio = 1
    else:
        print("downsampling image by", ratio)

    ratio = input_ratio or ratio

    reduced_data = block_reduce(image, ratio)

    reduced_data = (reduced_data - clim[0]) / (clim[1] - clim[0])

    reduced_data[reduced_data > 1] = 1
    reduced_data[reduced_data < 0] = 0
    
    # Divide by the square of the ratio to keep the flux the same in the
    # reduced image. We do *not* want to do this for images which are
    # masks, since their values should be zero or one.
    # reduced_data = reduced_data / ratio**2

    # Of course, now that we have downsampled, the axis limits are changed to
    # match the smaller image size. Setting the extent will do the trick to
    # change the axis display back to showing the actual extent of the image.
    extent = [0, image.shape[1], 0, image.shape[0]]


    im = ax.imshow(reduced_data, origin='lower',
                   extent=extent, aspect='equal')


    if not show_ticks:
        ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)



def color_clip_image(img, vmin=0, vmax=None, norm=np.arcsinh, scale=1):
    img_clipped = img.copy()
    I = np.sum(img, axis=2) / 3
    if vmax is None:
        h = np.quantile(I[np.isfinite(I)], 0.99)
    else:
        h = vmax
    l = vmin

    def f(I):
        x = (norm(I/scale) - norm(l/scale)) / (norm(h/scale) - norm(l/scale))
        return x
    for i in range(3):
        img_clipped[:, :, i] *= f(I) / I

    return img_clipped


def _get_img_data(img):
    img_data = np.moveaxis(img.data, 0, -1)[:, :, ::-1]

    return img_data




def density_scatter(
    x, y, bins=300, xlim=None, ylim=None, ax=None,
    cmap=None, norm="log", scatter_below=3,
    s=0.5, point_color="k", colorbar=False, cbar_label="count",
    xlabel=None, ylabel=None, **hist_kwargs,
):
    """
    2D-histogram density plot, with optional scatter points in sparse regions.

    Parameters
    ----------
    x, y : array-like
        Data. Non-finite values are dropped.
    bins : int, [int, int], or array(s) of edges
        Passed to hist2d.
    xlim, ylim : (min, max), optional
        Histogram range and axis limits. Defaults to the data range.
    ax : matplotlib Axes, optional
    cmap, norm : colormap and normalization ("log", "linear", or a Normalize).
    scatter_below : int, optional
        If given, bins with fewer than this many points are left blank and
        their points are drawn individually instead.
    s, point_color : marker size and color for those scatter points.
    colorbar : bool
        Add a colorbar.
    **hist_kwargs : passed to ax.hist2d (e.g. vmin, vmax, rasterized).

    Returns
    -------
    counts, xedges, yedges, image
    """
    if ax is None:
        ax = plt.gca()

    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if x.shape != y.shape:
        raise ValueError(f"x and y must have the same length ({x.size} vs {y.size})")

    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if x.size == 0:
        raise ValueError("no finite data points to plot")

    if xlim is None:
        xlim = (x.min(), x.max())
    if ylim is None:
        ylim = (y.min(), y.max())

    in_range = (x >= xlim[0]) & (x <= xlim[1]) & (y >= ylim[0]) & (y <= ylim[1])
    x, y = x[in_range], y[in_range]

    rng = [xlim, ylim]

    if scatter_below is not None:
        # Find the count of the bin each point falls in
        H, xe, ye = np.histogram2d(x, y, bins=bins, range=rng)
        ix = np.clip(np.searchsorted(xe, x, side="right") - 1, 0, H.shape[0] - 1)
        iy = np.clip(np.searchsorted(ye, y, side="right") - 1, 0, H.shape[1] - 1)
        sparse = H[ix, iy] < scatter_below
        hist_kwargs.setdefault("cmin", scatter_below)
        ax.scatter(x[sparse], y[sparse], s=s, c=point_color,
                   lw=0, rasterized=True, zorder=0)

    counts, xedges, yedges, im = ax.hist2d(
        x, y, bins=bins, range=rng, cmap=cmap, norm=norm, zorder=1, **hist_kwargs
    )

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if colorbar:
        plt.colorbar(im, ax=ax, label=cbar_label)

    #return counts, xedges, yedges, im
    return
