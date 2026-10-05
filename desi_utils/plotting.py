from matplotlib import pyplot as plt
import numpy as np
from .spectrum import Spectrum
from .spectral_lines import get_wavelength, format_line


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
        raise Error(f"Unknown spectral type: {type(spec)}")

    return x, y, w


def plot_wrapped_spectrum(inputs, styles=dict(), *, 
        xlim = None,
        N_chunks = 8,
        redshift=None,
        label_lines = [],
        xlabel = "wavelength", ylabel="flux", 
        show_uncertainty = False,
        **kwargs):

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
            style = styles.get(label, dict())
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

    specwave, flux, w = _get_spectra_xy(spec)
        
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
