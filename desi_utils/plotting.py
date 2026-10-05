from matplotlib import pyplot as plt
import numpy as np
from .spectrum import Spectrum


def get_spectra_xy(spec):
    if isinstance(spec, tuple):
        x, y = spec
    elif isinstance(spec, Spectrum):
        x, y = spec.wavelength, spec.flux
    else:
        raise Error(f"Unknown spectral type: {type(spec)}")

    return x, y


def plot_wrapped_spectrum(inputs, styles=dict(), *, 
        xlim = None,
        N_chunks = 8,
        xlabel = "wavelength", ylabel="flux", **kwargs):

    if xlim is None:
        first_spectrum = next(iter(inputs.values()))
        x = get_spectra_xy(first_spectrum)[0]
        xmin = np.min(x) 
        xmax = np.max(x)
        dx_tot = xmax - xmin
        pad = 0.001
        xmin -= pad * dx_tot
        xmax += pad*dx_tot
        xlim = (xmin, xmax)


    dx = np.diff(xlim)[0] / N_chunks

    fig, axs = plt.subplots(N_chunks, 1, figsize=(16, 1*N_chunks))

    for i in range(N_chunks):
        ax = axs[i]

        for label, spec in inputs.items():
            x, y = get_spectra_xy(spec)
            style = styles.get(label, dict())
            xlim_frame = xlim[0] + i*dx, xlim[0] + (1+i) * dx
            filt = x >= xlim_frame[0]
            filt &= x < xlim_frame[1]
            ax.plot(x[filt], y[filt], label=label, **style)
        
        if i == N_chunks - 1:
            xlabel_this = xlabel
        else:
            xlabel_this = ""
        ax.set(
            xlim = xlim_frame,
            xlabel = xlabel_this, 
            ylabel = ylabel, 
            **kwargs
        )

        if i == N_chunks - 1:
            ax.legend()
        

    fig.tight_layout()
