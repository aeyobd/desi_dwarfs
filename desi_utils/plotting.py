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
        xlim = (3550, 9850), 
        N_chunks = 4,
        xlabel = "wavelength", ylabel="flux", **kwargs):

    dx = np.diff(xlim)[0] / N_chunks

    fig, axs = plt.subplots(N_chunks, 1, figsize=(8, 2*N_chunks))

    for i in range(N_chunks):
        ax = axs[i]

        for label, spec in inputs.items():
            x, y = get_spectra_xy(spec)
            style = styles.get(label, dict())
            xlim_frame = xlim[0] + i*dx, xlim[0] + (1+i) * dx
            filt = x >= xlim_frame[0]
            filt &= x < xlim_frame[1]
            ax.plot(x[filt], y[filt], label=label, **style)
        
        ax.set(
            xlim = xlim_frame,
            xlabel = xlabel, 
            ylabel = ylabel, 
            **kwargs
        )

        if i == N_chunks - 1:
            ax.legend()
        
