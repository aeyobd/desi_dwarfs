from __future__ import annotations
import math
from dataclasses import dataclass, field

import numpy as np
from matplotlib import pyplot as plt
from astropy.units import Quantity

from ppxf import ppxf_util

from astropy.convolution import convolve, Gaussian1DKernel

@np.vectorize
def inv_var_to_std(inv_var):
    if inv_var == 0:
        return np.inf
    else:
        return 1 / np.sqrt(inv_var)


@dataclass 
class Spectrum:
    wavelength: numpy.ndarray
    flux: numpy.ndarray

    uncertainty: numpy.ndarray
    wave_dispersion: numpy.ndarray = None
    redshift: float = math.nan
    velocity_scale: float = 1
    meta: dict = field(default_factory=dict)

    flux_unit: u.Quantity = 1 
    wavelength_unit: u.Quantity  = 1

    @classmethod
    def from_desi(cls, spec):
        uncertainty = inv_var_to_std(spec.uncertainty.array)

        return cls(
            wavelength = spec.spectral_axis.value,
            flux = spec.flux.value,
            uncertainty = uncertainty,
            wave_dispersion = spec.meta["wave_sigma"],
            redshift = spec.meta["redshift"],
            wavelength_unit = spec.spectral_axis.unit,
            flux_unit = spec.flux.unit,
            meta = spec.meta
        )


    def plot(self, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()

        ax.plot(self.wavelength, self.flux, **kwargs)

    def __add__(self, other):
        assert np.all(self.wavelength == other.wavelength), "error: Subtraction only works on spectra with same axis"
        return Spectrum(
            wavelength = self.wavelength,
            flux = self.flux + other.flux,
            uncertainty = np.sqrt(self.uncertainty**2 + other.uncertainty**2),
            wave_dispersion = self.wave_dispersion,
            redshift = self.redshift, 
            velocity_scale = self.velocity_scale,
            meta = self.meta,
            flux_unit = self.flux_unit,
            wavelength_unit = self.wavelength_unit,
        )


    def __sub__(self, other):
        assert np.all(self.wavelength == other.wavelength), "error: Subtraction only works on spectra with same axis"
        return self + -1*other


    def __mul__(self, a: float):
        return Spectrum(
            wavelength = self.wavelength,
            flux = self.flux * a,
            uncertainty = self.uncertainty * np.abs(a),
            redshift = self.redshift, 
            velocity_scale = self.velocity_scale,
            wave_dispersion = self.wave_dispersion,
            meta = self.meta
        )


    def __div__(self, a: float):
        return self * (1/a)




def log_rebin_spectrum(spec: Spectrum, flux=False):
    flux_resam, log_lambda_resam, velscale_resam = ppxf_util.log_rebin(spec.wavelength, spec.flux, flux=flux)
    lambda_resam = np.exp(log_lambda_resam)
    wave_dispersion_resam = np.interp(lambda_resam, spec.wavelength, spec.wave_dispersion)

    uncertainty_resam = np.interp(lambda_resam, spec.wavelength, spec.uncertainty / np.median(flux_resam))
    flux_resam /= np.median(flux_resam)

    return Spectrum(
        wavelength = lambda_resam,
        flux = flux_resam,
        wave_dispersion =  wave_dispersion_resam,
        uncertainty = uncertainty_resam,
        velocity_scale = velscale_resam,
        redshift = spec.redshift,
        meta = spec.meta,
        flux_unit = spec.flux_unit, 
        wavelength_unit = spec.wavelength_unit,
    )



def conv_kernel_std(kern):
    """compute the standard deviation of a convolution kernel from astropy.convolution"""
    values = np.arange(len(kern.array))
    weights = kern.array

    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return np.sqrt(variance)


def smooth_spectrum(spec, kernel = Gaussian1DKernel(5)):
    flux_new = convolve(spec.flux, kernel)

    return Spectrum(
        wavelength=spec.wavelength,
        flux = convolve(spec.flux, kernel),
        uncertainty = np.sqrt(convolve(spec.uncertainty**2, kernel)),
        wave_dispersion = np.sqrt(spec.wave_dispersion**2 + conv_kernel_std(kernel)**2), # TODO: this is probably wrong
        redshift = spec.redshift,
        velocity_scale = spec.velocity_scale,
        meta = spec.meta
    )
