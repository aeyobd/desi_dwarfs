from __future__ import annotations
import math
from dataclasses import dataclass, field

import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
import warnings

from ppxf import ppxf_util

from astropy.convolution import convolve, Gaussian1DKernel


def inv_var_to_std(inv_var: np.ndarray):
    """
        inv_var_to_std(inv_var: np.ndarray)
    
    Given an array of inverse variances, returns the standard deviation
    """
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="divide by zero encountered in divide")
        std = np.where(inv_var == 0, np.inf, 1 / np.sqrt(inv_var))
    
    return std


def to_units(a:u.Quantity, unit):
    """
    Converts an astropy quantity to the given units and returns the 
    value only
    """
    decomposed = (a / unit).decompose()
    assert decomposed.unit == u.Unit(), "units incompatable"
    return decomposed.value


@dataclass 
class Spectrum:
    wavelength: np.ndarray
    flux: np.ndarray

    uncertainty: np.ndarray
    wave_dispersion: np.ndarray = None

    redshift: float = math.nan
    velocity_scale: float = 1
    meta: dict = field(default_factory=dict)

    flux_unit: u.Quantity = 1 
    wavelength_unit: u.Quantity  = 1

    @classmethod
    def from_desi(cls, spec):
        flux_unit = 1e-17 * u.erg / u.s / u.cm**2 / u.Angstrom
        wavelength_unit = u.Angstrom
        
        uncertainty = inv_var_to_std(to_units(spec.uncertainty.array * spec.uncertainty.unit, flux_unit**-2))
        
        return cls(
            wavelength = to_units(spec.spectral_axis, wavelength_unit),
            flux = to_units(spec.flux, flux_unit),
            uncertainty = uncertainty,
            wave_dispersion = to_units(spec.meta["wave_sigma"], wavelength_unit),
            redshift = spec.meta["redshift"],
            wavelength_unit = wavelength_unit, 
            flux_unit = flux_unit,
            meta = spec.meta
        )


    def plot(self, ax=None, restframe=False, plot_uncertainty=False, xlim=None, lw=0.5, **kwargs):
        if ax is None:
            ax = plt.gca()

        wave = self.get_wavelength(restframe=restframe)
        flux = self.flux
        err = self.uncertainty

        if xlim is not None:
            filt = np.full_like(wave, True, dtype=bool)
            if xlim[0] is not None:
                filt &= (wave >= xlim[0])
            if xlim[1] is not None:
                filt &= (wave <= xlim[1])

            wave = wave[filt]
            flux = flux[filt]
            err = err[filt]
            ax.set_xlim(xlim)


        p = ax.plot(wave, flux, lw=0.5, **kwargs)

        if plot_uncertainty:
            color = p[0].get_color()
            ax.fill_between(wave, flux - err, flux + err, 
                            color=color, alpha=0.2)


        if ax.get_xlabel() == '':
            if restframe:
                suffix = "rest"
            else:
                suffix = "obs"

            units = str(self.wavelength_unit)
            units = units.replace("Angstrom", r"\AA")
            ax.set_xlabel(f"$\\lambda_\\text{{{suffix}}}\ / \ {units}$")
        if ax.get_ylabel() == '':

            units = str(self.flux_unit)
            units = units.replace("Angstrom", r"\AA{}")
            ax.set_ylabel(f"flux / {units}")



    def get_wavelength(self, restframe=False):
        if restframe:
            wave = self.wavelength / (1 + self.redshift)
        else:
            wave = self.wavelength

        return wave

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
    flux_scale = np.median(flux_resam)

    uncertainty_resam = np.interp(lambda_resam, spec.wavelength, spec.uncertainty) / flux_scale
    flux_resam /= flux_scale

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
    ), flux_scale



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
        flux_unit = spec.flux_unit,
        wavelength_unit = spec.wavelength_unit,
        meta = spec.meta
        )


def clip_spectrum(spectrum, waverange):
    lambda_min = waverange[0] * (1 + spectrum.redshift)
    lambda_max = waverange[1] * (1 + spectrum.redshift)


    imin = np.where(spectrum.wavelength >= lambda_min)[0][0]
    imax = np.where(spectrum.wavelength <= lambda_max)[0][-1]
    idxs = np.arange(imin, imax)
    flux = spectrum.flux
    wave = spectrum.wavelength
    redshift = spectrum.redshift

    
    spec2 = Spectrum(
        flux = flux[idxs], 
        wavelength = wave[idxs], 
        uncertainty = spectrum.uncertainty[idxs], 
        wave_dispersion = None if spectrum.wave_dispersion is None else spectrum.wave_dispersion[idxs], 
        flux_unit = spectrum.flux_unit,
        wavelength_unit = spectrum.wavelength_unit,
        redshift = redshift, 
        velocity_scale = spectrum.velocity_scale, 
        meta = spectrum.meta
    )

    return spec2
