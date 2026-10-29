import os
from copy import copy
from dataclasses import dataclass, asdict, field
import math

import numpy as np

import ppxf.sps_util as ppxf_lib
from ppxf.ppxf import ppxf
from ppxf.ppxf_util import gaussian

from .spectrum import Spectrum, C_KMS
from .spectral_lines import EMISSION_LINES, get_wavelength
import tomli_w

SPS_DIR = os.path.join(os.path.dirname(__file__), "..", "spectral_fitting/sps_templates")


def redshift_to_vel(redshift: float):
    """
        redshift_to_vel(redshift: float)

    Return the los velocity in km/s given a redshift. 
    """
    vel = C_KMS * np.log1p(redshift)   # eq.(8) of Cappellari (2017)
    return vel


def vel_to_redshift(vel: float):
    """
        convert a LOS velocity to a redshift
    """
    redshift = np.expm1(vel / C_KMS)
    return redshift


def load_sps(sps_file: str, spec:Spectrum, norm_range=(5070, 5950)):
    """
        load_sps(sps_file, redshift, lam_gal, fwhm_gal, 
        norm_range)

    Load the Stellar Population file, resampling on the given lambda and fwhm of the spectra, 
    scaling by the redshift, and 
    """
    fwhm_gal_dic = fwhm_dictionary(spec)
    sps = ppxf_lib.sps_lib(os.path.join(SPS_DIR, sps_file), spec.velocity_scale, fwhm_gal_dic, norm_range=norm_range)
    return sps


def fwhm_dictionary(spec:Spectrum):
    d = {"lam": spec.wavelength, "fwhm": sigma_to_fwhm(spec.wave_dispersion)}
    return d



def sigma_to_fwhm(sigma):
    # see https://mathworld.wolfram.com/GaussianFunction.html
    return 2 * np.sqrt(2*np.log(2)) * sigma

@dataclass
class SEDModel:
    sps_file: str

    moments: list = field(default_factory=lambda: [4, 2, 2])
    ic_single: list = field(default_factory=lambda: [math.nan, 180.])
    ic: list = field(default_factory=list)
    gas_reddening: float = 0.0
    tie_balmer: bool = False
    regul_err: float = 0.01


    _templates = None
    _spec: Spectrum = None

    def get_ic(self, redshift):
        vel = redshift_to_vel(redshift)
        if len(self.ic) > 0:
            ic = self.ic
        else:
            if math.isnan(self.ic_single[0]):
                ic_single = copy(self.ic_single)
                ic_single[0] = float(vel)
            ic = [ic_single, ic_single, ic_single]

        return ic


    def _is_compatible_spectrum(self, spec):
        spec2 = self._spec
        return (np.all(np.isclose(spec.wavelength, spec2.wavelength)) 
                and np.isclose(spec.redshift, spec2.redshift)
                and np.all(np.isclose(spec.wave_dispersion, spec2.wave_dispersion))
            )


    def load_templates(self, spectrum):
        if self._spec is not None and self._is_compatible_spectrum(spectrum):
            return
        sps = load_sps(self.sps_file, spectrum)
        stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)
        redshift = spectrum.redshift

        lam_range_gal = np.array([np.min(spectrum.wavelength), np.max(spectrum.wavelength)])/(1 + redshift)
        fwhm_gal_dic = fwhm_dictionary(spectrum)

        gas_templates, gas_names, line_wave = emission_lines(
            sps.ln_lam_temp, lam_range_gal, fwhm_gal_dic)

        self._templates = np.column_stack([stars_templates, gas_templates])
        self._stars_templates = stars_templates
        self._gas_names = gas_names
        self._line_wave = line_wave
        self._sps = sps
        self._spec = spectrum


    @property
    def components(self):
        n_temps = self._stars_templates.shape[1]
        n_balmer = np.sum(["r_" in a for a in self._gas_names])
        n_forbidden = len(self._gas_names) - n_balmer

        component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
        return component

    @property
    def gas_components(self):
        return np.array(self.components) > 0

    @property
    def uncertainty(self):
        uncertainty = copy(self._spec.uncertainty)
        uncertainty[~np.isfinite(uncertainty)] = 1e10 * np.max(uncertainty[np.isfinite(uncertainty)])
        uncertainty[uncertainty == 0] =  1e10 * np.max(uncertainty[np.isfinite(uncertainty)])
        return uncertainty


    @property
    def reg_dim(self):
        return self._sps.templates.shape[1:]


    def fit(self, spectrum):
        ic = self.get_ic(spectrum.redshift)
        self.load_templates(spectrum)
        sps = self._sps


        kwargs = {
            "moments": self.moments,
            "degree": -1,
            "mdegree": 10,
            "lam": spectrum.wavelength,
            "lam_temp": sps.lam_temp,
            "component": self.components,
            "gas_component": self.gas_components,
            "gas_names": self._gas_names,
            "gas_reddening": self.gas_reddening,
        }

        if self.regul_err > 0:
            kwargs["regul"] = 1/self.regul_err
            kwargs["reg_dim"] = self.reg_dim

        p = ppxf(self._templates, spectrum.flux, self.uncertainty, 
                 spectrum.velocity_scale, ic, **kwargs)

        self.p_fit = p

        return p



    def __str__(self):
        d = {
            "SEDModel": asdict(self)
        }
        return tomli_w.dumps(d)

    def __repr__(self):
        return str(self)





def emission_lines(ln_lam_temp, lam_range_gal, FWHM_gal, pixel=True, vacuum=True):
    """Generates Gaussian emission line templates for pPXF.

    This function creates an array of Gaussian templates for gas emission
    lines, intended for use with pPXF.

    The templates typically represent the instrumental Line Spread Function
    (LSF) at the wavelength of each emission line. When using these templates,
    pPXF fits for the intrinsic (astrophysical) velocity dispersion of the gas.
    Alternatively, if `FWHM_gal=0` is provided, the lines are treated as
    delta functions, and pPXF returns the observed dispersion, which is a
    combination of the intrinsic and instrumental broadening.

    The function includes options to handle common physical constraints:
    - The [OI], [OIII], and [NII] doublets are fixed at theoretical flux ratios.
    - The [OII] and [SII] doublets can be constrained to physically plausible
      flux ratios via the `limit_doublets` parameter.
    - The Balmer series can be fixed to a theoretical decrement via the
      `tie_balmer` parameter.

    Parameters
    ----------
    ln_lam_temp : array_like
        Natural logarithm of the template wavelengths in Angstroms. This
        should match the wavelength grid of the stellar templates.
    lam_range_gal : array_like
        A 2-element array specifying the estimated rest-frame wavelength
        range of the galaxy spectrum. Typically calculated as::

            lam_range_gal = np.array([np.min(wave), np.max(wave)]) / (1 + z)
    FWHM_gal : float, callable, or dict
        Instrumental resolution (FWHM) in Angstroms. This can be:
        - A scalar value for constant FWHM.
        - A function `f(wavelength)` that returns the FWHM for given wavelengths.
        - A dictionary `{'lam': lam, 'fwhm': fwhm}` specifying the FWHM at
          each pixel of the galaxy spectrum.
    pixel : bool, default: True
        If True, analytically integrates the Gaussian LSF over each pixel for
        higher accuracy.
    vacuum : bool, default: False
        If True, assumes the provided line wavelengths are in vacuum. By
        default, they are assumed to be in air.

    Returns
    -------
    emission_lines : numpy.ndarray
        An array of shape `(ln_lam_temp.size, n_lines)` containing the
        gas emission line templates.
    line_names : numpy.ndarray
        An array of strings with the name for each gas template.
    line_wave : numpy.ndarray
        An array of the central wavelength for each gas template.

    """

    if isinstance(FWHM_gal, dict):
        FWHM_gal1 = lambda lam: np.interp(lam, FWHM_gal["lam"], FWHM_gal["fwhm"])
    else:
        FWHM_gal1 = FWHM_gal

    line_names = np.array(EMISSION_LINES)
    line_wave = np.array([get_wavelength(line, vacuum=vacuum) for line in
        line_names])
    emission_lines = gaussian(ln_lam_temp, line_wave, FWHM_gal1, pixel)


    # Only include lines falling within the estimated fitted wavelength range.
    w = (lam_range_gal[0] < line_wave) & (line_wave < lam_range_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave
