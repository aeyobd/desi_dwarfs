import os
from copy import copy
from dataclasses import dataclass, asdict, field
import math

import numpy as np
from astropy import constants

import ppxf.ppxf_util as ppxf_utils
import ppxf.sps_util as ppxf_lib
from ppxf.ppxf import ppxf

from .spectrum import Spectrum
import tomli_w


SPS_DIR = os.path.join(os.path.dirname(__file__), "..", "sps_templates")




def redshift_to_vel(redshift: float):
    """
        redshift_to_vel(redshift: float)

    Return the los velocity in km/s given a redshift. 
    """
    c = constants.c.to("km/s").value
    vel = c*np.log(1 + redshift)   # eq.(8) of Cappellari (2017)
    return vel


def load_sps(sps_file: str, spec:Spectrum, norm_range=[5070, 5950]):
    """
        load_sps(sps_file, redshift, lam_gal, fwhm_gal, 
        norm_range)

    Load the Stellar Population file, resampling on the given lambda and fwhm of the spectra, 
    scaling by the redshift, and 
    """
    fwhm_gal_dic = {"lam": spec.wavelength, "fwhm": spec.wave_dispersion}
    sps = ppxf_lib.sps_lib(os.path.join(SPS_DIR, sps_file), spec.velocity_scale, fwhm_gal_dic, norm_range=norm_range)
    return sps


@dataclass
class SEDModel:
    sps_file: str

    moments: list = field(default_factory=lambda: [4, 2, 2])
    ic_single: list = field(default_factory=lambda: [math.nan, 180.])
    ic: list = field(default_factory=list)
    gas_reddening: float = 0.0
    tie_balmer: bool = False
    regul_err: float = 0.01


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

    def fit(self, spectrum):
        sps = load_sps(self.sps_file, spectrum)
        reg_dim = sps.templates.shape[1:]
        stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)
        redshift = spectrum.redshift

        lam_range_gal = np.array([np.min(spectrum.wavelength), np.max(spectrum.wavelength)])/(1 + redshift)
        fwhm_gal_dic = {"lam": spectrum.wavelength, "fwhm": spectrum.wave_dispersion}

        gas_templates, gas_names, line_wave = ppxf_utils.emission_lines(
            sps.ln_lam_temp, lam_range_gal, fwhm_gal_dic, tie_balmer=self.tie_balmer,
            limit_doublets=False)

        templates = np.column_stack([stars_templates, gas_templates])

        n_temps = stars_templates.shape[1]
        n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
        n_balmer = len(gas_names) - n_forbidden

        component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
        gas_component = np.array(component) > 0  # gas_component=True for gas templates

        ic = self.get_ic(spectrum.redshift)

        uncertainty = copy(spectrum.uncertainty.value)
        uncertainty[~np.isfinite(uncertainty)] = 1e10 * np.max(uncertainty[np.isfinite(uncertainty)])
        uncertainty[uncertainty == 0] =  1e10 * np.max(uncertainty[np.isfinite(uncertainty)])

        p = ppxf(templates, spectrum.flux, uncertainty, spectrum.velocity_scale, ic, 
            moments = self.moments,
            degree = -1, mdegree = 10, 
            lam = spectrum.wavelength, 
            lam_temp =  sps.lam_temp,
            regul = 1/self.regul_err, 
            reg_dim = reg_dim, 
            component = component,
            gas_component = gas_component, 
            gas_names = gas_names,
            gas_reddening = self.gas_reddening
        )

        self.gas_component = gas_component
        self.n_temps = n_temps
        self.reg_dim = reg_dim
        self.n_forbidden = n_forbidden
        self.n_balmer = n_balmer

        self.p_fit = p

        return p


    def __str__(self):
        d = {
            "SEDModel": asdict(self)
        }
        return tomli_w.dumps(d)

    def __repr__(self):
        return str(self)




def fit_local_continuum(spec, windows):
    pass


def emission_lines(ln_lam_temp, lam_range_gal, FWHM_gal, pixel=True,
                   tie_balmer=False, limit_doublets=False, vacuum=False):
    """Generates Gaussian emission line templates for pPXF.

    This function creates an array of Gaussian templates for gas emission
    lines, intended for use with pPXF.

    .. note::
        This routine is a template. Users are welcome to copy, modify, and
        distribute it to accommodate their specific needs for different
        emission lines.

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
    tie_balmer : bool, default: False
        If True, ties the Balmer lines to a theoretical decrement (Case B
        recombination, T=1e4 K, n=100 cm^-3).

        .. important::
            This option assumes the input spectrum has flux units
            proportional to `erg/(cm**2 s A)`.
    limit_doublets : bool, default: False
        If True, constrains the [OII] and [SII] doublet ratios to physically
        allowed ranges. This is done by modeling each doublet as a linear
        combination of two templates representing the minimum and maximum
        allowed ratios. An alternative is to use the `constr_templ` keyword
        in pPXF.

        .. important::
            When using this keyword, the two output fluxes (`flux_1`
            and `flux_2`) for a doublet do not represent the actual fluxes of the
            two lines, but the weights of the two doublet templates. If the two
            templates have line ratios `rat_1` and `rat_2`, the actual fitted
            ratio and total flux are::

                flux_total = flux_1 + flux_2
                ratio_fit = (rat_1*flux_1 + rat_2*flux_2)/flux_total

            EXAMPLE: For the [SII] doublet, the adopted ratios for the templates are::

                ratio_d1 = flux([SII]6716/6731) = 0.44
                ratio_d2 = flux([SII]6716/6731) = 1.43.

            When pPXF prints (and returns in pp.gas_flux)::

                flux([SII]6731_d1) = flux_1
                flux([SII]6731_d2) = flux_2

            the total flux and true lines ratio of the [SII] doublet are::

                flux_total = flux_1 + flux_2
                ratio_fit([SII]6716/6731) = (0.44*flux_1 + 1.43*flux_2)/flux_total

            Similarly, for [OII], the adopted ratios for the templates are::

                ratio_d1 = flux([OII]3729/3726) = 0.28
                ratio_d2 = flux([OII]3729/3726) = 1.47.

            When pPXF prints (and returns in pp.gas_flux)::

                flux([OII]3726_d1) = flux_1
                flux([OII]3726_d2) = flux_2

            the total flux and true lines ratio of the [OII] doublet are::

                flux_total = flux_1 + flux_2
                ratio_fit([OII]3729/3726) = (0.28*flux_1 + 1.47*flux_2)/flux_total

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

    #        Balmer:     H10       H9         H8        Heps    Hdelta    Hgamma    Hbeta     Halpha
    balmer = np.array([3798.983, 3836.479, 3890.158, 3971.202, 4102.899, 4341.691, 4862.691, 6564.632])  # vacuum wavelengths

    if tie_balmer:

        # Balmer decrement for Case B recombination (T=1e4 K, ne=100 cm^-3)
        # from Storey & Hummer (1995) https://ui.adsabs.harvard.edu/abs/1995MNRAS.272...41S
        # In electronic form https://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/64
        # See Table B.7 of Dopita & Sutherland (2003) https://www.amazon.com/dp/3540433627
        # Also see Table 4.2 of Osterbrock & Ferland (2006) https://www.amazon.co.uk/dp/1891389343/
        wave = balmer
        if not vacuum:
            wave = vac_to_air(wave)
        gauss = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel)
        ratios = np.array([0.0530, 0.0731, 0.105, 0.159, 0.259, 0.468, 1, 2.86])
        # Account for varying log-sampled pixel size in Angstrom
        ratios *= wave[-2]/wave
        emission_lines = gauss @ ratios
        line_names = ['Balmer']
        w = (lam_range_gal[0] < wave) & (wave < lam_range_gal[1])
        line_wave = np.mean(wave[w]) if np.any(w) else np.mean(wave)

    else:

        line_wave = balmer
        if not vacuum:
            line_wave = vac_to_air(line_wave)
        line_names = ['H10', 'H9', 'H8', 'Heps', 'Hdelta', 'Hgamma', 'Hbeta', 'Halpha']
        emission_lines = gaussian(ln_lam_temp, line_wave, FWHM_gal1, pixel)
    if limit_doublets:

        # The line ratio of this doublet lam3727/lam3729 is constrained by
        # atomic physics to lie in the range 0.28--1.47 (e.g. fig.5.8 of
        # Osterbrock & Ferland (2006) https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #       -----[OII]-----
        wave = [3727.092, 3729.875]    # vacuum wavelengths
        if not vacuum:
            wave = vac_to_air(wave)
        names = ['[OII]3726_d1', '[OII]3726_d2']
        gauss = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel)
        doublets = gauss @ [[1, 1], [0.28, 1.47]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

        # The line ratio of this doublet lam6717/lam6731 is constrained by
        # atomic physics to lie in the range 0.44--1.43 (e.g. fig.5.8 of
        # Osterbrock & Ferland (2006) https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #        -----[SII]-----
        wave = [6718.294, 6732.674]    # vacuum wavelengths
        if not vacuum:
            wave = vac_to_air(wave)
        names = ['[SII]6731_d1', '[SII]6731_d2']
        gauss = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel)
        doublets = gauss @ [[0.44, 1.43], [1, 1]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    else:

        # Here the two doublets are free to have any ratio
        #         -----[OII]-----     -----[SII]-----
        wave = [3727.092, 3729.875, 6718.294, 6732.674]  # vacuum wavelengths
        if not vacuum:
            wave = vac_to_air(wave)
        names = ['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731']
        gauss = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel)
        emission_lines = np.column_stack([emission_lines, gauss])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    # Here the lines are free to have any ratio
    #       -----[NeIII]-----    HeII      HeI
    wave = [3968.59, 3869.86, 4687.015, 5877.243]  # vacuum wavelengths
    if not vacuum:
        wave = vac_to_air(wave)
    names = ['[NeIII]3968', '[NeIII]3869', 'HeII4687', 'HeI5876']
    gauss = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, wave)


    ######### Doublets with fixed ratios #########

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OIII]-----
    wave = [4960.295, 5008.240]    # vacuum wavelengths
    if not vacuum:
        wave = vac_to_air(wave)
    doublet = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])
    # single template for this doublet
    line_names = np.append(line_names, '[OIII]5007_d')
    line_wave = np.append(line_wave, wave[1])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OI]-----
    wave = [6302.040, 6365.535]    # vacuum wavelengths
    if not vacuum:
        wave = vac_to_air(wave)
    doublet = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel) @ [1, 0.33]
    emission_lines = np.column_stack([emission_lines, doublet])
    # single template for this doublet
    line_names = np.append(line_names, '[OI]6300_d')
    line_wave = np.append(line_wave, wave[0])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #       -----[NII]-----
    wave = [6549.860, 6585.271]    # air wavelengths
    if not vacuum:
        wave = vac_to_air(wave)
    doublet = gaussian(ln_lam_temp, wave, FWHM_gal1, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])

    # single template for this doublet
    line_names = np.append(line_names, '[NII]6583_d')
    line_wave = np.append(line_wave, wave[1])

    # Only include lines falling within the estimated fitted wavelength range.
    #
    w = (lam_range_gal[0] < line_wave) & (line_wave < lam_range_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave