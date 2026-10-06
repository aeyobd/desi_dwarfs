"""Bayesian fitting of (possibly blended) Gaussian emission lines with PyMC.

Typical use
-----------
>>> fitter = LineFitter(line_names=["H1r_4683A", "O3_5007A"])
>>> result = fitter.fit(spec)
>>> result.summary()
>>> result.line_fluxes()
>>> result.plot_fit(show_components=True)

`LineFitter` holds only configuration and is immutable, so one instance can
be reused on many spectra. All fit products (trace, data, plotting, derived
quantities) live on the returned `LineFitResult`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Sequence

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm
import pytensor.tensor as pt
from scipy.special import erf as np_erf

from .spectral_lines import get_wavelength

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from .spectrum import Spectrum


_SQRT2 = np.sqrt(2.0)
_SQRT_HALF_PI = np.sqrt(np.pi / 2.0)
_SQRT_2PI = np.sqrt(2.0 * np.pi)

_ARVIZ_SUMMARY_KWARGS = dict(kind = "all_median", ci_prob=0.68)

# --------------------------------------------------------------------------
# Pure helpers (shared by the PyMC model and the numpy post-processing so the
# two can never drift apart).
# --------------------------------------------------------------------------
def _pixel_edges(wave: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Lower/upper pixel edges: midpoints between samples, ends extrapolated.

    Works for non-uniform grids (e.g. log-lambda) as long as `wave` is
    strictly ascending.
    """
    if wave.size < 2:
        raise ValueError("Need at least two wavelength samples.")
    if np.any(np.diff(wave) <= 0):
        raise ValueError("`wave` must be strictly ascending.")
    mid = 0.5 * (wave[:-1] + wave[1:])
    edges = np.concatenate(([2 * wave[0] - mid[0]], mid, [2 * wave[-1] - mid[-1]]))
    return edges[:-1], edges[1:]


def _gaussian_flux(
    x,
    amp,
    mu,
    sigma,
    *,
    edges: tuple | None,
    exp: Callable,
    erf: Callable,
):
    """Gaussian flux density, all arguments pre-broadcast by the caller.

    If `edges` is None the Gaussian is evaluated at `x`. Otherwise it is
    averaged analytically over each pixel [lo, hi] using the error function,
    which conserves line flux even when a line spans only a few pixels.
    """
    if edges is None:
        return amp * exp(-0.5 * ((x - mu) / sigma) ** 2)
    lo, hi = edges
    s2 = sigma * _SQRT2
    integral = amp * sigma * _SQRT_HALF_PI * (erf((hi - mu) / s2) - erf((lo - mu) / s2))
    return integral / (hi - lo)


def _per_line(value, n: int, name: str) -> np.ndarray:
    """Broadcast a scalar or length-n sequence to a float array of length n."""
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr))
    if arr.shape != (n,):
        raise ValueError(f"`{name}` must be a scalar or have length {n}, got shape {arr.shape}.")
    return arr


# --------------------------------------------------------------------------
# LineFitter
# --------------------------------------------------------------------------
@dataclass(frozen=True, kw_only=True)
class LineFitter:
    """Fit `len(line_names)` Gaussian emission lines simultaneously.

    Parameters
    ----------
    line_names : sequence of str
        Lines to fit; rest wavelengths come from `get_wavelength(name)` and
        are shifted by the spectrum's redshift.
    amp_guess : float or sequence of float, optional
        Initial amplitudes (default: max flux within 5 pixels of each line).
    sigma_guess : float or sequence of float, optional
        Initial line widths (default: 3x median pixel spacing).
    center_prior_width : float
        Std dev (wavelength units) of the Normal prior on each line centre.
        Keep it tight enough that blended lines cannot swap identities.
    sigma_bounds : (low, high)
        Line widths get a HalfNormal prior truncated to this range.
    fit_continuum : bool
        Fit a linear residual continuum c0 + c1 * (wave - wave0).
    exact_gaussian : bool
        Average the Gaussian analytically over each pixel (erf) instead of
        evaluating it at pixel centres. Requires ascending wavelengths.
    draws, tune, chains, cores, target_accept, random_seed, progressbar :
        PyMC sampling controls.
    """

    line_names: Sequence[str]
    amp_guess: Sequence[float] | float | None = None
    sigma_guess: Sequence[float] | float | None = None
    center_prior_width: float = 3.0
    sigma_bounds: tuple[float, float] = (0.3, 15.0)
    fit_continuum: bool = True
    exact_gaussian: bool = True

    draws: int = 2000
    tune: int = 2000
    chains: int = 4
    cores: int = 1
    target_accept: float = 0.9
    random_seed: int = 42
    progressbar: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "line_names", tuple(self.line_names))
        if not self.line_names:
            raise ValueError("`line_names` must not be empty.")
        lo, hi = self.sigma_bounds
        if not 0 < lo < hi:
            raise ValueError("`sigma_bounds` must satisfy 0 < low < high.")
        if self.center_prior_width <= 0:
            raise ValueError("`center_prior_width` must be positive.")

    # ---- public API -------------------------------------------------------
    def fit(self, spec: Spectrum) -> LineFitResult:
        """Fit a `Spectrum` (uses its flux, uncertainty and redshift)."""
        return self.fit_arrays(
            spec.wavelength, spec.flux, spec.uncertainty, redshift=spec.redshift
        )

    def fit_arrays(
        self,
        wave,
        flux,
        flux_err=None,
        *,
        redshift: float = 0.0,
    ) -> LineFitResult:
        """Fit raw arrays.

        `flux_err` may be an array, a scalar, or None (a single noise scale
        is then inferred from the data).
        """
        wave, flux, flux_err = self._clean_data(wave, flux, flux_err)
        n_lines = len(self.line_names)
        centers = np.array([get_wavelength(n) for n in self.line_names]) * (1.0 + redshift)
        amp0, sigma0 = self._initial_guesses(wave, flux, centers)
        wave0 = float(wave.mean())

        model = self._build_model(wave, flux, flux_err, centers, amp0, sigma0, wave0)
        with model:
            trace = pm.sample(
                draws=self.draws,
                tune=self.tune,
                chains=self.chains,
                cores=self.cores,
                target_accept=self.target_accept,
                random_seed=self.random_seed,
                progressbar=self.progressbar,
            )

        return LineFitResult(
            trace=trace,
            model=model,
            wave=wave,
            flux=flux,
            flux_err=flux_err,
            line_names=self.line_names,
            line_centers=tuple(float(c) for c in centers),
            redshift=float(redshift),
            wave0=wave0,
            fit_continuum=self.fit_continuum,
            exact_gaussian=self.exact_gaussian,
        )

    # ---- internals --------------------------------------------------------
    @staticmethod
    def _clean_data(wave, flux, flux_err):
        wave = np.asarray(wave, dtype=float)
        flux = np.asarray(flux, dtype=float)
        good = np.isfinite(wave) & np.isfinite(flux)
        if flux_err is not None:
            flux_err = np.broadcast_to(np.asarray(flux_err, dtype=float), flux.shape)
            good &= np.isfinite(flux_err) & (flux_err > 0)
        wave, flux = wave[good], flux[good]
        flux_err = None if flux_err is None else flux_err[good]

        order = np.argsort(wave)
        wave, flux = wave[order], flux[order]
        flux_err = None if flux_err is None else flux_err[order]
        if wave.size < 2:
            raise ValueError("Fewer than two usable pixels after removing non-finite values.")
        return wave, flux, flux_err

    def _initial_guesses(self, wave, flux, centers):
        n = len(self.line_names)
        dw = np.median(np.diff(wave))

        if self.sigma_guess is None:
            sigma0 = np.full(n, 3.0 * dw)
        else:
            sigma0 = _per_line(self.sigma_guess, n, "sigma_guess")

        if self.amp_guess is None:
            amp0 = np.array(
                [
                    flux[m].max() if (m := np.abs(wave - c) < 5 * dw).any() else flux.max()
                    for c in centers
                ]
            )
        else:
            amp0 = _per_line(self.amp_guess, n, "amp_guess")
        return np.clip(amp0, 1e-6, None), sigma0

    def _build_model(self, wave, flux, flux_err, centers, amp0, sigma0, wave0) -> pm.Model:
        sigma_lo, sigma_hi = self.sigma_bounds
        flux_std = float(np.nanstd(flux)) or 1.0
        wave_span = float(np.ptp(wave))
        edges = None
        if self.exact_gaussian:
            edges = tuple(e[:, None] for e in _pixel_edges(wave))

        with pm.Model(coords={"line": list(self.line_names)}) as model:
            amp = pm.Normal("amp", mu=0.0, sigma=3.0 * amp0, dims="line", initval=amp0)
            mu = pm.Normal(
                "mu", mu=centers, sigma=self.center_prior_width, dims="line", initval=centers
            )
            # HalfNormal truncated to the allowed range. (Clipping instead
            # would give the sampler zero gradient outside the bounds.)
            sigma = pm.Truncated(
                "sigma",
                pm.HalfNormal.dist(sigma=sigma0),
                lower=sigma_lo,
                upper=sigma_hi,
                dims="line",
                initval=np.clip(sigma0, sigma_lo, sigma_hi),
            )

            profiles = _gaussian_flux(
                wave[:, None], amp[None, :], mu[None, :], sigma[None, :],
                edges=edges, exp=pt.exp, erf=pt.erf,
            )
            model_flux = profiles.sum(axis=1)

            if self.fit_continuum:
                c0 = pm.Normal("c0", mu=0.0, sigma=flux_std)
                c1 = pm.Normal("c1", mu=0.0, sigma=flux_std / wave_span)
                model_flux = model_flux + c0 + c1 * (wave - wave0)

            pm.Deterministic("line_flux", amp * sigma * _SQRT_2PI, dims="line")

            noise = flux_err if flux_err is not None else pm.HalfNormal("flux_scale", sigma=flux_std)
            pm.Normal("obs", mu=model_flux, sigma=noise, observed=flux)
        return model


# --------------------------------------------------------------------------
# LineFitResult
# --------------------------------------------------------------------------
@dataclass(frozen=True, eq=False)
class LineFitResult:
    """Posterior and data from a `LineFitter` run, with plotting helpers.

    All wavelengths stored here are in the observed frame. `restframe=True`
    on plotting methods divides by (1 + redshift) for display only.
    """

    trace: az.InferenceData
    model: pm.Model
    wave: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray | None
    line_names: tuple[str, ...]
    line_centers: tuple[float, ...]  # observed frame
    redshift: float
    wave0: float  # reference wavelength of the linear continuum
    fit_continuum: bool
    exact_gaussian: bool

    # ---- posterior access ------------------------------------------------
    def _posterior(self):
        post = self.trace.posterior
        return post.to_dataset() if hasattr(post, "to_dataset") else post

    def _samples(self, max_samples: int | None) -> dict[str, np.ndarray]:
        """Posterior draws as arrays: line params (L, S), continuum (S,)."""
        post = self._posterior().stack(sample=("chain", "draw"))
        n = post.sizes["sample"]
        if max_samples is not None and n > max_samples:
            post = post.isel(sample=np.linspace(0, n - 1, max_samples).astype(int))
        out = {k: post[k].transpose("line", "sample").values for k in ("amp", "mu", "sigma")}
        if self.fit_continuum:
            out["c0"] = post["c0"].values
            out["c1"] = post["c1"].values
        return out

    def _means(self) -> dict[str, np.ndarray]:
        """Posterior means in the same layout as `_samples` with S = 1."""
        post = self._posterior().mean(dim=("chain", "draw"))
        out = {k: post[k].values[:, None] for k in ("amp", "mu", "sigma")}
        if self.fit_continuum:
            out["c0"] = np.atleast_1d(post["c0"].values)
            out["c1"] = np.atleast_1d(post["c1"].values)
        return out

    # ---- model evaluation ------------------------------------------------
    def _evaluate(
        self, wave: np.ndarray, params: dict[str, np.ndarray]
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return per-line profiles (n_wave, L, S) and continuum (n_wave, S)."""
        edges = None
        if self.exact_gaussian:
            edges = tuple(e[:, None, None] for e in _pixel_edges(wave))
        comps = _gaussian_flux(
            wave[:, None, None],
            params["amp"][None],
            params["mu"][None],
            params["sigma"][None],
            edges=edges,
            exp=np.exp,
            erf=np_erf,
        )
        n_samples = comps.shape[2]
        if self.fit_continuum:
            cont = params["c0"][None, :] + params["c1"][None, :] * (wave[:, None] - self.wave0)
        else:
            cont = np.zeros((wave.size, n_samples))
        return comps, cont

    def best_fit_curve(
        self,
        wave: np.ndarray | None = None,
        interval: float | None = 0.94,
        max_samples: int | None = 1000,
    ):
        """Posterior-mean model and an equal-tailed credible interval.

        Parameters
        ----------
        wave : observed-frame wavelengths (ascending); defaults to the fit grid.
        interval : credible mass (e.g. 0.94), or None to skip the interval.
        max_samples : evenly thin the posterior to this many draws to bound
            memory (n_wave x n_lines x n_draws floats). None uses all draws.

        Returns
        -------
        mean, lo, hi : arrays (lo/hi are None if `interval` is None)
        """
        wave = self.wave if wave is None else np.asarray(wave, dtype=float)
        comps, cont = self._evaluate(wave, self._samples(max_samples))
        model = comps.sum(axis=1) + cont  # (n_wave, S)
        mean = model.mean(axis=1)
        if interval is None:
            return mean, None, None
        tail = 50.0 * (1.0 - interval)
        lo, hi = np.percentile(model, [tail, 100.0 - tail], axis=1)
        return mean, lo, hi

    # ---- tabular summaries -----------------------------------------------
    def summary(self, var_names: Sequence[str] | None = None):
        """ArviZ summary table of the posterior."""
        if var_names is None:
            var_names = ["amp", "mu", "sigma", "line_flux"]
            if self.fit_continuum:
                var_names += ["c0", "c1"]
        return az.summary(self.trace, var_names=list(var_names), **_ARVIZ_SUMMARY_KWARGS)

    def line_fluxes(self):
        """Integrated flux per line, amp * sigma * sqrt(2 pi), from the posterior."""
        return az.summary(self.trace, var_names=["line_flux"],
                **_ARVIZ_SUMMARY_KWARGS)

    # ---- plotting --------------------------------------------------------
    def _x(self, wave: np.ndarray, restframe: bool) -> np.ndarray:
        # Assumes rest = observed / (1 + z), the usual convention.
        return wave / (1.0 + self.redshift) if restframe else wave

    @staticmethod
    def _get_ax(ax: Axes | None) -> Axes:
        return plt.subplots(figsize=(9, 5))[1] if ax is None else ax

    def _plot_data_and_fit(self, ax, mask, restframe, components, interval,
            show_legend=True):
        x = self._x(self.wave, restframe)
        ax.errorbar(
            x[mask], self.flux[mask],
            yerr=None if self.flux_err is None else self.flux_err[mask],
            fmt=".", color="k", ms=12, label="data", zorder=3,
        )
        mean, lo, hi = self.best_fit_curve(interval=interval)
        ax.plot(x[mask], mean[mask], color="crimson", lw=2, label="total fit", zorder=2)
        if lo is not None:
            ax.fill_between(x[mask], lo[mask], hi[mask], color="crimson", alpha=0.2, zorder=1)
        for name, comp in components:
            ax.plot(x[mask], comp[mask], "--", lw=1.3, label=name, zorder=2)
        ax.set_xlabel("Rest wavelength" if restframe else "Wavelength")
        ax.set_ylabel("Flux")
        if show_legend:
            ax.legend(fontsize=8)

    def _mean_components(self, indices: Sequence[int]):
        comps, cont = self._evaluate(self.wave, self._means())
        return [(self.line_names[i], comps[:, i, 0] + cont[:, 0]) for i in indices]

    def plot_fit(
        self,
        ax: Axes | None = None,
        show_components: bool = False,
        restframe: bool = False,
        interval: float | None = 0.94,
        show_legend: bool = True,
    ) -> Axes:
        ax = self._get_ax(ax)
        components = (
            self._mean_components(range(len(self.line_names))) if show_components else []
        )
        self._plot_data_and_fit(
            ax, np.ones(self.wave.size, dtype=bool), restframe, components,
            interval, show_legend=show_legend
        )
        return ax

    def plot_line_fit(
        self,
        line: str,
        ax: Axes | None = None,
        wave_range: float = 30.0,
        restframe: bool = False,
        interval: float | None = 0.94,
        show_legend: bool = True,
    ) -> Axes:
        """Zoom on one line; `wave_range` is in the displayed frame."""
        if line not in self.line_names:
            raise ValueError(f"Unknown line {line!r}; fitted lines: {list(self.line_names)}")
        idx = self.line_names.index(line)
        ax = self._get_ax(ax)
        center = float(self._x(np.asarray(self.line_centers[idx]), restframe))
        x = self._x(self.wave, restframe)
        mask = np.abs(x - center) < wave_range
        self._plot_data_and_fit(ax, mask, restframe,
                self._mean_components([idx]), interval, show_legend =
                show_legend)
        return ax

    def plot_corner(self, var_names: Sequence[str] | None = None):
        var_names = list(var_names) if var_names is not None else ["amp", "mu", "sigma"]
        return az.plot_pair(
            self.trace, var_names=var_names, kind="kde", marginals=True, figsize=(9, 9)
        )




def fit_lines_exact(wave, flux, flux_err, line_names: Sequence[str], **kwargs) -> LineFitResult:
    """Backwards-compatible wrapper: pixel-integrated fit of raw arrays (z = 0)."""
    return LineFitter(line_names=line_names, exact_gaussian=True, **kwargs).fit_arrays(
        wave, flux, flux_err
    )
