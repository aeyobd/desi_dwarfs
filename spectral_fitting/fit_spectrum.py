import numpy as np
from pathlib import Path
from copy import deepcopy

import tomli_w
from desiutil.dust import dust_transmission
import fitsio

import matplotlib.pyplot as plt
import arya
arya.use("mnras")

import sys
sys.path.append("..")

import desi_utils
from desi_utils import  get_wavelength
from desi_utils import LineFitter



LINE_WAVES = np.array([desi_utils.get_wavelength(l) for l in desi_utils.EMISSION_LINES])

def tomlize(obj):
    # numpy scalars (np.float32, np.int64, np.bool_, ...) -> native Python types
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, dict):
        return {str(k): tomlize(v) for k, v in obj.items()}
    if hasattr(obj, "to_dict"):                 # pandas DataFrame / Series
        return tomlize(obj.to_dict(orient="list") if hasattr(obj, "columns") else obj.to_dict())
    if isinstance(obj, (list, tuple)) or hasattr(obj, "__array__"):
        return [tomlize(x) for x in np.asarray(obj).tolist()] if hasattr(obj, "__array__") else [tomlize(x) for x in obj]
    raise TypeError(f"cannot convert {type(obj)} to TOML")



def reduce_spectrum(galspec_raw):
    if len(galspec_raw.flux.shape) == 2:
        print("spectra has multiple elements, returning first")
        galspec_new= galspec_raw[0]
        for key in galspec_raw.meta:
            if key in ["original_wcs"]:
                continue
            galspec_new.meta[key] = galspec_raw.meta[key][0]
    else:
        galspec_new = galspec_raw

    return desi_utils.Spectrum.from_desi(galspec_new)

def correct_mw_dust(spec, EBV):
    dust_trans = dust_transmission(spec.wavelength, EBV)
    spec_corr = deepcopy(spec)
    spec_corr.flux /= dust_trans
    spec_corr.uncertainty /= dust_trans

    return spec_corr


def get_fits_row(filename, hdu_id, targetid, targetid_col="TARGETID"):
    with fitsio.FITS(filename) as f:
        hdu = f[hdu_id]
        idx = hdu.where(f"{targetid_col} == {targetid}")[0]
        row = hdu[idx]

        df = {col: val for col, val in zip(hdu.get_colnames(), row)}

    return df


def get_properties(targetid):
    dwarf_file = "../data/desi_dr1_dwarf_catalog.fits"

    props_main = get_fits_row(dwarf_file, 1, targetid)
    props_phot = get_fits_row(dwarf_file, 3, targetid)
    props_fsf = get_fits_row(dwarf_file, 4, targetid)
    props_mel = get_fits_row("../data/stellar_mass_emline_dwarfs.fits", 1, targetid)

    return {
        "main": props_main,
        "phot": props_phot,
        "fsf_dwarf": props_fsf,
        "mel": props_mel,
    }

def get_spectrum(props):
    targetid = props["main"]["TARGETID"]
    galspec_raw = desi_utils.get_data.get_spectra(targetid)
    galspec_redux = reduce_spectrum(galspec_raw, )
    EBV = props["phot"]["EBV"]

    spec = correct_mw_dust(galspec_redux, EBV)

    return spec

 
def fit_region(flat_spec, region):
    line_filt = LINE_WAVES > region[0]
    line_filt &= LINE_WAVES < region[1]
    line_names = np.array(desi_utils.EMISSION_LINES)[line_filt]
    print("fitting lines: ", line_names)
    
    subspec = desi_utils.spectrum.clip_spectrum(flat_spec, region)
    line_fitter = LineFitter(line_names=line_names, progressbar=False)
    fit = line_fitter.fit(subspec)    

    return fit


def extract_region(lines, wavelength_range, width = 100):
    if isinstance(lines, str):
        line = lines
        lambda0 = get_wavelength(line)
        if lambda0 < wavelength_range[0] or lambda0 > wavelength_range[1]:
            return None
        region = get_wavelength(line) - width / 2, get_wavelength(line) + width / 2

    else:
        lambdas = [get_wavelength(line) for line in lines]
        good_lambdas = [wavelength_range[0] < lam < wavelength_range[1] for lam in lambdas]
        lambdas = np.array(lambdas)[good_lambdas]
        
        region = np.min(lambdas) - width/2, np.max(lambdas) + width / 2

    return region


def construct_regions(spec):
    multiplits = [
        ("O2_3726A", "O2_3729A"),
        ("Ne3_3968A", "H1r_3970A"),
        ("H1r_4341A", "O3_4363A"),
        ("He2r_4686A", "Ar4_4711A", "Ar4_4740A"),
        ("O3_4959A", "O3_5007A"),
        ("O1_6300A", "S3_6312A", "O1_6364A"), 
        ('N2_6548A', 'H1r_6563A', 'N2_6584A'),
        ("S2_6716A", "S2_6731A"),
        ("He1r_7281A", "O2_7320A", "O2_7331A"),
    ]

    all_multiplit_lines = np.concat(multiplits)
    single_lines = [line for line in desi_utils.EMISSION_LINES if line not in all_multiplit_lines]
    assert len(single_lines) + len(all_multiplit_lines) == len(desi_utils.EMISSION_LINES)

    wave_range = (spec.get_wavelength(restframe=True).min(), spec.get_wavelength(restframe=True).max())
    single_regions = [extract_region(line, wave_range) for line in single_lines]
    multi_regions = [extract_region(line, wave_range) for line in multiplits]
    all_regions = [*single_regions, *multi_regions]
    return all_regions

def fit_regions(flat_spec):
    all_regions = construct_regions(flat_spec)
    region_fits = []
    for region in all_regions:
        if region is None:
            continue
        fit = fit_region(flat_spec, region)
        region_fits.append(fit)

    return region_fits



def summarize_fits(region_fits):
    mc_fluxes = {}
    for fit in region_fits:
        for line in fit.line_names:
            df = fit.summary()["median"]
            df_err = fit.summary()["mad"]
            
            flux = fit.line_fluxes()["median"][f"line_flux[{line}]"]
            flux_err = fit.line_fluxes()["mad"][f"line_flux[{line}]"]
            
            mu = df[f"mu[{line}]"]
            mu_err = df_err[f"mu[{line}]"]
            
            
            
            sigma = df[f"sigma[{line}]"]
            sigma_err = df_err[f"sigma[{line}]"]

            mc_fluxes[line + "_flux"] = flux
            mc_fluxes[line + "_flux_err"] = flux_err
            mc_fluxes[line + "_sigma"] = sigma
            mc_fluxes[line + "_sigma_err"] = sigma_err
            mc_fluxes[line + "_cen"] = mu
            mc_fluxes[line + "_cen_err"] = mu_err

    return mc_fluxes

def to_pyframe(df):
    return {k: float(v) for k, v in df.items()}

def summarize_ppxf(ppxf_fit):
    df_ppxf = {
        "weights": ppxf_fit.weights,
        "gas_bestfit": ppxf_fit.gas_bestfit,
        "bestfit": ppxf_fit.bestfit,
        "sol": ppxf_fit.sol,
        "gas_names": ppxf_fit.gas_names,
        "chi2": ppxf_fit.chi2,
        "component": ppxf_fit.component,
        "gas_component": ppxf_fit.gas_component,
        "gas_reddening": ppxf_fit.gas_reddening,
        "reg_dim": ppxf_fit.reg_dim
    }

    return df_ppxf


def write_toml(df, filename):
    with open(filename, "wb") as f:
        tomli_w.dump(tomlize(df), f)



def plot_fit(out_dir, fit):
    for line in fit.line_names:
        fig, ax = plt.subplots(1)
        
        fit.plot_line_fit(line, wave_range=10, ax=ax, show_legend=False)
        ax.set(title = desi_utils.spectral_lines.format_line(line))
        
        ax.legend()

        plt.savefig(out_dir / f"lines/{line}.png")




def main(targetid):
    out_dir = Path(".", "example_fits", str(targetid))
    if not out_dir.exists():
        out_dir.mkdir()

    print("reading properties")
    props = get_properties(targetid)
    write_toml(props, out_dir / "properties.toml")

    print("loading spectrum")
    spec = get_spectrum(props)
    spec_resamp, flux_scale = desi_utils.spectrum.log_rebin_spectrum(spec)

    print("fitting continuum")
    sedmodel = desi_utils.fit_continuum.SEDModel(
        sps_file="spectra_fsps_9.0.npz"
    )
    ppxf_fit = sedmodel.fit(spec_resamp)

    stars_bestfit = ppxf_fit.bestfit - ppxf_fit.gas_bestfit
    flat_spec = desi_utils.Spectrum(spec_resamp.wavelength, (spec_resamp.flux - stars_bestfit) * flux_scale, 
                         uncertainty=spec_resamp.uncertainty * flux_scale, redshift=spec_resamp.redshift)

    print("fitting spectral lines")
    region_fits = fit_regions(flat_spec)
    mc_fluxes = summarize_fits(region_fits)

    print("writing results")
    write_toml(to_pyframe(mc_fluxes), out_dir / "mc_fluxes.toml")
    write_toml(summarize_ppxf(ppxf_fit), out_dir / "ppxf_fit.toml")
    
    print("plotting line fits")
    for fit in region_fits:
        plot_fit(out_dir, fit)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: ", sys.argv[0], "targetid")

    targetid = int(sys.argv[1])
    main(targetid)
