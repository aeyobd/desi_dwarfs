# SPARCL

SPARCL is the easiest means of recovering DESI spectra. Otherwise, we need to retrieve the specific healpix file. While SPARCL returns most all the information we need, it notably excludes the resolution matrix and only returns the wavelength dispersion. I am not sure yet how much this detail matters yet. 

| ?    | Field name                   | Units                   | Description                                                  |
| ---- | ---------------------------- | ----------------------- | ------------------------------------------------------------ |
|      | `tsnr2_qso`                  | None                    | QSO template (S/N)^2 summed over B,R,Z                       |
|      | `sv2_scnd_target`            | None                    | Secondary target selection bitmask for SV2                   |
| *    | `survey`                     | None                    | Survey name                                                  |
|      | `main_nspec`                 | None                    | Number of coadded spectra for this TARGETID in Main survey   |
|      | `bgs_target`                 | None                    | BGS (Bright Galaxy Survey) target selection bitmask          |
|      | `program`                    | None                    | DESI program type - BRIGHT, DARK, BACKUP, OTHER              |
|      | `zcat_nspec`                 | None                    | Number of coadded spectra for this TARGETID in this zcatalog |
|      | `desi_target`                | None                    | DESI (dark time program) target selection bitmask            |
|      | `sv3_scnd_target`            | None                    | Secondary target selection bitmask for SV3                   |
|      | `tsnr2_gpbbackup`            | None                    | GPBBACKUP template (S/N)^2 summed over BRZ                   |
|      | `sv_nspec`                   | None                    | Number of coadded spectra for this TARGETID in SV (SV1+2+3)  |
|      | `coadd_fiberstatus`          | None                    | Bitwise AND of the FIBERSTATUS values of the spectra used in the coadd |
|      | `pmdec`                      | mas year-1              | Reference catalog proper motion in the Dec direction         |
|      | `std_fiber_ra`               | arcsec                  | Standard deviation (over exposures) of RA of actual fiber position |
|      | `sv2_mws_target`             | None                    | MWS (bright time program) target selection bitmask for SV2   |
|      | `tsnr2_gpbbright`            | None                    | GPBBRIGHT template (S/N)^2 summed over B,R,Z                 |
|      | `sv1_mws_target`             | None                    | MWS (bright time program) target selection bitmask for SV1   |
|      | `objtype`                    | None                    | Object type: TGT, SKY, NON, BAD                              |
|      | `plate_ra`                   | deg                     | Right Ascension to be used by PlateMaker                     |
|      | `coadd_numexp`               | None                    | Number of exposures included in the coadd for this target    |
|      | `sv3_desi_target`            | None                    | DESI (dark time program) target selection bitmask for SV3    |
|      | `rms_delta_y`                | mm                      | RMS (over exposures) of the fiber difference between measured and requested CS5 Y location on focal plane |
|      | `ref_epoch`                  | year                    | Reference epoch for Gaia/Tycho astrometry. Typically 2015.5 for Gaia |
|      | `rms_delta_x`                | mm                      | RMS (over exposures) of the fiber difference between measured and requested CS5 X location on focal plane |
|      | `sv3_mws_target`             | None                    | MWS (bright time program) target selection bitmask for SV3   |
|      | `sv2_bgs_target`             | None                    | BGS (bright time program) target selection bitmask for SV2   |
|      | `tsnr2_elg`                  | None                    | ELG template (S/N)^2 summed over B,R,Z                       |
|      | `mean_fiber_ra`              | deg                     | Mean (over exposures) RA of actual fiber position            |
|      | `sv1_desi_target`            | None                    | DESI (dark time program) target selection bitmask for SV1    |
|      | `mean_delta_x`               | mm                      | Mean (over exposures) fiber difference between measured and requested CS5 X location on focal plane |
|      | `healpix`                    | None                    | HEALPixel containing this location at NSIDE=64 in the NESTED scheme |
|      | `fa_target`                  | None                    | Targeting bit internally used by fiberassign (linked with FA_TYPE) |
|      | `tsnr2_bgs`                  | None                    | BGS template (S/N)^2 summed over B,R,Z                       |
|      | `mean_mjd`                   | days                    | Mean (over exposures) of the Modified Julian Date (MJD) when the shutter was opened |
|      | `tsnr2_lya`                  | None                    | LYA template (S/N)^2 summed over B,R,Z                       |
|      | `deltachi2`                  | None                    | Delta-chi-squared for template fit from Redrock              |
|      | `scnd_target`                | None                    | Target selection bitmask for secondary programs              |
|      | `coadd_numtile`              | None                    | Number of different tiles included in the coadd for this target |
|      | `coeff`                      | None                    | Redrock template coefficients                                |
|      | `coadd_numnight`             | None                    | Number of different nights included in the coadd for this target |
|      | `tsnr2_gpbdark`              | None                    | GPBDARK template (S/N)^2 summed over BRZ                     |
|      | `plate_dec`                  | deg                     | Declination to be used by PlateMaker                         |
|      | `main_primary`               | None                    | Boolean flag (True/False) for the primary coadded spectrum in Main survey |
|      | `mean_fiber_dec`             | deg                     | Mean (over exposures) DEC of actual fiber position           |
|      | `subpriority`                | None                    | Random subpriority [0-1] to break assignment ties            |
|      | `mean_psf_to_fiber_specflux` | None                    | Mean of input exposures fraction of light from point-like source captured by 1.5 arcsec diameter fiber given atmospheric seeing |
|      | `desiname`                   | None                    | DESI JXXX.XXXX[+/-]YY.YYYY where X, Y are TARGET_{RA, DEC} truncated to 4 decimals |
|      | `spgrpval`                   | None                    | Value corresponding to the grouping type indicated in the spgrp column. This is the same as the healpix column |
|      | `npixels`                    | None                    | Number of unmasked wavelengths that contributed to the redshift fit |
|      | `obsconditions`              | None                    | Bit-coded of allowed observing conditions                    |
|      | `sv1_bgs_target`             | None                    | BGS (bright time program) target selection bitmask for SV1   |
|      | `cmx_target`                 | None                    | Target selection bitmask for commissioning                   |
|      | `pmra`                       | mas year-1              | Reference catalog proper motion in the RA direction          |
|      | `ncoeff`                     | None                    | Number of Redrock template coefficients                      |
|      | `sv3_bgs_target`             | None                    | BGS (bright time program) target selection bitmask for SV3   |
|      | `chi2`                       | None                    | Best fit chi-squared                                         |
|      | `sv1_scnd_target`            | None                    | Secondary target selection bitmask for SV1                   |
|      | `sv2_desi_target`            | None                    | DESI (dark time program) target selection bitmask for SV2    |
|      | `sv_primary`                 | None                    | Boolean flag (True/False) for the primary coadded spectrum in SV (SV1+2+3) |
|      | `mean_delta_y`               | mm                      | Mean (over exposures) fiber difference requested - actual CS5 Y location on focal plane |
|      | `fa_type`                    | None                    | Fiberassign internal target type (science, standard, sky, safe, suppsky) |
|      | `std_fiber_dec`              | arcsec                  | Standard deviation (over exposures) of DEC of actual fiber position |
|      | `tsnr2_lrg`                  | None                    | LRG template (S/N)^2 summed over B,R,Z                       |
|      | `numobs_init`                | None                    | Initial number of observations for target calculated across target selection bitmasks and OBSCONDITIONS |
|      | `subtype`                    | None                    | Spectral subtype                                             |
|      | `priority_init`              | None                    | Target initial priority from target selection bitmasks and OBSCONDITIONS |
|      | `mws_target`                 | None                    | Milky Way Survey targeting bits                              |
|      | `updated`                    | None                    | Date and time when the record was created/updated            |
|      | `wavemax`                    | AA                      | Maximum value of wavelength coverage in Angstroms            |
| *    | `redshift`                   | None                    | Measured redshift, or radial velocity in units of redshift   |
|      | `redshift_warning`           | None                    | Redshift warning bitmask measured by Redrock                 |
|      | `site`                       | None                    | [Categorical] Site of observations                           |
|      | `instrument`                 | None                    | [Categorical] Instrument used for observation                |
|      | `file`                       | None                    | File used for this record                                    |
|      | `telescope`                  | None                    | [Categorical] Telescope used for the observation             |
|      | `extra_files`                | None                    | Additional files read for ingest                             |
|      | `ra`                         | deg                     | Target right ascension                                       |
|      | `dateobs_center`             | None                    | Central value of the dateobs range                           |
|      | `data_release`               | None                    | [Categorical] Data release of a given record (with a common data structure, e.g., SDSS-DR16, BOSS-DR16) |
|      | `wavemin`                    | AA                      | Minimum value of wavelength coverage in Angstroms            |
|      | `datasetgroup`               | None                    | [Categorical] Common name for a coherent set of spectra published as a single release (e.g., SDSS_BOSS) |
|      | `spectype`                   | None                    | [Categorical] Spectral type from Redrock file (STAR, GALAXY, QSO) |
|      | `sparcl_id`                  | None                    | Universally Unique Identifier for spectrum in SPARCL         |
|      | `dateobs`                    | None                    | Range of observational dates for a given spectrum (single value for a single date) |
|      | `specid`                     | None                    | Dataset-specific spectrum identifier (may not be unique)     |
|      | `dec`                        | deg                     | Target declination                                           |
|      | `specprimary`                | None                    | Boolean flag (True/False) for the primary coadded spectrum for a given target object |
|      | `redshift_err`               | None                    | Uncertainty on the measured redshift, or radial velocity in units of redshift |
| *    | `targetid`                   | None                    | Unique DESI target ID                                        |
|      | `exptime`                    | second                  | Exposure time in seconds                                     |
| ?    | `model`                      | 1e-17 erg cm-2 s-1 AA-1 | Pipeline best model fit                                      |
| *    | `mask`                       | None                    | Mask value                                                   |
| *    | `wavelength`                 | AA                      | Wavelength in units of Angstroms                             |
| *    | `ivar`                       | 1e+34 cm4 s2 AA2 erg-2  | Inverse variance of flux                                     |
| *    | `wave_sigma`                 | AA                      | Line-spread function width (sigma of approximating Gaussian) in units of Angstroms |
| *    | `flux`                       | 1e-17 erg cm-2 s-1 AA-1 | Coadded calibrated flux                                      |