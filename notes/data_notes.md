| Storage | Field name | Units | Description | Default | All | Provenance field name | File | HDU |
|---|---|---|---|---|---|---|---|---|
| A | `tsnr2_qso` | None | QSO template (S/N)^2 summed over B,R,Z | False | True | `TSNR2_QSO` | zall-pix-iron | 1 |
| A | `sv2_scnd_target` | None | Secondary target selection bitmask for SV2 | False | True | `SV2_SCND_TARGET` | zall-pix-iron | 1 |
| A | `survey` | None | Survey name | False | True | `SURVEY` | zall-pix-iron | 1 |
| A | `main_nspec` | None | Number of coadded spectra for this TARGETID in Main survey | False | True | `MAIN_NSPEC` | zall-pix-iron | 1 |
| A | `bgs_target` | None | BGS (Bright Galaxy Survey) target selection bitmask | False | True | `BGS_TARGET` | zall-pix-iron | 1 |
| A | `program` | None | DESI program type - BRIGHT, DARK, BACKUP, OTHER | False | True | `PROGRAM` | zall-pix-iron | 1 |
| A | `zcat_nspec` | None | Number of coadded spectra for this TARGETID in this zcatalog | False | True | `ZCAT_NSPEC` | zall-pix-iron | 1 |
| A | `desi_target` | None | DESI (dark time program) target selection bitmask | False | True | `DESI_TARGET` | zall-pix-iron | 1 |
| A | `sv3_scnd_target` | None | Secondary target selection bitmask for SV3 | False | True | `SV3_SCND_TARGET` | zall-pix-iron | 1 |
| A | `tsnr2_gpbbackup` | None | GPBBACKUP template (S/N)^2 summed over BRZ | False | True | `TSNR2_GPBBACKUP` | zall-pix-iron | 1 |
| A | `sv_nspec` | None | Number of coadded spectra for this TARGETID in SV (SV1+2+3) | False | True | `SV_NSPEC` | zall-pix-iron | 1 |
| A | `coadd_fiberstatus` | None | Bitwise AND of the FIBERSTATUS values of the spectra used in the coadd | False | True | `COADD_FIBERSTATUS` | zall-pix-iron | 1 |
| A | `pmdec` | mas year-1 | Reference catalog proper motion in the Dec direction | False | True | `PMDEC` | zall-pix-iron | 1 |
| A | `std_fiber_ra` | arcsec | Standard deviation (over exposures) of RA of actual fiber position | False | True | `STD_FIBER_RA` | zall-pix-iron | 1 |
| A | `sv2_mws_target` | None | MWS (bright time program) target selection bitmask for SV2 | False | True | `SV2_MWS_TARGET` | zall-pix-iron | 1 |
| A | `tsnr2_gpbbright` | None | GPBBRIGHT template (S/N)^2 summed over B,R,Z | False | True | `TSNR2_GPBBRIGHT` | zall-pix-iron | 1 |
| A | `sv1_mws_target` | None | MWS (bright time program) target selection bitmask for SV1 | False | True | `SV1_MWS_TARGET` | zall-pix-iron | 1 |
| A | `objtype` | None | Object type: TGT, SKY, NON, BAD | False | True | `OBJTYPE` | zall-pix-iron | 1 |
| A | `plate_ra` | deg | Right Ascension to be used by PlateMaker | False | True | `PLATE_RA` | zall-pix-iron | 1 |
| A | `coadd_numexp` | None | Number of exposures included in the coadd for this target | False | True | `COADD_NUMEXP` | zall-pix-iron | 1 |
| A | `sv3_desi_target` | None | DESI (dark time program) target selection bitmask for SV3 | False | True | `SV3_DESI_TARGET` | zall-pix-iron | 1 |
| A | `rms_delta_y` | mm | RMS (over exposures) of the fiber difference between measured and requested CS5 Y location on focal plane | False | True | `RMS_DELTA_Y` | zall-pix-iron | 1 |
| A | `ref_epoch` | year | Reference epoch for Gaia/Tycho astrometry. Typically 2015.5 for Gaia | False | True | `REF_EPOCH` | zall-pix-iron | 1 |
| A | `rms_delta_x` | mm | RMS (over exposures) of the fiber difference between measured and requested CS5 X location on focal plane | False | True | `RMS_DELTA_X` | zall-pix-iron | 1 |
| A | `sv3_mws_target` | None | MWS (bright time program) target selection bitmask for SV3 | False | True | `SV3_MWS_TARGET` | zall-pix-iron | 1 |
| A | `sv2_bgs_target` | None | BGS (bright time program) target selection bitmask for SV2 | False | True | `SV2_BGS_TARGET` | zall-pix-iron | 1 |
| A | `tsnr2_elg` | None | ELG template (S/N)^2 summed over B,R,Z | False | True | `TSNR2_ELG` | zall-pix-iron | 1 |
| A | `mean_fiber_ra` | deg | Mean (over exposures) RA of actual fiber position | False | True | `MEAN_FIBER_RA` | zall-pix-iron | 1 |
| A | `sv1_desi_target` | None | DESI (dark time program) target selection bitmask for SV1 | False | True | `SV1_DESI_TARGET` | zall-pix-iron | 1 |
| A | `mean_delta_x` | mm | Mean (over exposures) fiber difference between measured and requested CS5 X location on focal plane | False | True | `MEAN_DELTA_X` | zall-pix-iron | 1 |
| A | `healpix` | None | HEALPixel containing this location at NSIDE=64 in the NESTED scheme | False | True | `HEALPIX` | zall-pix-iron | 1 |
| A | `fa_target` | None | Targeting bit internally used by fiberassign (linked with FA_TYPE) | False | True | `FA_TARGET` | zall-pix-iron | 1 |
| A | `tsnr2_bgs` | None | BGS template (S/N)^2 summed over B,R,Z | False | True | `TSNR2_BGS` | zall-pix-iron | 1 |
| A | `mean_mjd` | days | Mean (over exposures) of the Modified Julian Date (MJD) when the shutter was opened | False | True | `MEAN_MJD` | zall-pix-iron | 1 |
| A | `tsnr2_lya` | None | LYA template (S/N)^2 summed over B,R,Z | False | True | `TSNR2_LYA` | zall-pix-iron | 1 |
| A | `deltachi2` | None | Delta-chi-squared for template fit from Redrock | False | True | `DELTACHI2` | zall-pix-iron | 1 |
| A | `scnd_target` | None | Target selection bitmask for secondary programs | False | True | `SCND_TARGET` | zall-pix-iron | 1 |
| A | `coadd_numtile` | None | Number of different tiles included in the coadd for this target | False | True | `COADD_NUMTILE` | zall-pix-iron | 1 |
| A | `coeff` | None | Redrock template coefficients | False | True | `COEFF` | zall-pix-iron | 1 |
| A | `coadd_numnight` | None | Number of different nights included in the coadd for this target | False | True | `COADD_NUMNIGHT` | zall-pix-iron | 1 |
| A | `tsnr2_gpbdark` | None | GPBDARK template (S/N)^2 summed over BRZ | False | True | `TSNR2_GPBDARK` | zall-pix-iron | 1 |
| A | `plate_dec` | deg | Declination to be used by PlateMaker | False | True | `PLATE_DEC` | zall-pix-iron | 1 |
| A | `main_primary` | None | Boolean flag (True/False) for the primary coadded spectrum in Main survey | False | True | `MAIN_PRIMARY` | zall-pix-iron | 1 |
| A | `mean_fiber_dec` | deg | Mean (over exposures) DEC of actual fiber position | False | True | `MEAN_FIBER_DEC` | zall-pix-iron | 1 |
| A | `subpriority` | None | Random subpriority [0-1] to break assignment ties | False | True | `SUBPRIORITY` | zall-pix-iron | 1 |
| A | `mean_psf_to_fiber_specflux` | None | Mean of input exposures fraction of light from point-like source captured by 1.5 arcsec diameter fiber given atmospheric seeing | False | True | `MEAN_PSF_TO_FIBER_SPECFLUX` | zall-pix-iron | 1 |
| A | `desiname` | None | DESI JXXX.XXXX[+/-]YY.YYYY where X, Y are TARGET_{RA, DEC} truncated to 4 decimals | False | True | `DESINAME` | zall-pix-iron | 1 |
| A | `spgrpval` | None | Value corresponding to the grouping type indicated in the spgrp column. This is the same as the healpix column | False | True | `SPGRPVAL` | zall-pix-iron | 1 |
| A | `npixels` | None | Number of unmasked wavelengths that contributed to the redshift fit | False | True | `NPIXELS` | zall-pix-iron | 1 |
| A | `obsconditions` | None | Bit-coded of allowed observing conditions | False | True | `OBSCONDITIONS` | zall-pix-iron | 1 |
| A | `sv1_bgs_target` | None | BGS (bright time program) target selection bitmask for SV1 | False | True | `SV1_BGS_TARGET` | zall-pix-iron | 1 |
| A | `cmx_target` | None | Target selection bitmask for commissioning | False | True | `CMX_TARGET` | zall-pix-iron | 1 |
| A | `pmra` | mas year-1 | Reference catalog proper motion in the RA direction | False | True | `PMRA` | zall-pix-iron | 1 |
| A | `ncoeff` | None | Number of Redrock template coefficients | False | True | `NCOEFF` | zall-pix-iron | 1 |
| A | `sv3_bgs_target` | None | BGS (bright time program) target selection bitmask for SV3 | False | True | `SV3_BGS_TARGET` | zall-pix-iron | 1 |
| A | `chi2` | None | Best fit chi-squared | False | True | `CHI2` | zall-pix-iron | 1 |
| A | `sv1_scnd_target` | None | Secondary target selection bitmask for SV1 | False | True | `SV1_SCND_TARGET` | zall-pix-iron | 1 |
| A | `sv2_desi_target` | None | DESI (dark time program) target selection bitmask for SV2 | False | True | `SV2_DESI_TARGET` | zall-pix-iron | 1 |
| A | `sv_primary` | None | Boolean flag (True/False) for the primary coadded spectrum in SV (SV1+2+3) | False | True | `SV_PRIMARY` | zall-pix-iron | 1 |
| A | `mean_delta_y` | mm | Mean (over exposures) fiber difference requested - actual CS5 Y location on focal plane | False | True | `MEAN_DELTA_Y` | zall-pix-iron | 1 |
| A | `fa_type` | None | Fiberassign internal target type (science, standard, sky, safe, suppsky) | False | True | `FA_TYPE` | zall-pix-iron | 1 |
| A | `std_fiber_dec` | arcsec | Standard deviation (over exposures) of DEC of actual fiber position | False | True | `STD_FIBER_DEC` | zall-pix-iron | 1 |
| A | `tsnr2_lrg` | None | LRG template (S/N)^2 summed over B,R,Z | False | True | `TSNR2_LRG` | zall-pix-iron | 1 |
| A | `numobs_init` | None | Initial number of observations for target calculated across target selection bitmasks and OBSCONDITIONS | False | True | `NUMOBS_INIT` | zall-pix-iron | 1 |
| A | `subtype` | None | Spectral subtype | False | True | `SUBTYPE` | zall-pix-iron | 1 |
| A | `priority_init` | None | Target initial priority from target selection bitmasks and OBSCONDITIONS | False | True | `PRIORITY_INIT` | zall-pix-iron | 1 |
| A | `mws_target` | None | Milky Way Survey targeting bits | False | True | `MWS_TARGET` | zall-pix-iron | 1 |
| C | `updated` | None | Date and time when the record was created/updated | False | False | `N/A` | N/A | N/A |
| C | `wavemax` | AA | Maximum value of wavelength coverage in Angstroms | False | True | `N/A` | N/A | N/A |
| C | `redshift` | None | Measured redshift, or radial velocity in units of redshift | False | True | `Z` | zall-pix-iron | 1 |
| C | `redshift_warning` | None | Redshift warning bitmask measured by Redrock | False | True | `ZWARN` | zall-pix-iron | 1 |
| C | `site` | None | [Categorical] Site of observations | False | True | `N/A` | N/A | N/A |
| C | `instrument` | None | [Categorical] Instrument used for observation | False | True | `N/A` | N/A | N/A |
| C | `file` | None | File used for this record | False | False | `N/A` | N/A | N/A |
| C | `telescope` | None | [Categorical] Telescope used for the observation | False | True | `N/A` | N/A | N/A |
| C | `extra_files` | None | Additional files read for ingest | False | False | `N/A` | N/A | N/A |
| C | `ra` | deg | Target right ascension | True | True | `TARGET_RA` | camcoadd | 1 |
| C | `dateobs_center` | None | Central value of the dateobs range | False | True | `MIN_MJD, MAX_MJD` | zall-pix-iron | 1 |
| C | `data_release` | None | [Categorical] Data release of a given record (with a common data structure, e.g., SDSS-DR16, BOSS-DR16) | False | True | `N/A` | N/A | N/A |
| C | `wavemin` | AA | Minimum value of wavelength coverage in Angstroms | False | True | `N/A` | N/A | N/A |
| C | `datasetgroup` | None | [Categorical] Common name for a coherent set of spectra published as a single release (e.g., SDSS_BOSS) | False | True | `N/A` | N/A | N/A |
| C | `spectype` | None | [Categorical] Spectral type from Redrock file (STAR, GALAXY, QSO) | False | True | `SPECTYPE` | zall-pix-iron | 1 |
| C | `sparcl_id` | None | Universally Unique Identifier for spectrum in SPARCL | True | True | `N/A` | N/A | N/A |
| C | `dateobs` | None | Range of observational dates for a given spectrum (single value for a single date) | False | True | `MIN_MJD, MAX_MJD` | zall-pix-iron | 1 |
| C | `specid` | None | Dataset-specific spectrum identifier (may not be unique) | True | True | `TARGETID` | camcoadd | 1 |
| C | `dec` | deg | Target declination | True | True | `TARGET_DEC` | camcoadd | 1 |
| C | `specprimary` | None | Boolean flag (True/False) for the primary coadded spectrum for a given target object | False | True | `ZCAT_PRIMARY` | zall-pix-iron | 1 |
| C | `redshift_err` | None | Uncertainty on the measured redshift, or radial velocity in units of redshift | False | True | `ZERR` | zall-pix-iron | 1 |
| C | `targetid` | None | Unique DESI target ID | False | True | `TARGETID` | camcoadd | 1 |
| C | `exptime` | second | Exposure time in seconds | False | True | `COADD_EXPTIME` | camcoadd | 1 |
| S | `model` | 1e-17 erg cm-2 s-1 AA-1 | Pipeline best model fit | False | True | `BRZ_MODEL` | camcoadd | 8 |
| S | `mask` | None | Mask value | False | True | `BRZ_MASK` | camcoadd | 6 |
| S | `wavelength` | AA | Wavelength in units of Angstroms | True | True | `BRZ_WAVELENGTH` | camcoadd | 3 |
| S | `ivar` | 1e+34 cm4 s2 AA2 erg-2 | Inverse variance of flux | False | True | `BRZ_IVAR` | camcoadd | 5 |
| S | `wave_sigma` | AA | Line-spread function width (sigma of approximating Gaussian) in units of Angstroms | False | True | `BRZ_wavedisp` | camcoadd | 10 |
| S | `flux` | 1e-17 erg cm-2 s-1 AA-1 | Coadded calibrated flux | True | True | `BRZ_FLUX` | camcoadd | 4 |
