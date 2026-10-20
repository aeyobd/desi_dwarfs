#!/usr/bin/env python

# # Selecting a high-quality sample of dwarf galaxies
# 
# In this notebook, we apply basic selection cuts to our dwarf galaxy catalogue. 
# 
# Here, we further filter by some basic quality cuts:
# - Including only main zurvey targets (`zcat_primary`)
# - `snr_median` > 2 (i.e. some measured continuum)
# - SNR  > 3 for the BPT lines (\[O III\] 5007$\AA$, \[N II\] 6583$\AA$, H$\alpha$ and H$\beta$
# - H$\alpha$ equivalent width greater than 14 (avoiding diffuse intergalactic medium contamination)

import os

from pathlib import Path
from astropy.table import Table, unique
import numpy as np


outfile = Path("../data/dwarfs_good.fits")

if outfile.exists():
    os.remove(outfile)


# ## Load the catalogue
dwarf_cat = Table.read("../data/stellar_mass_emline_dwarfs.fits")


# ## Remove duplicated targets
filt_primary = dwarf_cat["zcat_primary"] == b't'
print(len(dwarf_cat))
print(np.sum(filt_primary))
primary_cat = dwarf_cat[filt_primary]
primary_cat.sort("snr_med", reverse=True)
primary_cat = unique(primary_cat, "targetid")
print("after removing duplicates", len(primary_cat))



filt = primary_cat["zcat_primary"] == b't'

filt &= primary_cat["fracflux_g"] < 0.3
filt &= primary_cat["fracflux_r"] < 0.3
filt &= primary_cat["fracflux_z"] < 0.3
print("after fracflux cut", np.sum(filt))

required_lines = ["nii6583", "halpha", "oiii5007", "hbeta"]
for line in required_lines:
    filt &= primary_cat[f"{line}_flux"] / primary_cat[f"{line}_fluxerr"] > 3
    print(line, np.sum(filt))

filt &= primary_cat["snr_med"] > 2
print("SNR", np.sum(filt))

filt &= primary_cat["halpha_ew"] > 14
print("Halpha", np.sum(filt))

print("Halpha EW cut", sum(filt))


filt &= primary_cat["logm_cigale"] < 9.0
print("mass cut", np.sum(filt))
good_dwarfs = primary_cat[filt]

print("writing catalogue")
good_dwarfs.write(outfile)

