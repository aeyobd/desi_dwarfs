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
from astropy.io import fits
import numpy as np


outfile = Path("../data/dwarfs_good.fits")

if outfile.exists():
    os.remove(outfile)


# ## Load the catalogue
dwarf_main = Table.read("../data/desi_dr1_dwarf_catalog.fits", hdu=1)
dwarf_z = Table.read("../data/desi_dr1_dwarf_catalog.fits", hdu=2)
dwarf_phot = Table.read("../data/desi_dr1_dwarf_catalog.fits", hdu=3)
dwarf_fsf = Table.read("../data/desi_dr1_dwarf_catalog.fits", hdu=4)


assert np.all(dwarf_fsf["TARGETID"] == dwarf_z["TARGETID"])
assert np.all(dwarf_fsf["TARGETID"] == dwarf_main["TARGETID"])
assert np.all(dwarf_fsf["TARGETID"] == dwarf_phot["TARGETID"])


# check all primary catalogue and no duplicates
filt_primary = dwarf_z["ZCAT_PRIMARY"]
assert np.all(filt_primary)
assert len(filt_primary) == len(np.unique(dwarf_main["TARGETID"]))


filt = filt_primary
filt &= dwarf_main["DWARF_MASKBIT"] == 0
filt &= dwarf_main["MSTAR_MASKBIT"] == 0


filt &= dwarf_phot["FRACFLUX_G"] < 0.3
filt &= dwarf_phot["FRACFLUX_R"] < 0.3
filt &= dwarf_phot["FRACFLUX_Z"] < 0.3


print("after fracflux cut", np.sum(filt))

required_lines = ["NII_6584", "HALPHA", "OIII_5007", "HBETA"]
for line in required_lines:
    filt &= dwarf_fsf[f"{line}_FLUX"] / dwarf_fsf[f"{line}_FLUX_IVAR"]**-0.5 > 3
    print(line, np.sum(filt))

filt &= dwarf_fsf["SNR_R"] > 2
print("SNR", np.sum(filt))


good_main = dwarf_main[filt]
good_z = dwarf_z[filt]
good_phot = dwarf_phot[filt]
good_fsf = dwarf_fsf[filt]

print("writing catalogue")
hdr = fits.Header()
hdus = fits.HDUList([
    fits.PrimaryHDU(header=hdr),
    fits.BinTableHDU(good_main),
    fits.BinTableHDU(good_z),
    fits.BinTableHDU(good_phot),
    fits.BinTableHDU(good_fsf),
])


hdus.writeto(outfile)
