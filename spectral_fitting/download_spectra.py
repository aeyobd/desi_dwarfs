import tomllib
from astropy.table import Table
import numpy as np
import urllib.request
from pathlib import Path

def get_fsf_filename(props):
    survey = props["SURVEY"]
    program = props["PROGRAM"]
    healpix = props["HEALPIX"]
    group = healpix // 100
    
    url = f"https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v3.0/healpix/{survey}/{program}/{group}/{healpix}/fastspec-{survey}-{program}-{healpix}.fits.gz"

    return url



def get_desi_filename(props):
    survey = props["SURVEY"]
    program = props["PROGRAM"]
    healpix = props["HEALPIX"]
    group = healpix // 100
    
    url = f"https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/{survey}/{program}/{group}/{healpix}/coadd-{survey}-{program}-{healpix}.fits"

    return url

def get_props(dwarf_main, targetid):
    idx = np.where(dwarf_main["TARGETID"] == targetid)[0][0]

    props = dwarf_main[idx]
    return props


def main():
    with open("../analysis/selected_targetid.txt", "rb") as f:
        targetids = tomllib.load(f)["targetid"]

    dwarf_main = Table.read("../data/desi_dr1_dwarf_catalog.fits", hdu=1)

    for targetid in targetids:
        print("downloading", targetid)
        out_dir = Path(".", "example_fits", str(targetid))

        if not out_dir.exists():
            out_dir.mkdir()

        props = get_props(dwarf_main, targetid)

        desi_file = out_dir / "coadd.fits"
        if not desi_file.exists():
            url = get_desi_filename(props)
            urllib.request.urlretrieve(url, desi_file)

        fsf_file = out_dir / "fsf.fits"
        if not fsf_file.exists():
            url = get_fsf_filename(props)
            urllib.request.urlretrieve(url, fsf_file)




if __name__ == "__main__":
    main()
