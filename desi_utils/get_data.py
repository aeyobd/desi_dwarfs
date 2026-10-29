from astropy.utils.data import download_file
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import CCDData
from sparcl.client import SparclClient
# from . import Spectrum


def download_ls_image(ra, dec):
    "retrieve the LSS image at default size located at the given ra and dec in giz colours"

    link = f"https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer=ls-dr9&pixscale=0.262&bands=grz"

    print(link)
    filename = download_file(link, cache=True, show_progress=False, timeout=120)
    hdu = fits.open(filename)[0]
    image = hdu.data
    hdr = hdu.header
    wcs = WCS(hdr)

    return CCDData(image, header=hdr, wcs=wcs, unit="nanomaggy")


def get_spectra(targetid):
    if hasattr(targetid, "__len__"):
        targetid = [int(i) for i in targetid]
    else:
        targetid = [int(targetid)]


    inc = "ALL"


    client = SparclClient()

    galspec_raw = client.retrieve_by_specid(specid_list = targetid,
        include = inc,
        dataset_list = ['DESI-DR1'], 
        fmt="specutils")

    return galspec_raw
