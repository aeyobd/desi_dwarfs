import re
import numpy as np
import pyneb

EMISSION_LINES = [
    'Ne5_3346A', 
    'Ne5_3426A',
    'O2_3726A',
    'O2_3729A',
    'Ne3_3869A',
    'Ne3_3968A',
    'H1r_3970A',
    'He1r_4026A', 
    'H1r_4102A',
    'H1r_4341A', 
    'O3_4363A',
    'He1r_4471A', 
    'He2r_4686A', 
    'Ar4_4711A', # new
    'Ar4_4740A', # new
    'H1r_4861A',
    'O3_4959A',
    'O3_5007A',
    'N2_5755A',
    'He1r_5876A', 
    'O1_6300A', 
    'S3_6312A',
    'O1_6364A', # new
    'N2_6548A',
    'H1r_6563A',  
    'N2_6584A',
    'S2_6716A', 
    'S2_6731A',  
    'He1r_7065A',
    'Ar3_7136A',
    'He1r_7281A',
    'O2_7320A', 
    'O2_7331A',
    'Ar3_7751A',
    'S3_9069A', 
    'S3_9531A',
]





# Helium I wavelengths from NIST as pyneb is not setup properly
HE_WAVELENGTHS = {
        'He1r_4026A': 4026.1914, 
        'He1r_4471A': 4471.4802, 
        'He1r_5876A': 5875.621, 
        'He1r_7065A': 7065.190,
        'He1r_7281A': 7281.349,
}


def get_line_element(emline):
    return emline.split("_")[0]


def get_wavelength(emline, *, vacuum=True, wavetol=0.005):
    """Return the (vacuum) wavelength of a pyneb-format emline """
    line = get_pyneb_line(emline)
    if emline.startswith("He1r"):
        wave = pyneb.utils.physics.airtovac(HE_WAVELENGTHS[emline])
    else:
        if is_recomb_line(line):
            atom = pyneb.RecAtom(line.elem, line.spec)
        else:
            atom = pyneb.Atom(line.elem, line.spec)

        i, j = atom.getTransition(line.wave, maxErrorA=wavetol)
        Ei = atom.getEnergy(i)

        Ej = atom.getEnergy(j)
        assert Ei > Ej
        wave = 1 / (Ei - Ej)
        assert np.abs(1 - wave / line.wave) < wavetol

    if not vacuum:
        wave = pyneb.utils.physics.vactoair(wave)

    return wave

def is_recomb_line(line):
    return line.atom.endswith("r")


def get_named_wavelength(emline):
    if get_line_element(emline) == "h":
        return  balmer_line_wavelengths[emline]
        
    matches = re.findall(r"\d{4}", emline)
    if len(matches) == 1:
        return int(matches[0])
    else:
        print("bad format for line: ", emline)
        print("we expect the line to be of the form `el_iii_1000` or a balmer line")


def roman_num_to_int(roman):
    return {
        "i": 1,
        "ii": 2,
        "iii": 3,
        "iv": 4,
        "v": 5,
        "vi": 6,
        "vii": 7,
        "viii": 8,
        "ix": 9,
        "x": 10
    }[roman]
    



def get_pyneb_line(emline: str):
    "Return the pyneb.EmissionLine object from the line name"
    return pyneb.EmissionLine(label=emline)

def to_elsm_format(emline):
    "Convert pyneb formatted lines to ELSM formatted lines"

    return emline.replace("_", "")



def get_good_lines_elsm(galaxy_properties_elsm, snr_min=3):
    """
    Get the lines in an emission line stellar mass format dictionary-like
    object which have the given SNR.
    """
    good_lines = []
    for line in EMISSION_LINES:
        line_elsm = to_elsm_format(line)
        snr = galaxy_properties_elsm[line_elsm + "_flux"] / galaxy_properties_elsm[line_elsm + "_fluxerr"]
        if snr > snr_min:
            good_lines.append(line)

    return good_lines


def int_to_roman(i):
    return ["i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x"][i-1]
            
def format_line(emline):

    line = get_pyneb_line(emline)
    roman = int_to_roman(int(line.spec))
    ele = line.elem
    wave = get_named_wavelength(emline)

    species  = f"{ele.title()} {roman.upper()}"
    if not is_recomb_line(line):
        species = f"[{species}]"


    return f"{species} $\\lambda${wave}"
