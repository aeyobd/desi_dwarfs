import numpy as np

from .spectral_lines import get_named_wavelength, BALMER_WAVELENGTHS



MEL_LINES = [
    'nev3346', 
    'nev3426',
    'oii3726',
    'oii3729',
    'neiii3869',
    'neiii3967',
    'hepsilon',
    'hei4026', 
    'hdelta',
    'hgamma', 
    'oiii4363',
    'hei4471', 
    'heii4686', 
    'hbeta',
    'oiii4959',
    'oiii5007',
    'nii5755',
    'hei5876', 
    'oi6300', 
    'siii6312',
    'nii6548',
    'halpha',  
    'nii6583',
    'sii6716', 
    'sii6731',  
    'hei7065',
    'ariii7136',
    'hei7281',
    'oii7320', 
    'oii7331',
    'ariii7751',
    'siii9071', 
    'siii9533',
]

FSF_LINES = ['LYALPHA',
     'NV_1240',
     'OI_1304',
     'SILIV_1396',
     'CIV_1549',
     'HEII_1640',
     'ALIII_1857',
     'SILIII_1892',
     'CIII_1908',
     'MGII_2796',
     'MGII_2803',
     'NEV_3346',
     'NEV_3426',
     'OII_3726',
     'OII_3729',
     'NEIII_3869',
     'H6',
     'H6_BROAD',
     'HEPSILON',
     'HEPSILON_BROAD',
     'HDELTA',
     'HDELTA_BROAD',
     'HGAMMA_BROAD',
     'HGAMMA',
     'OIII_4363',
     'HEI_4471',
     'HEII_4686',
     'HBETA',
     'HBETA_BROAD',
     'OIII_4959',
     'OIII_5007',
     'NII_5755',
     'HEI_5876',
     'OI_6300',
     'SIII_6312',
     'NII_6548',
     'HALPHA',
     'HALPHA_BROAD',
     'NII_6584',
     'SII_6716',
     'SII_6731',
     'ARIII_7135',
     'OII_7320',
     'OII_7330',
     'SIII_9069',
     'SIII_9532'
]

FSF_NAMED_LINES = {
    "LYALPHA": 1215.7,
    "H6": 3889.064
}

def to_mel_line(line):
    wave = get_named_wavelength(line)
    mel_waves = [get_named_wavelength(line) for line in MEL_LINES]
    idx = np.argmin(np.abs(np.array(mel_waves) - wave))

    mel_wave = mel_waves[idx]
    if abs(mel_wave - wave) / wave < 0.01:
        mel_line = MEL_LINES[idx]
    else:
        print("no match for ", line)
        return None

    return mel_line



def get_fsf_wave(line):
    line = line.replace("_BROAD", "")
    if line in FSF_NAMED_LINES:
        return FSF_NAMED_LINES[line]
    elif line.lower() in BALMER_WAVELENGTHS:
        return BALMER_WAVELENGTHS[line.lower()]
    else:
        return get_named_wavelength(line)


def to_fsf_line(line):
    wave = get_named_wavelength(line)
    fsf_waves = [get_fsf_wave(line) for line in FSF_LINES]
    idx = np.argmin(np.abs(np.array(fsf_waves) - wave))

    fsf_wave = fsf_waves[idx]
    if abs(fsf_wave - wave) / wave < 0.01:
        fsf_line = FSF_LINES[idx]
    else:
        print("no match for ", line)
        return None

    return fsf_line
