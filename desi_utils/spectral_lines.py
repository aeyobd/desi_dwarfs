import re
import pyneb

EMISSION_LINES = [
    'Ne5_3346A', 
    'Ne5_3426A',
    'O2_3726A',
    'O2_3729A',
    'Ne3_3869A',
    'Ne3_3967A',
    'H1r_3970A',
    'He1_4026A', 
    'H1r_4102A',
    'H1r_4340A', 
    'O3_4363A',
    'He1_4471A', 
    'He2_4686A', 
    'H1r_4861A',
    'O3_4959A',
    'O3_5007A',
    'N2_5755A',
    'He1_5876A', 
    'O1_6300A', 
    'S3_6312A',
    'N2_6548A',
    'H1r_6563A',  
    'N2_6583A',
    'S2_6716A', 
    'S2_6731A',  
    'He1_7065A',
    'Ar3_7136A',
    'He1_7281A',
    'O2_7320A', 
    'O2_7331A',
    'Ar3_7751A',
    'S3_9071A', 
    'S3_9533A',
]





balmer_line_wavelengths = {
    "h_alpha": 6563,
    "h_beta": 4861,
    "h_gamma": 4340,
    "h_delta": 4102,
    "h_epsilon": 3970,
}



line_labels = {
    "h_alpha": r"H $\alpha$",
    "h_beta": r"H $\beta$",
    "h_gamma": r"H $\gamma$",
    "h_delta": r"H $\delta$",
    "h_epsilon": r"H $\epsilon$",
}


def get_line_element(emline):
    return emline.split("_")[0]


def get_wavelength(emline):
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
    

def get_species(emline):
    if emline in balmer_line_wavelengths.keys():
        return "ii"
    else:
        return emline.split("_")[1]


def to_elsm_format(emline):
    return emline.replace("_", "")


def get_pyneb_line(emline: str):
    if emline in balmer_line_wavelengths.keys():
        wave = get_named_wavelength(emline)
        Ha = pyneb.EmissionLine(label=f"H1r_{wave}A")

    ele = get_line_element(emline).title()
    wave = get_line_named_wavelength(emline)
    spec = roman_num_to_int(get_species(emline))

    return pyneb.EmissionLine(ele, spec, wave)





def retrieve_good_lines(galaxy, snr_min=3):
    good_lines = []
    for line in emission_lines:
        line_elsm = to_elsm_format(line)
        snr = galaxy[line_elsm + "_flux"] / galaxy[line_elsm + "_fluxerr"]
        if snr > snr_min:
            good_lines.append(line)

    return good_lines



def nicer_label(line):
    if line in line_labels:
        return line_labels[line]
    
    else: 
        λ = retrieve_wavelength(line)
        species = line.replace(str(λ), "")
        species = line_labels[species]
        return species + r" $\lambda$" + str(λ)