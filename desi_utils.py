import re
import pyneb

emission_lines = [
    'nev_3346', 
    'ne_v_3426',
    'o_ii_3726',
    'o_ii_3729',
    'ne_iii_3869',
    'ne_iii_3967',
    'h_epsilon',
    'he_i_4026', 
    'h_delta',
    'h_gamma', 
    'o_iii_4363',
    'he_i_4471', 
    'he_ii_4686', 
    'h_beta',
    'o_iii_4959',
    'o_iii_5007',
    'n_ii_5755',
    'he_i_5876', 
    'o_i_6300', 
    's_iii_6312',
    'n_ii_6548',
    'h_alpha',  
    'n_ii_6583',
    's_ii_6716', 
    's_ii_6731',  
    'he_i_7065',
    'ar_iii_7136',
    'he_i_7281',
    'o_ii_7320', 
    'o_ii_7331',
    'ar_iii_7751',
    's_iii_9071', 
    's_iii_9533',
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




def plot_lines(meas, z=0):
    good_lines = retrieve_good_lines(meas)

    xlim = plt.gca().get_xlim()
    for line in emission_lines:
        λ = retrieve_wavelength(line)
        λ = λ * (1 + z)

        if xlim[0] <= λ <= xlim[1]:

            if line in good_lines:
                color = "k"
            else:
                color = "grey"
                
            plt.plot([λ, λ], [-10, -8], color=color)
            offset = line_position_shifts.get(line, (0,0))
            plt.annotate(nicer_label(line), (λ, -10), offset, 
                         va="top", ha="center", rotation=90,
                        xycoords = "data",
                        textcoords="offset fontsize", color=color)
