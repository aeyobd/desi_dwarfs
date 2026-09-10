import re

emline_cols = [
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


balmer_lines = {
    "halpha": 6565,
    "hbeta": 4861,
    "hgamma": 4340,
    "hdelta": 4102,
    "hepsilon": 3970,
}



line_labels = {
    "halpha": r"H $\alpha$",
    "hbeta": r"H $\beta$",
    "hgamma": r"H $\gamma$",
    "hdelta": r"H $\delta$",
    "hepsilon": r"H $\epsilon$",
    "nev": "[Ne V]",
    "oii": "[O II]",
    "neiii": "[Ne III]",
    "hei": "[He I]",
    "heii": "[He II]",
    "oiii": "[O III]",
    "nii": "[N II]",
    "oi": "[O I]",
    "siii": "[S III]",
    "sii": "[S II]", 
    "ariii": "[Ar III]",
}


def retrieve_wavelength(line):
    matches = re.findall(r"\d{4}", line)
    if len(matches) == 1:
        return int(matches[0])

    if line in balmer_lines.keys():
        return balmer_lines[line]


    raise Exception(f"Line not known {line}")





def retrieve_good_lines(galaxy, snr_min=3):
    good_lines = []
    for line in emline_cols:
        if galaxy[line + "_flux"] / galaxy[line + "_fluxerr"] > snr_min:
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
    for line in emline_cols:
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
