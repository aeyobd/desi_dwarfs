




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
