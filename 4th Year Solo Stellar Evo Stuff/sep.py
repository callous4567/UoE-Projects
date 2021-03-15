import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Patch
plt.rcParams["font.family"] = "serif"
import read_mist_models
from matplotlib import ticker as pltick
# Specify root directory for code. EEP should be under rootdir\\MIST
rootdir = "D:\\Prog\\pycharm\\StellarEvolutionProject"

# Specify phases. "extra" = something useful.
available_phases = [-1,0,2,3,4,5,6,9]
available_phases_str = [("{0:.0f}").format(d) for d in available_phases]
available_phases_extra = [-1,0,2,3,4,5,6,9,10]
phase_names = ["PMS", "MS", "RGB", "HEB", "EAB", "TPAB", "PAB", "WRS"]

# Colours for HR and linear colours for the param plots
colours = ['sienna','red','gold','darkorange','lime','aquamarine','mediumorchid','black', "hotpink", "deepskyblue", "olive"]
phase_colours_extralinear = ['thistle','plum','orchid','purple','deeppink','hotpink','crimson','sienna']
phase_colours_linear = ['lightcoral','indianred','brown','firebrick','maroon','darkred','orangered','black']
spancolours = ['']

# Phases in number form
"""
Pre-main-sequence
Main-sequence
Red-giant-branch
Helium-burning-core
Early-asymptotic-branch
Thermal-pulse-asymptotic-branch
Wolf-rayet-star
"""

# Grab EEP from directory
def model_grab(mass_int):
    if mass_int >= 100:
        mass_str = ("{0:.0f}").format(mass_int)
        mass_filename = rootdir + "\\" + "MIST" + "\\" + ("{}0000M.track.eep").format(mass_str)
        eep = read_mist_models.EEP(mass_filename)
        return eep
    if 100 > mass_int >= 10:
        mass_str = ("{0:.0f}").format(mass_int)
        mass_filename = rootdir + "\\" + "MIST" + "\\" + ("0{}0000M.track.eep").format(mass_str)
        eep = read_mist_models.EEP(mass_filename)
        return eep
    if 10 > mass_int >= 1:
        mass_str = ("0{0:.0f}").format(mass_int)
        mass_filename = rootdir + "\\" + "MIST" + "\\" + ("0{}0000M.track.eep").format(mass_str)
        eep = read_mist_models.EEP(mass_filename)
        return eep
    if 1 > mass_int >= 0:
        mass_int = mass_int*10
        mass_str = ("00{0:.0f}").format(mass_int)
        mass_filename = rootdir + "\\" + "MIST" + "\\" + ("0{}0000M.track.eep").format(mass_str)
        eep = read_mist_models.EEP(mass_filename)
        return eep


# Plot HR for EEP in directory
def hr_plot(mass_int, xlim, ylim):
    model = model_grab(mass_int)
    rcParams['figure.figsize'] = 20, 20  # this will set the figure size and keep edges from being cut off
    rcParams.update({'font.size': 22})  # this sets the fontsize in your plot to a reasonable value

    """
    # Grab ages + limiters + etc
    phase_delimiters = []
    # Identify the phase ranges
    for phaseval in available_phases:
        for num, val in enumerate(model.eeps['phase']):
            if (val+0.1) >= phaseval:
                phase_delimiters.append(num)
                print("yep")
                break
    phase_delimiters.append(len(model.eeps['phase'])-1)
    phase_len = len(phase_delimiters) - 1
    age_delimiters = [model.eeps['star_age'][d] for d in phase_delimiters]
    age_lengths = []
    for d in range(phase_len):
        age_lengths.append(age_delimiters[d+1] - age_delimiters[d])
    age_delimiters = [("{:.5e}").format(d) for d in age_delimiters]
    age_lengths = [("{:.5e}").format(d) for d in age_lengths]"""

    # Identify the phases present, their age limits, and create patch list.

    phase_strings = [("{0:.0f}").format(d) for d in model.eeps['phase']]
    phase_strings.append("99") # Random limiter for loops that operate for == statements.

    # First grab patches for all present phases.
    present_phases = []
    for num,phase in enumerate(available_phases_str):
        if phase in phase_strings:
            present_phases.append(phase)

    # Then estimate the ranges for each phase that is present, going left-to-right.
    # You can use find() to do this. This works, so we won't clean it up/change it, but find() would make it easier.
    phase_ranges = []
    for phase in present_phases:
        start, finish = None,None
        # Find start of phase
        for num,d in enumerate(phase_strings):
            if d == phase:
                start = num
                break
        # Find end of phase.
        for num,d in enumerate(phase_strings):
            if phase_strings[num+start] != phase:
                finish = num+start-1
                break
        phase_ranges.append([start,finish])

    # Estimate the length of time spent in each phase.]
    phase_times = [model.eeps['star_age'][d[1]] - model.eeps['star_age'][d[0]] for d in phase_ranges]
    phase_lengths_str = [("{:.5e}").format(d) for d in phase_times]

    # Now we iterate over the present phases and set up the necessary patches for plotting.
    patch_elements_phases = []
    patch_elements_ages = []
    for num,phase in enumerate(available_phases_str):
        for index, value in enumerate(present_phases):
            if value == phase:
                patch_elements_phases.append(Patch(edgecolor='black',facecolor=colours[num],label=("{0}").format(phase_names[num])))
                patch_elements_ages.append(Patch(edgecolor='black', facecolor=colours[num],label=phase_lengths_str[index]))



    # [elements_1, loc_1, bbox_tuple_1, title1]
    list1 = [patch_elements_phases, "lower left", (1,0.51), "Phases"]
    list2 = [patch_elements_ages, "upper left", (1,0.51), "Time (yr)"]
    resize_scale = 0.9
    full = [list1, list2, resize_scale]


    model.plot_HR(custom_arguments=full, phases=available_phases, phasecolor=colours[0:8])

    plt.savefig(str(mass_int) + "HR.png")
    plt.show()
#hr_plot(1, False, False)




# Multiple HR plot, specify list of masses
def hr_multiplot(mass_ints, xlim, ylim):
    rcParams['figure.figsize'] = 20, 20  # this will set the figure size and keep edges from being cut off
    rcParams.update({'font.size': 22})  # this sets the fontsize in your plot to a reasonable value
    for mass in mass_ints:
        model = model_grab(mass)
        model.plot_HR(phases=available_phases, phasecolor=colours)
    legend_elements = []
    for num, colour in enumerate(colours):
        legend_elements.append(Patch(facecolor=colour, edgecolor='black', label=phase_names[num]))
    plt.legend(handles=legend_elements, loc='lower left')
    if xlim != False:
        if ylim != False:
            plt.xlim(xlim), plt.ylim(ylim)
            plt.gca().invert_xaxis()
    plt.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
    plt.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
    plt.show()
#hr_multiplot([35,36,37,38,39,40], False, False)

# Burn fraction plot. This is generalised but currently set for phases.
# xparam is the x axis argument
# yparam is a list of y-axis parameters (i.e. mass fracts, log'l's, etc)
# xlabel = xaxis label
# ylabels are the PATCH labels for them
# axes_scales = False unless you know what the hell it does
# logx/logy = true or false: if your values are already logarithmic, don't change to true.
# loc = legend location
# finaly = ylabel for the y axis (not patches, just axis.)
# titler = title, obviously.
def burnplot(mass_int, xlim, ylim, xparam, yparam, xlabel, ylabels, axes_scales, logx, logy, loc, finaly, titler):
    model = model_grab(mass_int)
    model_xvalues = model.eeps[xparam]
    model_yvalues = [model.eeps[p] for p in yparam]
    ############## TEST #####################

    fig, axs = plt.subplots(1, constrained_layout=True, dpi=600)

    if xlim != False:
        axs.set(xlim=xlim)
    if ylim != False:
        axs.set(ylim=ylim)

    # Legend elements for the actual burns
    legend_elements = []
    for num, yvals in enumerate(model_yvalues):
        axs.plot(model_xvalues, yvals, color=colours[num], lw=1)
        legend_elements.append(Patch(facecolor=colours[num], label=ylabels[num]))

    if axes_scales != False:
        plt.axis(axes_scales)

    # Identify the phases present, their age limits, and create patch list.
    phase_strings = [("{0:.0f}").format(d) for d in model.eeps['phase']]
    phase_strings.append("99")  # Random limiter for loops that operate for == statements.

    # First grab patches for all present phases.
    present_phases = []
    for num, phase in enumerate(available_phases_str):
        if phase in phase_strings:
            present_phases.append(phase)

    # Then estimate the ranges for each phase that is present, going left-to-right.
    # You can use find() to do this. This works, so we won't clean it up/change it, but find() would make it easier.
    phase_ranges = []
    for phase in present_phases:
        start, finish = None, None
        # Find start of phase
        for num, d in enumerate(phase_strings):
            if d == phase:
                start = num
                break
        # Find end of phase.
        for num, d in enumerate(phase_strings):
            if phase_strings[num + start] != phase:
                finish = num + start - 1
                break
        phase_ranges.append([start, finish])
    phase_ranges_times = [model.eeps['star_age'][d] for d in phase_ranges]

    # Estimate the length of time spent in each phase.]
    phase_times = [model.eeps['star_age'][d[1]] - model.eeps['star_age'][d[0]] for d in phase_ranges]
    phase_lengths_str = [("{:.5e}").format(d) for d in phase_times]

    # Now we iterate over the present phases and set up the necessary patches for plotting. Also sort out the spans.
    patch_elements_phases = []
    patch_elements_ages = []
    for num, phase in enumerate(available_phases_str):
        for index, value in enumerate(present_phases):
            if value == phase:
                patch_elements_phases.append(
                    Patch(facecolor=phase_colours_extralinear[num], label=("{0}").format(phase_names[num]), alpha=0.5))
                patch_elements_ages.append(
                    Patch(facecolor=phase_colours_extralinear[num], label=phase_lengths_str[index], alpha=0.5))
                axs.axvspan(phase_ranges_times[index][0], phase_ranges_times[index][1], facecolor=phase_colours_extralinear[num], alpha=0.5)

    axs.set(xlabel=xlabel,
            ylabel=finaly,
            title=titler)


    ax2 = axs.twinx()
    ax2.get_yaxis().set_visible(False)


    # Shrink current axis by X%
    shrinker = 0.9
    box = axs.get_position()
    axs.set_position([box.x0, box.y0, box.width * shrinker, box.height])
    box2 = ax2.get_position()
    ax2.set_position([box2.x0, box2.y0, box2.width * shrinker, box2.height])

    axs.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 0.48), title="Fusion Process")
    ax2.legend(handles=patch_elements_phases,loc='lower left', bbox_to_anchor=(1, 0.48), title="Phases")



    plt.show()
def massplot(mass_int, xlim, ylim, xparam, yparam, xlabel, ylabels, axes_scales, logx, logy, loc, finaly, titler):
    model = model_grab(mass_int)
    model_xvalues = model.eeps[xparam]
    model_yvalues = [model.eeps[p] for p in yparam]
    ############## TEST #####################

    fig, axs = plt.subplots(1, constrained_layout=True, dpi=600)

    if xlim != False:
        axs.set(xlim=xlim)
    if ylim != False:
        axs.set(ylim=ylim)

    # Legend elements for the actual burns
    legend_elements = []
    for num, yvals in enumerate(model_yvalues):
        axs.plot(model_xvalues, yvals, color=colours[num], lw=1)
        legend_elements.append(Patch(facecolor=colours[num], label=ylabels[num]))

    if axes_scales != False:
        plt.axis(axes_scales)

    # Identify the phases present, their age limits, and create patch list.
    phase_strings = [("{0:.0f}").format(d) for d in model.eeps['phase']]
    phase_strings.append("99")  # Random limiter for loops that operate for == statements.

    # First grab patches for all present phases.
    present_phases = []
    for num, phase in enumerate(available_phases_str):
        if phase in phase_strings:
            present_phases.append(phase)

    # Then estimate the ranges for each phase that is present, going left-to-right.
    # You can use find() to do this. This works, so we won't clean it up/change it, but find() would make it easier.
    phase_ranges = []
    for phase in present_phases:
        start, finish = None, None
        # Find start of phase
        for num, d in enumerate(phase_strings):
            if d == phase:
                start = num
                break
        # Find end of phase.
        for num, d in enumerate(phase_strings):
            if phase_strings[num + start] != phase:
                finish = num + start - 1
                break
        phase_ranges.append([start, finish])
    phase_ranges_times = [model.eeps['star_age'][d] for d in phase_ranges]

    # Estimate the length of time spent in each phase.]
    phase_times = [model.eeps['star_age'][d[1]] - model.eeps['star_age'][d[0]] for d in phase_ranges]
    phase_lengths_str = [("{:.5e}").format(d) for d in phase_times]

    # Now we iterate over the present phases and set up the necessary patches for plotting. Also sort out the spans.
    patch_elements_phases = []
    patch_elements_ages = []
    for num, phase in enumerate(available_phases_str):
        for index, value in enumerate(present_phases):
            if value == phase:
                patch_elements_phases.append(
                    Patch(facecolor=phase_colours_extralinear[num], label=("{0}").format(phase_names[num]),
                          alpha=0.5))
                patch_elements_ages.append(
                    Patch(facecolor=phase_colours_extralinear[num], label=phase_lengths_str[index], alpha=0.5))
                axs.axvspan(phase_ranges_times[index][0], phase_ranges_times[index][1],
                            facecolor=phase_colours_extralinear[num], alpha=0.5)

    axs.set(xlabel=xlabel,
            ylabel=finaly,
            title=titler)

    ax2 = axs.twinx()
    ax2.get_yaxis().set_visible(False)

    # Shrink current axis by X%
    shrinker = 0.9
    box = axs.get_position()
    axs.set_position([box.x0, box.y0, box.width * shrinker, box.height])
    box2 = ax2.get_position()
    ax2.set_position([box2.x0, box2.y0, box2.width * shrinker, box2.height])

    axs.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 0.48))
    ax2.legend(handles=patch_elements_phases, loc='lower left', bbox_to_anchor=(1, 0.48), title="Phases")
    fig.set_size_inches(7.5,4)
    plt.savefig("Massplot_ending.png")
    plt.show()
#burnplot(1, False, [-50,5], 'star_age', ["pp", "cno", "tri_alfa", "burn_c", "burn_n", "burn_o"], "Age / Yr", ['P-P', "CNO", r'$3\alpha$', "C", "N", "O"],False, True, False, 'lower right', r'$\log{\frac{L}{L_X}}$', "             Logarithmic contribution to Luminosity for fusion processes, " + r'1 M$_\odot$') # [1.13e10,1.15e10] 1 sol mass
#massplot(1, False, [-50,5], 'star_age', ["pp", "cno", "tri_alfa", "burn_c", "burn_n", "burn_o"], "Age / Yr", ['P-P', "CNO", r'$3\alpha$', "C", "N", "O"],False, True, False, 'lower right', r'$\log{\frac{L}{L_X}}$', "             Logarithmic contribution to Luminosity for fusion processes, " + r'1 M$_\odot$') # [1.13e10,1.15e10] 1 sol mass
#massplot(1, [1.146e10,1.1465e10], [0.45,1], 'star_age', ["star_mass"], "Age / Yr", [r'$\frac{M}{M_{\odot}}$'],False, True, False, 'lower right', r'$\frac{M}{M_{\odot}}$', "Total stellar mass") # [1.13e10,1.15e10] 1 sol mass


# , r'$\Sigma$' + " H", r'$\Sigma$' + "He", r'$\Sigma$' + " Z" , "log_LH", "log_LHe", "log_LZ"

# This one will instead plot other metrics of convection/mixing
# Max grad T div grad A is ratio of temperature gradient (ambient) to adiabatic gradient
# If this exceeds 1, then convection is possible.
# Minimum ratio of gas pressure to the total pressure in the stellar interior is minpgasdivp
def convplot_alternative(mass_int, xlim, ylim, title, loc):
    # Set up model for mass_int + grab relevant parameters (hard coded)
    model = model_grab(mass_int)
    conv_labels = ['log_surf_z',
                   'min_Pgas_div_P',
                   'max_gradT_div_grada',
                   'phase',
                   'star_age']
    parameters = [model.eeps[d] for d in conv_labels]

    # Set up the multiplotter w/ plots.
    fig, axes = plt.subplots(1, dpi=450, sharex="col", constrained_layout=True)
    axes.plot(parameters[-1], parameters[0], color='black', lw=1)
    if title != False:
        fig.suptitle(title)
    if xlim != False:
        axes.set(xlim=xlim)
    if ylim != False:
        axes.set(ylim=ylim)

    phase_delimiters = []
    # Identify the phase ranges
    for phaseval in available_phases:
        for num, val in enumerate(parameters[-2]):
            print(phaseval, val)
            if (val+0.1) >= phaseval:
                phase_delimiters.append(num)
                print("yep")
                break
    phase_delimiters.append(len(parameters[-2])-1)
    phase_len = len(phase_delimiters) - 1
    age_delimiters = [parameters[-1][d] for d in phase_delimiters]
    print(age_delimiters)

    # Produce the spans necessary
    for val in range(phase_len):
        axes.axvspan(age_delimiters[val], age_delimiters[val + 1], facecolor=phase_colours_extralinear[val], alpha=0.5)


    legend_elements = []
    for num, colour in enumerate(phase_colours_extralinear):
        legend_elements.append(Patch(facecolor=colour, edgecolor='black', label=phase_names[num]))
    plt.legend(handles=legend_elements, loc=loc)

    plt.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
    plt.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
    axes.set(xlabel="Age / Yr",
             title="Surface Abundance Ratio",
             ylabel="logZ / dex")

    fig.set_size_inches(15,5)
    plt.savefig("metz.png")
    plt.show()
#convplot_alternative(1, [1.146e10,1.1465e10], False, False, "lower right")


# Generic plotter for solo SEP project
def sep_plot(mass_int, xlim, ylim, titles, save):
    model = model_grab(mass_int)
    model_x = model.eeps['star_age']

    # Identify the phases present, their age limits, and create patch list.
    phase_strings = [("{0:.0f}").format(d) for d in model.eeps['phase']]
    phase_strings.append("99")  # Random limiter for loops that operate for == statements.

    # First grab patches for all present phases.
    present_phases = []
    for num, phase in enumerate(available_phases_str):
        if phase in phase_strings:
            present_phases.append(phase)

    # Then estimate the ranges for each phase that is present, going left-to-right.
    # You can use find() to do this. This works, so we won't clean it up/change it, but find() would make it easier.
    # Now we iterate over the present phases and set up the necessary patches for plotting.

    phase_ranges = []
    for phase in present_phases:
        start, finish = None, None
        # Find start of phase
        for num, d in enumerate(phase_strings):
            if d == phase:
                start = num
                break
        # Find end of phase.
        for num, d in enumerate(phase_strings):
            if phase_strings[num + start] != phase:
                finish = num + start - 1
                break
        phase_ranges.append([start, finish])
    phase_ranges_times = [model.eeps['star_age'][d] for d in phase_ranges]



    # Estimate the length of time spent in each phase.]
    phase_times = [model.eeps['star_age'][d[1]] - model.eeps['star_age'][d[0]] for d in phase_ranges]
    phase_lengths_str = [("{:.5e}").format(d) for d in phase_times]


    # Now we iterate over the present phases and set up the necessary patches for plotting.
    patch_elements_phases = []
    patch_elements_ages = []
    for num,phase in enumerate(available_phases_str):
        for index, value in enumerate(present_phases):
            if value == phase:
                patch_elements_phases.append(Patch(facecolor=phase_colours_extralinear[num],label=("{0}").format(phase_names[num]), alpha=0.5))
                patch_elements_ages.append(Patch(alpha=0.5, facecolor=phase_colours_extralinear[num],label=phase_lengths_str[index]))



    CoreT = model.eeps['log_center_T']
    H, HE, Z = model.eeps['log_LH'],model.eeps['log_LHe'],model.eeps['log_LZ']
    fig, axs = plt.subplots(9, sharex="col",  constrained_layout=True, dpi=300)
    axs[-1].set(xlabel="Age / yr")
    if xlim != False:
        for axis in axs:
            axis.set(xlim=xlim)
    if ylim != False:
        for num,axis in enumerate(axs):
            axis.set(ylim=ylim[num])
    if titles != False:
        for num,axis in enumerate(axs):
            axis.text(.005, .92, titles[num],
                    horizontalalignment='left',
                    transform=axis.transAxes)
    fig.set_size_inches(7,7)
    axs[0].plot(model_x, CoreT, lw=1, label=r'$\log{T_{core}}$', color='black')
    # Legend elements for the actual burns
    shrinker = 0.8
    bboxta = (1.3,0.48)
    box = axs[0].get_position()
    axs[0].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    axs[0].legend(handles=patch_elements_phases,loc='center right', title="Phase",bbox_to_anchor=bboxta)
    axs[0].set(ylabel=r'$\log{T_{core}}$')



    axs[1].plot(model_x, H, lw=0.5, label="H", color="black")
    axs[1].plot(model_x, HE, lw=0.5, label="He", color="lime")
    axs[1].plot(model_x, Z, lw=0.5, label="Z", color='red')
    # Legend elements for the actual burns
    box = axs[1].get_position()
    axs[1].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    axs[1].legend(loc='center right', bbox_to_anchor=bboxta, title="Fusion Type")
    axs[1].set(ylabel=r'$\log\frac{L}{L_x}$')


    # BURN PHASES
    for axis in axs:
        # Now we iterate over the present phases and set up the necessary patches for plotting. Also sort out the spans.
        for num, phase in enumerate(available_phases_str):
            for index, value in enumerate(present_phases):
                if value == phase:
                    axis.axvspan(phase_ranges_times[index][0], phase_ranges_times[index][1],
                                facecolor=phase_colours_extralinear[num], alpha=0.5)

    yparam = ["pp", "cno", "tri_alfa", "burn_c", "burn_n", "burn_o"] # ,
    ylabels = ['P-P', "CNO", r'$3\alpha$', "C", "N", "O"]
    model_yvalues = [model.eeps[p] for p in yparam]


    # Legend elements for the actual burns
    legend_elements_two = []
    for num, yvals in enumerate(model_yvalues):
        axs[2].plot(model_x, yvals, color=colours[num], lw=1)
        legend_elements_two.append(Patch(facecolor=colours[num], label=ylabels[num]))
    box = axs[2].get_position()
    axs[2].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    axs[2].legend(handles=legend_elements_two, loc='center right', bbox_to_anchor=bboxta, title="Fusion Type")
    axs[2].set(ylabel=r'$\log\frac{L}{L_x}$')


    # axs[3] is the mass plots!!!
    mass_parameters = ['center_h1','center_he4','center_c12','center_n14','center_o16','center_ne20','center_mg24','center_si28']
    mass_parameters_labels = [r'${^1}H$',r'${^4}He$',r'$^{12}C$',r'$^{14}N$',r'$^{16}O$',r'$^{20}Ne$',r'$^{24}Mg$',r'$^{28}Si$']
    ymasses = [model.eeps[d] for d in mass_parameters]
    legend_elements_three = []
    for num, yvals in enumerate(ymasses):
        axs[3].plot(model_x, yvals, color=colours[num], lw=1)
        legend_elements_three.append(Patch(facecolor=colours[num], label=mass_parameters_labels[num]))
    box = axs[3].get_position()
    axs[3].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    axs[3].legend(handles=legend_elements_three, loc='center right', bbox_to_anchor=bboxta + (0,-0.8), title="Core Frac")
    axs[3].set(ylabel=r'$\frac{M}{M_x}$')

    # axs[4] is the pressure plots!!!
    mass_parameters = ['center_degeneracy']
    mass_parameters_labels = [r'$\mu_{e-}$' + " / " + r'$k_b{T}$']
    ymasses = [model.eeps[d] for d in mass_parameters]
    legend_elements_three = []
    for num, yvals in enumerate(ymasses):
        axs[4].plot(model_x, yvals, color=colours[num], lw=1)
        legend_elements_three.append(Patch(facecolor=colours[num], label=mass_parameters_labels[num]))
    box = axs[4].get_position()
    axs[4].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    #axs[4].legend(handles=legend_elements_three, loc='center right', bbox_to_anchor=bboxta, title="Core Frac")
    axs[4].set(ylabel=r'$\mu_{e}$' + " / " + r'$k_b{T}$')

    # axs[4] is the pressure plots!!!
    mass_parameters = ['mass_conv_core']
    mass_parameters_labels = [r'$\frac{M_{Conc}}{M_\odot$']
    ymasses = [model.eeps[d] for d in mass_parameters]
    legend_elements_three = []
    for num, yvals in enumerate(ymasses):
        axs[5].plot(model_x, yvals, color=colours[num], lw=1)
        legend_elements_three.append(Patch(facecolor=colours[num], label=mass_parameters_labels[num]))
    box = axs[5].get_position()
    axs[5].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    #axs[4].legend(handles=legend_elements_three, loc='center right', bbox_to_anchor=bboxta, title="Core Frac")
    axs[5].set(ylabel=r'$\frac{M_{conv}}{M_{\odot}}$')
    axs[5].legend(handles=patch_elements_phases,loc='center right', title="Phase",bbox_to_anchor=bboxta)

    # axs[4] is the pressure plots!!!
    mass_parameters = ['log_R']
    mass_parameters_labels = [r'$\frac{R}{R_{\odot}}$']
    ymasses = [model.eeps[d] for d in mass_parameters]
    legend_elements_three = []
    for num, yvals in enumerate(ymasses):
        axs[6].plot(model_x, yvals, color=colours[num], lw=1)
        legend_elements_three.append(Patch(facecolor=colours[num], label=mass_parameters_labels[num]))
    box = axs[6].get_position()
    axs[6].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    #axs[4].legend(handles=legend_elements_three, loc='center right', bbox_to_anchor=bboxta, title="Core Frac")
    axs[6].set(ylabel = r'$\frac{R}{R_{\odot}}$')

    # axs[4] is the pressure plots!!!
    mass_parameters = ['star_mass']
    mass_parameters_labels = [r'$\frac{M}{M_{\odot}}$']
    ymasses = [model.eeps[d] for d in mass_parameters]
    legend_elements_three = []
    for num, yvals in enumerate(ymasses):
        axs[7].plot(model_x, yvals, color=colours[num], lw=1)
        legend_elements_three.append(Patch(facecolor=colours[num], label=mass_parameters_labels[num]))
    box = axs[7].get_position()
    axs[7].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    #axs[4].legend(handles=legend_elements_three, loc='center right', bbox_to_anchor=bboxta, title="Core Frac")
    axs[7].set(ylabel = r'$\frac{M}{M_{\odot}}$')

    # axs[4] is the pressure plots!!!
    mass_parameters = ['log_L']
    mass_parameters_labels = [r'$\frac{L}{L_{\odot}}$']
    ymasses = [model.eeps[d] for d in mass_parameters]
    legend_elements_three = []
    for num, yvals in enumerate(ymasses):
        axs[8].plot(model_x, yvals, color=colours[num], lw=1)
        legend_elements_three.append(Patch(facecolor=colours[num], label=mass_parameters_labels[num]))
    box = axs[8].get_position()
    axs[8].set_position([box.x0, box.y0, box.width * shrinker, box.height])
    #axs[4].legend(handles=legend_elements_three, loc='center right', bbox_to_anchor=bboxta, title="Core Frac")
    axs[8].set(ylabel = r'$\log\frac{L}{L_{\odot}}$')

    fig.set_size_inches(10, 16)
    fig.savefig(save)
    plt.show()
sep_plot(1, [1.14615e10, 1.14655e10], [[5,8.4],[-10,5],[-10,5],[0,1.1],[0,50],[0,1.1], [0,2.5], [0,1.1],[-0.2,5]], ["Core Temperature", "Bulk log burn fraction", "Process-specific log burn fraction", "Elemental core mass fraction","Central e- chemical potential", "Convective core mass", "Log Radius", "Mass", "Log Luminosity"], "zoom3.png") # [1.144e10,1.1465e10]


