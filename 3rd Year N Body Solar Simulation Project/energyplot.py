import numpy as np
import matplotlib.pyplot as plt

def energyplot(energytension, timestep):
    energyfile = open(energytension)
    energylines = energyfile.readlines()
    energies = []
    for line in energylines:
        linesplit = line.split()
        energies.append(float(linesplit[0]))

    steplist = np.arange(0, len(energies), 1)
    timelist = [d*float(timestep) for d in steplist]

    fig = plt.figure(figsize=(10, 10))
    ax69 = fig.add_subplot(111)
    ax69.scatter(timelist, energies, label="Total Energy vs. Time for Simulation", lw=0.000001)

    ax69.set_xlabel("Time Elapsed, Days")
    ax69.set_ylabel("Total Energy, Joules")

    plt.legend()
    ax69.grid()

    plt.savefig("totenergy.png", dpi=600)
    plt.show()


def paramgetter(paramname_and_extension):
    owo = open(paramname_and_extension)
    # Timestep, Total Steps, Graviconstant
    uwu = owo.readlines()
    step = uwu[0].split(",")[0]
    owo.close()
    return step

energyplot("Energy.dat", paramgetter("parameters.dat"))

