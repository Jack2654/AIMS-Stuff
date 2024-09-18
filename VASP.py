import matplotlib.pyplot as plt

def plot_DOS(input):
    with open(input, "r") as f:
        lines = f.readlines()
    energies = []
    dos = []
    for count, line in enumerate(lines[6:1006]):
        temp = line.split()
        energies.append(float(temp[0]))
        dos.append(float(temp[1]))

    plt.fill_between(energies, 0, dos)
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS (/eV)")
    plt.show()


plot_DOS("../../Documents/24-25/Classes/ME490/tmp/DOSCAR")
