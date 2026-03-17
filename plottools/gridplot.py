import glob
import matplotlib.pyplot as plt
import numpy as np
import h5py
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator,LogLocator)
import matplotlib as mpl
import matplotlib.ticker
from matplotlib import cm

mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['axes.linewidth'] = 2

def read_input_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                key, value = line.strip().split('=')
                if key.startswith('coordinate'):
                    data.append(float(value))
                elif key in ['element_number', 'type', 'right_element', 'left_element', 'vertical_element']:
                    data.append(int(value))
    return data

if __name__ == "__main__":

    # Get all file names matching the pattern "element*.txt"
    file_names = glob.glob("grid/element_*.txt")

    # Remove the ".txt" extension from each file name and sort them
    sorted_files = sorted(file_names, key=lambda x: int(x.lstrip("grid/element_").rstrip(".txt")))

    fig, ax = plt.subplots()

    for dire in sorted_files:
        data = read_input_file(dire)
        x = [data[1],data[3],data[5],data[1]]
        y = [data[2],data[4],data[6],data[2]]
        ax.plot(x,y)
    fig.savefig('plottools/grid.pdf',bbox_inches='tight')
    plt.clf()