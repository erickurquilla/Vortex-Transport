import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np
import h5py
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator,LogLocator)
import matplotlib as mpl
import matplotlib.ticker
from matplotlib import cm
import matplotlib.tri as tri
from matplotlib.animation import FuncAnimation
import matplotlib.tri as tri

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

################################################################################
# PDF PLOTS
################################################################################

# Get directory names
dir_names = glob.glob("output/step_*")
dir_names = sorted(dir_names, key=lambda x: int(x.lstrip("output/step_")))

# Loop through each frame
for i in range(len(dir_names)):

    print(dir_names[i])

    # Load data for the current frame
    file_names = glob.glob(dir_names[i] + "/element_*.txt")
    # dir_names = sorted(file_names, key=lambda x: int(x.lstrip(dir_names[i] + "/element_").rstrip(".txt")))
    file_names = sorted(file_names, key=lambda x: int(x.split("/element_")[1].rstrip(".txt")))

    x, y, u1, u2, u3, u4, t = [], [], [], [], [], [], []

    for dire in file_names:
        df = pd.read_csv(dire, sep='\s+')
        x.append(df['x'])
        y.append(df['y'])
        u1.append(df['u0'])
        u2.append(df['u1'])
        u3.append(df['u2'])
        u4.append(df['u3'])
        t.append(list(df['time'])[0])

    # Flatten the data
    x = np.array(x).flatten()
    y = np.array(y).flatten()
    rho = np.array(u1).flatten()
    u = np.array(u2).flatten() / rho
    v = np.array(u3).flatten() / rho
    E = np.array(u4).flatten() / rho
    gamma = 1.4
    p = rho * (gamma - 1) * (E - (u**2 + v**2) / 2)
    H = E + p / rho

    # Create the triangular mesh
    triangles = tri.Triangulation(x, y)

    # Create subplots with shared axes
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12, 10))

    # Plot each variable in a separate subplot
    im1 = axs[0, 0].tripcolor(triangles, rho, cmap='viridis')
    axs[0, 0].set_title(r'$\rho$')
    fig.colorbar(im1, ax=axs[0, 0])

    im2 = axs[0, 1].tripcolor(triangles, u, cmap='viridis')
    axs[0, 1].set_title(r'$u$')
    fig.colorbar(im2, ax=axs[0, 1])

    im3 = axs[1, 0].tripcolor(triangles, p, cmap='viridis')
    axs[1, 0].set_title(r'$p$')
    fig.colorbar(im3, ax=axs[1, 0])

    im4 = axs[1, 1].tripcolor(triangles, v, cmap='viridis')
    axs[1, 1].set_title(r'$v$')
    fig.colorbar(im4, ax=axs[1, 1])

    plt.suptitle(f'Time: {t[0]}')

    axs[0, 0].set_ylim([-5, 5])
    axs[0, 0].set_xlim([-5, 5])
    axs[1, 0].set_ylim([-5, 5])
    axs[1, 0].set_xlim([-5, 5])
    axs[0, 1].set_ylim([-5, 5])
    axs[0, 1].set_xlim([-5, 5])
    axs[1, 1].set_ylim([-5, 5])
    axs[1, 1].set_xlim([-5, 5])

    axs[0, 0].yaxis.set_major_locator(MultipleLocator(2))
    axs[0, 0].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 0].xaxis.set_major_locator(MultipleLocator(2))
    axs[0, 0].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 0].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    axs[0, 1].yaxis.set_major_locator(MultipleLocator(2))
    axs[0, 1].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 1].xaxis.set_major_locator(MultipleLocator(2))
    axs[0, 1].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 1].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    axs[1, 0].yaxis.set_major_locator(MultipleLocator(2))
    axs[1, 0].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 0].xaxis.set_major_locator(MultipleLocator(2))
    axs[1, 0].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 0].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    axs[1, 1].yaxis.set_major_locator(MultipleLocator(2))
    axs[1, 1].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 1].xaxis.set_major_locator(MultipleLocator(2))
    axs[1, 1].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 1].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    plt.suptitle(f'$t$ = {t[0]}')

    axs[0, 0].set_ylabel(r'$y$')
    axs[1, 0].set_xlabel(r'$x$')
    axs[1, 0].set_ylabel(r'$y$')
    axs[1, 1].set_xlabel(r'$x$')

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot
    plt.savefig('plottools/step_' + dir_names[i].lstrip("output/step_") + '.pdf', bbox_inches='tight')

    # Close the plot to avoid memory leaks
    plt.close(fig)
    plt.clf()

################################################################################
# MP4 ANIMATION
################################################################################

# Initialize variables for color bars
cbar1 = None
cbar2 = None
cbar3 = None
cbar4 = None

# Get directory names
dir_names = glob.glob("output/step_*")
dir_names = sorted(dir_names, key=lambda x: int(x.lstrip("output/step_")))

# Create a function to update the plot for each frame
def update_4fig(frame):

    global cbar1, cbar2, cbar3, cbar4

    print(f'frame: {frame}')

    # Load data for the current frame
    file_names = glob.glob(dir_names[frame] + "/element_*.txt")
    file_names = sorted(file_names, key=lambda x: int(x.split("/element_")[1].rstrip(".txt")))

    x, y, u1, u2, u3, u4, t = [], [], [], [], [], [], []

    for dire in file_names:
        df = pd.read_csv(dire, sep='\s+')
        x.append(df['x'])
        y.append(df['y'])
        u1.append(df['u0'])
        u2.append(df['u1'])
        u3.append(df['u2'])
        u4.append(df['u3'])
        t.append(list(df['time'])[0])

    # Flatten the data
    x = np.array(x).flatten()
    y = np.array(y).flatten()
    rho = np.array(u1).flatten()
    u = np.array(u2).flatten() / rho
    v = np.array(u3).flatten() / rho
    E = np.array(u4).flatten() / rho
    gamma = 1.4
    p = rho * (gamma - 1) * (E - (u**2 + v**2) / 2)
    H = E + p / rho

    # Clear the axes
    for ax in axs.flat:
        ax.clear()

    # Create the triangular mesh
    triangles = tri.Triangulation(x, y)

    # Plot each variable in a separate subplot
    im1 = axs[0, 0].tripcolor(triangles, rho, cmap='viridis')
    axs[0, 0].set_title(r'$\rho$')
    if cbar1 is None:
        cbar1 = fig.colorbar(im1, ax=axs[0, 0])
    else:
        cbar1.update_normal(im1)

    im2 = axs[0, 1].tripcolor(triangles, u, cmap='viridis')
    axs[0, 1].set_title(r'$u$')
    if cbar2 is None:
        cbar2 = fig.colorbar(im2, ax=axs[0, 1])
    else:
        cbar2.update_normal(im2)

    im3 = axs[1, 0].tripcolor(triangles, p, cmap='viridis')
    axs[1, 0].set_title(r'$p$')
    if cbar3 is None:
        cbar3 = fig.colorbar(im3, ax=axs[1, 0])
    else:
        cbar3.update_normal(im3)

    im4 = axs[1, 1].tripcolor(triangles, v, cmap='viridis')
    axs[1, 1].set_title(r'$v$')
    if cbar4 is None:
        cbar4 = fig.colorbar(im4, ax=axs[1, 1])
    else:
        cbar4.update_normal(im4)

    axs[0, 0].set_ylim([-5, 5])
    axs[0, 0].set_xlim([-5, 5])
    axs[1, 0].set_ylim([-5, 5])
    axs[1, 0].set_xlim([-5, 5])
    axs[0, 1].set_ylim([-5, 5])
    axs[0, 1].set_xlim([-5, 5])
    axs[1, 1].set_ylim([-5, 5])
    axs[1, 1].set_xlim([-5, 5])

    axs[0, 0].yaxis.set_major_locator(MultipleLocator(2))
    axs[0, 0].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 0].xaxis.set_major_locator(MultipleLocator(2))
    axs[0, 0].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 0].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    axs[0, 1].yaxis.set_major_locator(MultipleLocator(2))
    axs[0, 1].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 1].xaxis.set_major_locator(MultipleLocator(2))
    axs[0, 1].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[0, 1].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    axs[1, 0].yaxis.set_major_locator(MultipleLocator(2))
    axs[1, 0].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 0].xaxis.set_major_locator(MultipleLocator(2))
    axs[1, 0].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 0].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    axs[1, 1].yaxis.set_major_locator(MultipleLocator(2))
    axs[1, 1].yaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 1].xaxis.set_major_locator(MultipleLocator(2))
    axs[1, 1].xaxis.set_minor_locator(MultipleLocator(2/10))
    axs[1, 1].tick_params(axis='both', which='both', direction='in', right=True, top=True)

    plt.suptitle(f'$t$ = {t[0]}')

    axs[0, 0].set_ylabel(r'$y$')
    axs[1, 0].set_xlabel(r'$x$')
    axs[1, 0].set_ylabel(r'$y$')
    axs[1, 1].set_xlabel(r'$x$')

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot
    # plt.savefig('plottools/step_' + dir_names[frame].lstrip("output/step_") + '.pdf', bbox_inches='tight')

# Create the figure and axes outside the functions
fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12, 10))

# Create the animation using the specified figure and update function
ani = FuncAnimation(fig, update_4fig, frames=len(dir_names), interval=400)

# Save the animation as a video file
ani.save('plottools/animation.mp4', writer='ffmpeg', dpi=200)

# Close the plot to avoid memory leaks
plt.close(fig)