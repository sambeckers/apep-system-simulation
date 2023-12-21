"""
plotting_routine
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie

Several plotting routines to plot the Apep system/SPH particles in 2D/3D.
"""
# Importing modules
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import glob
from amuse.plot import _plot
from amuse.io import read_set_from_file
from amuse.units import units
import numpy as np

# Setting frame size
frame_size = 5000

def plot_apep(hdf5_files, timestep_supernova, plot_supernova = False, min_size=100, max_size=10000):
    """Plots the Apep system in 2D at each timestep, with a black background. 
    You can choose to plot the supernova at a certain timestep.

    Args:
        hdf5_files (str) : list of hdf5 files for each timestep
        timestep_supernova (int, optional): timestep at which supernova occurs
        plot_supernova (bool, optional): whether to plot the supernova or not. Defaults to False.
        min_size (int, optional): minimum size of SPH particles, in pixel**2. Defaults to 100.
        max_size (int, optional): maximum size of SPH particles, in pixel**2. Defaults to 10000.
    """
    with plt.rc_context({'axes.edgecolor': 'white',
                         'xtick.color': 'white',
                         'ytick.color': 'white',
                         'figure.facecolor': 'black',
                         'axes.facecolor': 'black',
                         'axes.labelcolor': 'white',
                         'axes.titlecolor': 'white'}):
        for i, filename in enumerate(hdf5_files):
            print(filename)

            #Load in the particles and separate the stellar wind particles from the stars
            particles = read_set_from_file(filename, "hdf5", copy_history=False, close_file=False)
            stellar_wind = particles[3:]
            apep = particles[0:3]

            positions = stellar_wind.position
            h_smooths = stellar_wind.h_smooth
            x, y = positions.x, positions.y

            pos = apep.position
            x_a, y_a = pos.x, pos.y

            # Set up the plot
            plt.figure(dpi=450)
            n_pixels = plt.gcf().get_dpi() * plt.gcf().get_size_inches()
            current_axes = plt.gca()
            current_axes.set_aspect("equal", adjustable="datalim")
            phys_to_pix2 = n_pixels[0] * n_pixels[1] / ((max(x) - min(x)) ** 2 + (max(y) - min(y)) ** 2)
            sizes = np.minimum(np.maximum((h_smooths ** 2 * phys_to_pix2), min_size), max_size)

            plt.scatter(x.value_in(units.AU), y.value_in(units.AU), s=10, c='orange', alpha=0.05) #plotting the stellar wind particles, s=10 but can be set to 'sizes'
            if plot_supernova:
                if i<timestep_supernova-1:
                    for x_i, y_i, color in zip(x_a, y_a, ['blue', 'cyan', 'red']):
                        plt.scatter(x_i.value_in(units.AU), y_i.value_in(units.AU), s=100, c=color, marker='*')
            else:
                for x_i, y_i, color in zip(x_a, y_a, ['blue', 'red', 'cyan']):
                    plt.scatter(x_i.value_in(units.AU), y_i.value_in(units.AU), s=100, c=color, marker='*')
            plt.xlabel('x [AU]')
            plt.ylabel('y [AU]')
            plt.text(0.1, 0.94, f'  timestep = {i} ',
                     horizontalalignment='center',
                     verticalalignment='center',
                     transform=current_axes.transAxes, color='white', fontsize=12) #adding the timestep to the plot

            # Adding custom legend:
            for i in plt.legend(handles=[plt.scatter([], [], marker=".", color='orange', label='Stellar wind'),
                                        plt.scatter([], [], marker="*", color='blue', label='WC8'),
                                        plt.scatter([], [], marker="*", color='cyan', label='WN46b'),
                                        plt.scatter([], [], marker="*", color='red', label='O8')],
                                loc='upper center', bbox_to_anchor=(0.5, 1.05),
                                ncol=2, fancybox=True, shadow=True).get_texts():
                i.set_color("white")
            plt.tight_layout()
            plt.xlim(-frame_size, frame_size)
            plt.ylim(-frame_size, frame_size)
            plt.savefig(filename.replace('.hdf5', '.png'))
            plt.close()

def plot_sph_particles_2D(filename):
    """Default plotting routine to plot the SPH particles in 2D, adapted from AMUSE plotting routine.

    Args:
        filename (str) : hdf5 file to be plotted
    """
    particles = read_set_from_file(filename, "hdf5", copy_history=False,
                                   close_file=False)  # reading the saved hdf5 files
    _plot.sph_particles_plot(particles)  # plotting routine to plot the SPH particles
    plt.savefig(filename.replace('.hdf5', '.png'))

def smart_length_units_for_vector_quantity(quantity):
    """
    Given a vector quantity, return a length unit that is appropriate for. Dependency for plot_sph_particles_3D.
    Adapted from AMUSE plotting routine.
    """
    length_units = [units.Mpc, units.kpc, units.parsec, units.AU, units.RSun, units.km]
    total_size = max(quantity) - min(quantity)
    for length_unit in length_units:
        if total_size > (1 | length_unit):
            return length_unit
    return units.m

def plot_sph_particles_3D(filename, u_range=None, min_size=100, max_size=10000,
                          alpha=0.1, gd_particles=None, width=None, view=None):
    """
    Very simple and fast procedure to make a plot of the hydrodynamics state of
    a set of SPH particles. The particles must have the following attributes defined:
    position, u, h_smooth. Adapted from AMUSE plotting routine.

    Args: 
        particles: the SPH particles to be plotted
        u_range: range of internal energy for color scale [umin, umax]
        min_size: minimum size to use for plotting particles, in pixel**2
        max_size: maximum size to use for plotting particles, in pixel**2
        alpha: the opacity of each particle
        gd_particles: non-SPH particles can be indicated with white circles
        view: the (physical) region to plot [xmin, xmax, ymin, ymax]
    """
    particles = read_set_from_file(filename, "hdf5", copy_history=False, close_file=False)
    positions = particles.position
    us = particles.u
    h_smooths = particles.h_smooth
    x, y, z = positions.x, positions.y, positions.z
    z, x, y, us, h_smooths = z.sorted_with(x, y, us, h_smooths)

    if u_range:
        u_min, u_max = u_range
    else:
        u_min, u_max = min(us), max(us)
    log_u = np.log((us / u_min)) / np.log((u_max / u_min))
    clipped_log_u = np.minimum(np.ones_like(log_u), np.maximum(np.zeros_like(log_u), log_u))

    red = 1.0 - clipped_log_u ** 4
    blue = clipped_log_u ** 4
    green = np.minimum(red, blue)

    colors = np.transpose(np.array([red, green, blue]))
    n_pixels = plt.gcf().get_dpi() * plt.gcf().get_size_inches()

    current_axes = plt.gca()
    try:
        current_axes.set_facecolor('#101010')
    except:
        current_axes.set_axis_bgcolor('#101010')
    if width is not None:
        view = width * [-0.5, 0.5, -0.5, 0.5]

    if view:
        current_axes.set_aspect("equal", adjustable="box")
        length_unit = smart_length_units_for_vector_quantity(view)
        current_axes.set_xlim(view[0].value_in(length_unit),
                              view[1].value_in(length_unit), emit=True, auto=False)
        current_axes.set_ylim(view[2].value_in(length_unit),
                              view[3].value_in(length_unit), emit=True, auto=False)
        phys_to_pix2 = n_pixels[0] * n_pixels[1] / ((view[1] - view[0]) ** 2 + (view[3] - view[2]) ** 2)
    else:
        current_axes.set_aspect("equal", adjustable="datalim")
        length_unit = smart_length_units_for_vector_quantity(x)
        phys_to_pix2 = n_pixels[0] * n_pixels[1] / ((max(x) - min(x)) ** 2 + (max(y) - min(y)) ** 2)
    sizes = np.minimum(np.maximum((h_smooths ** 2 * phys_to_pix2), min_size), max_size)

    x = x.as_quantity_in(length_unit)
    y = y.as_quantity_in(length_unit)
    z = z.as_quantity_in(length_unit)

    fig = plt.figure(dpi=100)
    ax = plt.axes(projection='3d')
    ax.scatter3D(x.value_in(units.AU), y.value_in(units.AU), z.value_in(units.AU), s=sizes, c=colors, edgecolors="none",
                 alpha=alpha)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.tight_layout()
    # plt.show()
    plt.savefig(filename.replace('.hdf5', '.png'))
    plt.close()
    # if gd_particles:
    #     scatter(gd_particles.x, gd_particles.y, c='w', marker='o')

def main():
    hdf5_files = sorted(glob.glob("stellargravhydro_supernova_*.hdf5"), key=lambda x: int(x.split('_')[-1].split('.')[0]))
    # for filename in hdf5_files:
    #     plot_sph_particles_2D(filename)
    #     plot_sph_particles_3D(filename)
    plot_apep(hdf5_files,timestep_supernova = 180,plot_supernova=True)


if __name__ == "__main__":
    main()
