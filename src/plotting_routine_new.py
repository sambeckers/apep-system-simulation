"""
plotting_routine
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import glob
from amuse.plot import _plot
#import new_plotting as _plot
from amuse.io import read_set_from_file
from amuse.units import units
import numpy as np
import re

def plot_sph_particles_2D(filename):
    #particles = read_set_from_file(filename, "hdf5")
    particles = read_set_from_file(filename, "hdf5", copy_history = False, close_file = False) #reading the saved hdf5 files
    #print(particles)
    #_plot.plot(particles)
    xmin = -1000 | units.AU
    xmax = 1000 | units.AU
    ymin = -1000 | units.AU
    ymax = 1000 | units.AU
    _plot.sph_particles_plot(particles[3:], gd_particles=particles[0:3],  min_size = 10, max_size= 100, alpha = .05, view=[xmin, xmax, ymin, ymax]) #plotting routine to plot the SPH particles
    plt.savefig(filename.replace('.hdf5', '.png'))

def smart_length_units_for_vector_quantity(quantity):
    length_units = [units.Mpc, units.kpc, units.parsec, units.AU, units.RSun, units.km]
    total_size = max(quantity) - min(quantity)
    for length_unit in length_units:
        if total_size > (1 | length_unit):
            return length_unit
    return units.m

def plot_sph_particles_3D(filename, u_range = None, min_size = 100, max_size = 10000,
        alpha = 0.1, gd_particles=None, width=None, view=None):
    """
    Very simple and fast procedure to make a plot of the hydrodynamics state of
    a set of SPH particles. The particles must have the following attributes defined:
    position, u, h_smooth.

    :argument particles: the SPH particles to be plotted
    :argument u_range: range of internal energy for color scale [umin, umax]
    :argument min_size: minimum size to use for plotting particles, in pixel**2
    :argument max_size: maximum size to use for plotting particles, in pixel**2
    :argument alpha: the opacity of each particle
    :argument gd_particles: non-SPH particles can be indicated with white circles
    :argument view: the (physical) region to plot [xmin, xmax, ymin, ymax]
    """
    particles = read_set_from_file(filename, "hdf5", copy_history = False, close_file = False)[2:]
    positions = particles.position
    us        = particles.u
    h_smooths = particles.h_smooth
    x, y, z = positions.x, positions.y, positions.z
    z, x, y, us, h_smooths = z.sorted_with(x, y, us, h_smooths)

    if u_range:
        u_min, u_max = u_range
    else:
        u_min, u_max = min(us), max(us)
    log_u = np.log((us / u_min)) / np.log((u_max / u_min))
    clipped_log_u = np.minimum(np.ones_like(log_u), np.maximum(np.zeros_like(log_u), log_u))

    red   = 1.0 - clipped_log_u**4
    blue  = clipped_log_u**4
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
        current_axes.set_aspect("equal", adjustable = "box")
        length_unit = smart_length_units_for_vector_quantity(view)
        current_axes.set_xlim(view[0].value_in(length_unit),
            view[1].value_in(length_unit), emit=True, auto=False)
        current_axes.set_ylim(view[2].value_in(length_unit),
            view[3].value_in(length_unit), emit=True, auto=False)
        phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((view[1]-view[0])**2 + (view[3]-view[2])**2)
    else:
        current_axes.set_aspect("equal", adjustable = "datalim")
        length_unit = smart_length_units_for_vector_quantity(x)
        phys_to_pix2 = n_pixels[0]*n_pixels[1] / ((max(x)-min(x))**2 + (max(y)-min(y))**2)
    sizes = np.minimum(np.maximum((h_smooths**2 * phys_to_pix2), min_size), max_size)

    x = x.as_quantity_in(length_unit)
    y = y.as_quantity_in(length_unit)
    z = z.as_quantity_in(length_unit)

    fig = plt.figure(dpi=100)
    ax = plt.axes(projection = '3d')
    ax.scatter(x.value_in(units.AU), y.value_in(units.AU), z.value_in(units.AU), s=sizes, c=colors, edgecolors="none", alpha=alpha)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.tight_layout()
    #plt.show()
    plt.savefig(filename.replace('.hdf5', '.png'))
    plt.close()
    # if gd_particles:
    #     scatter(gd_particles.x, gd_particles.y, c='w', marker='o')

def extract_timestep_from_filename(filename):
    # Extract the timestep (last numeric part of the filename)
    match = re.search(r'(\d+)\.hdf5$', filename)
    if match:
        return int(match.group(1))
    else:
        return 53

def plot_particles(sph_particles, stars, timestep, filename):
    plt.figure(figsize=(8, 8))
    plt.scatter(sph_particles.x.value_in(units.au), sph_particles.y.value_in(units.au), c='red', marker='o', s=1.5)
    plt.scatter(stars.x.value_in(units.au), stars.y.value_in(units.au), c='yellow', marker='o', s=.8)

    plt.gca().set_facecolor('black')
    plt.title(f"timestep = {timestep}")
    plt.xlabel('X (AU)')
    plt.ylabel('Y (AU)')
    plt.xlim(-1000, 1000)
    plt.ylim(-1000, 1000)
    plt.savefig(filename.replace('.hdf5', '.png'))
    plt.close()

def main():
    hdf5_files = glob.glob("snewstellargravhydro_*.hdf5")
    #print(hdf5_files)
    for filename in hdf5_files:
        particles = read_set_from_file(filename, "hdf5", copy_history=False, close_file=False)
        #plot_sph_particles_2D(filename)
        #plot_sph_particles_3D(filename)
        timestep = extract_timestep_from_filename(filename)
        plot_particles(particles[3:],particles[0:3],timestep, filename)

if __name__ == "__main__":
    main()
