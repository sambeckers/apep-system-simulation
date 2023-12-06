"""
plotting_routine
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
import matplotlib.pyplot as plt
import glob
from amuse.plot import _plot
from amuse.io import read_set_from_file

def plot_sph_particles(filename):
    #particles = read_set_from_file(filename, "hdf5")
    particles = read_set_from_file(filename, "hdf5", copy_history = False, close_file = False) #reading the saved hdf5 files
    _plot.sph_particles_plot(particles) #plotting routine to plot the SPH particles
    plt.savefig(filename.replace('.hdf5', '.png'))
    plt.close()

def main():
    hdf5_files = glob.glob("hydro_outflow_step_*.hdf5")
    for filename in hdf5_files:
        plot_sph_particles(filename)

if __name__ == "__main__":
    main()
