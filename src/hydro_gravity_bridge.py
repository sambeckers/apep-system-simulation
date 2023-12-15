"""
hydro_sph
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
import numpy
from amuse.lab import *
from amuse import datamodel
from amuse.io import write_set_to_file
from amuse.units import units, nbody_system
from amuse.lab import Particles, ParticlesSuperset
from amuse.units.constants import G
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.community.seba.interface import SeBa
from amuse.community.fi.interface import Fi
import os
from amuse.community.ph4.interface import ph4
from amuse.couple import bridge
from amuse.ext.composition_methods import *

# Own modules
from initialize_apep import Initialize_apep
from initialize_apep import M_loss_WN, M_loss_WC, v_inf_wind_WN, v_inf_wind_WC, P_binary, d_WR_binary_to_SG
# from plotting_routine import plot_sph_particles


set_printing_strategy(
    "custom",  # nbody_converter = converter,
    preferred_units=[units.MSun, units.AU, units.Myr],
    precision=5,
    prefix="",
    separator=" [",
    suffix="]",
)


def new_sph_particles_from_stellar_wind(
    stars, mgas
):  # was in the example codes and also in the textbook
    new_sph = datamodel.Particles(0)
    for si in stars:
        Ngas = int(si.Mwind / mgas)
        print(Ngas)
        if Ngas == 0:
            continue
        add = datamodel.Particles(Ngas)
        add.mass = mgas
        add.h_smooth = 0.0 | units.parsec

        dx, dy, dz = uniform_unit_sphere(Ngas).make_xyz()
        add.x = si.x + (dx * si.radius)
        add.y = si.y + (dy * si.radius)
        add.z = si.z + (dz * si.radius)
        for ri in range(len(add)):
            r = add[ri].position - si.position
            r = r / r.length()
            v_wind = (G * si.mass / (add[ri].position - si.position).length()).sqrt()
            add.u = 0.5 * (v_wind) ** 2
            add.vx = si.vx + r[0] * si.terminal_wind_velocity
            add.vy = si.vy + r[1] * si.terminal_wind_velocity
            add.vz = si.vz + r[2] * si.terminal_wind_velocity
        new_sph.add_particles(add)
    return new_sph


def main():
    apep = Initialize_apep()
    inner_binary = Particles(particles=[apep[apep.name=="WC8"], apep[apep.name=="WN46b"]])

    dt = 2 | units.day
    mgas = 0.1 * abs(apep.dmdt.sum() * dt)  # mass of gas lost through stellar wind

    # Setting up gravity
    converter = nbody_system.nbody_to_si(apep.mass.sum(), d_WR_binary_to_SG)
    gravity = ph4(converter)
    gravity.particles.add_particles(apep)
    channel = {"from apep:": apep.new_channel_to(gravity.particles),
                "to_apep": gravity.particles.new_channel_to(apep)}

    # Setting up the SPH particles
    wind = Particles(0)
    wind.mass = mgas
    wind.position = (0, 0, 0) | units.AU
    wind.velocity = (0, 0, 0) | units.kms
    wind.u = 0 | units.m**2 * units.s**-2
    wind.h_smooth = 0.01 * d_WR_binary_to_SG

    # Setting up the hydrodynamics
    hydro = Fi(converter, redirection="none", mode="openmp")
    if len(wind) > 0:
        hydro.gas_particles.add_particles(wind)
        hydro.dm_particles.add_particles(inner_binary.as_set())
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.timestep = dt
    hydro.parameters.periodic_box_size = 1000 * d_WR_binary_to_SG
    hydro_to_framework = hydro.gas_particles.new_channel_to(wind)
    
    channel.update({"from_wind": wind.new_channel_to(hydro.gas_particles)})
    channel.update({"to_wind": hydro.gas_particles.new_channel_to(wind)})
    channel.update({"from_stars": inner_binary.new_channel_to(hydro.dm_particles)})
    channel.update({"to_stars": hydro.dm_particles.new_channel_to(inner_binary)})

    # Setting up the bridge between gravity and hydrodynamics
    gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,))
    gravhydro.timestep = 0.2*P_binary

    # moving_wind = ParticlesSuperset([apep, wind])

    # Evolving the system
    filename = "hydro_outflow.hdf5"
    istep = 0
    while (
        hydro.model_time < 200 | units.day
    ):  # evolving for 2 days just to see if this works
        apep.Mwind += apep.dmdt * dt
        new_sph = new_sph_particles_from_stellar_wind(inner_binary, mgas)

        if len(new_sph) > 0:
            wind.add_particles(new_sph)
            wind.synchronize_to(hydro.gas_particles)
        print("time=", hydro.model_time, "Ngas=", len(wind), mgas * len(wind))
        if len(wind) > 100:
            gravhydro.evolve_model(hydro.model_time + dt)
            hydro_to_framework.copy()
            channel["to_apep"].copy()
            channel["to_wind"].copy()
            channel["to_stars"].copy()
            
            # Saving the stellar wind and stars separately at each timestep
            if istep % 1 == 0:
                filename_h = f"hydro_outflow_step_{istep}.hdf5"  
                if os.path.exists(filename_h):
                    os.remove(filename_h)
                write_set_to_file(hydro.gas_particles, filename, "hdf5", append_to_file=False)

                filename_g = f"gravity_apep_step_{istep}.hdf5"
                if os.path.exists(filename_g):
                    os.remove(filename_g)
                write_set_to_file(apep, filename_g, "hdf5", append_to_file=False)

            istep += 1
    gravity.stop()
    hydro.stop()

if __name__ in ("__main__", "__plot__"):
    main()
