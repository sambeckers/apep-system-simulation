"""
gravity_hydro_bridge
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

- Creates stellar wind from SPH particles based on WR binary properties
- Includes a Bridge between the hydro and gravity codes.
- Includes rotation of the WC8 star
- Evolves the model and saves selected timesteps to a hdf5 file
"""
# Importing modules
from __future__ import print_function
import os
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.units.constants import G
from amuse.io import write_set_to_file
from initialize_apep import Initialize_apep, d_WR_binary_to_SG, P_binary # own module


def new_sph_particles_from_stellar_wind(stars, mgas):
    """ 
    Creates new SPH particles that represent the stellar wind of the stars in the binary system.
    Adapted from AMUSE example codes and textbook.
    Args:
        stars (Particles): the WR binary
        mgas (Particles): mass of the gas lost through stellar wind in the WR binary

    Returns:
        new_sph (Particles): the new SPH particles representing the stellar wind
    """
    new_sph = datamodel.Particles(0)
    for si in stars:
        Ngas = int(si.Mwind / mgas)
        si.Mwind -= Ngas * mgas
        print(f"Number of new particles for {si.name}", Ngas)
        if Ngas == 0:
            print("No new particles")
            continue
        add = datamodel.Particles(Ngas)
        add.mass = mgas
        add.h_smooth = 0.0 | units.parsec

        # Adding rotational velocity to WC8 star only
        if si.name == "WC8":
            vrot = 500 | units.kms
        elif si.name == "WN46b":
            vrot = 0.0 | units.kms
        else:
            raise ValueError("No star found")

        dx, dy, dz = uniform_unit_sphere(Ngas).make_xyz()
        add.x = si.x + (dx * si.radius)
        add.y = si.y + (dy * si.radius)
        add.z = si.z + (dz * si.radius)
        for ri in range(len(add)):
            r = add[ri].position - si.position
            r = r / r.length()
            v_wind = (G * si.mass / (add[ri].position - si.position).length()).sqrt()
            add[ri].u = 0.5 * (v_wind) ** 2
            add[ri].vx = si.vx + r[0] * si.terminal_wind_velocity + r[1] * vrot
            add[ri].vy = si.vy + r[1] * si.terminal_wind_velocity + r[0] * vrot
            add[ri].vz = si.vz + r[2] * si.terminal_wind_velocity
        new_sph.add_particles(add)
    
    return new_sph


def gravity_hydro_bridge():
    # Setting up the system
    apep = Initialize_apep()
    central_engine = apep[0:2]

    dt = 0.1 | units.yr
    central_engine_mass_loss_rate = abs(central_engine.dmdt.sum() * dt)
    mgas = 0.01 * central_engine_mass_loss_rate

    apep.h_smooth = 0.0 * d_WR_binary_to_SG
    apep.u = 0 | units.kms**2

    # Setting up gravity
    converter = nbody_system.nbody_to_si(apep.mass.sum(), d_WR_binary_to_SG)
    gravity = ph4(converter, redirection="none")
    gravity.particles.add_particles(apep)
    gravity.parameters.epsilon_squared = (10 | units.RSun) ** 2

    channel = {"from apep:": apep.new_channel_to(gravity.particles),
               "to_apep": gravity.particles.new_channel_to(apep),}

    # Setting up the SPH particles
    wind = Particles(0)
    wind.mass = mgas
    wind.position = (0, 0, 0) | units.AU
    wind.velocity = (0, 0, 0) | units.kms
    wind.u = 0 | units.kms**2
    wind.h_smooth = 0.01 * d_WR_binary_to_SG

    # Setting up the hydrodynamics
    hydro = Fi(converter, redirection="none")
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.timestep = dt / 8.0
    hydro.parameters.radiation_flag = False
    hydro.parameters.self_gravity_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.gamma = 1.0
    hydro.parameters.isothermal_flag = True
    hydro.parameters.epsilon_squared = (10 | units.RSun) ** 2

    # Solution to broken Fi for larger timesteps (credit Steven Rieder):
    hydro.parameters.timestep = 1 | units.s
    print(hydro.model_time)
    hydro.evolve_model(hydro.model_time)
    hydro.parameters.timestep = dt/8.
    hydro.evolve_model(0 | units.yr)
    if len(wind) > 0:
        hydro.gas_particles.add_particles(wind)

    channel.update({"from_wind": wind.new_channel_to(hydro.gas_particles)})
    channel.update({"to_wind": hydro.gas_particles.new_channel_to(wind)})
    channel.update({"from_apep": central_engine.new_channel_to(hydro.dm_particles)})
    channel.update({"to_apep": hydro.dm_particles.new_channel_to(central_engine)})

    moving_bodies = ParticlesSuperset([apep, wind])
    model_time = 0 | units.yr
    filename = "apep_rot.hdf5"
    if len(wind) > 0:
        write_set_to_file(moving_bodies, filename, "hdf5")

    # Setting up the bridge between gravity and hydrodynamics
    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,), True)
    gravhydro.timestep = min(dt, 2 * hydro.parameters.timestep)

    # Evolving the system
    istep = 0
    save_every = 1
    first_time = True

    while model_time < 150 | units.yr:
        model_time += gravhydro.timestep
        central_engine.Mwind += central_engine.dmdt * dt
        print("Wind mass loss: ", central_engine.Mwind)

        # First time we must guarantee some sph particles (otherwise it crashes)
        if first_time:
            mass_sph = 0.1 * central_engine_mass_loss_rate

        new_sph = new_sph_particles_from_stellar_wind(central_engine, mass_sph)

        if first_time:
            mass_sph = mgas
            first_time = False

        print("Total number of particles", len(wind))
        if len(new_sph) > 0:
            wind.add_particles(new_sph)
            wind.synchronize_to(hydro.gas_particles)
        
        gravhydro.evolve_model(model_time)
        channel["to_apep"].copy() 
        channel["to_wind"].copy() 
        channel["to_apep"].copy()
        channel["to_wind"].copy_attributes(["u"])

        if istep % 1 / save_every == 0:
            filename = f"apep_rot_{int(istep/save_every)}.hdf5"
            if os.path.exists(filename):
                os.remove(filename)
            write_set_to_file(moving_bodies, filename, "hdf5")
        istep += 1

    gravity.stop()
    hydro.stop()

if __name__ in ("__main__", "__plot__"):
    gravity_hydro_bridge()
