"""
hydro_gravity_bridge
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
import os
from amuse.lab import *
from amuse import datamodel
from amuse.io import write_set_to_file
from amuse.units import units, nbody_system
from amuse.lab import Particles, ParticlesSuperset
from amuse.units.constants import G
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.community.fi.interface import Fi
from amuse.community.ph4.interface import ph4
from amuse.couple import bridge

# Own modules
from initialize_apep import Initialize_apep, Initialize_inner_binary
from initialize_apep import d_WR_binary_to_SG

def new_sph_particles_from_stellar_wind(stars, mgas):
    """ 
    Creates new SPH particles that represent the stellar wind of the stars in the binary system.
    Args:
        stars (Particles): the WR binary
        mgas (Particles): mass of the gas lost through stellar wind in the WR binary

    Returns:
        new_sph (Particles): the new SPH particles representing the stellar wind
    """
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


def gravity_hydro_bridge():
    # Setting up the system
    apep = Initialize_apep()
    inner_binary = Initialize_inner_binary()
    dt = .1 | units.yr
    mgas = 0.1 * abs(inner_binary.dmdt.sum() * dt)  # mass of gas lost through stellar wind

    apep.h_smooth = 0.0 * d_WR_binary_to_SG
    apep.u = 0.0 | units.kms**2

    # Setting up gravity
    converter = nbody_system.nbody_to_si(apep.mass.sum(), d_WR_binary_to_SG)
    gravity = ph4(converter)
    gravity.particles.add_particles(apep)
    gravity.paramaters.epsilon_squared = (10|units.RSun)**2
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
    hydro.parameters.radiation_flag = False
    hydro.parameters.self_gravity_flag = True
    hydro.parameters.gamma = 1
    hydro.parameters.isothermal_flag = True
    hydro.parameters.epsilon_squared = (10|units.RSun)**2
    # hydro.parameters.periodic_box_size = 1000 * d_WR_binary_to_SG

    # Solution to broken Fi for larger timesteps (credit Steven Rieder):
    hydro.parameters.timestep = 1 | units.s
    hydro.evolve_model(hydro.model_time)
    hydro.parameters.timestep = dt
    hydro.evolve_model(0 | units.yr)
    if len(wind) > 0:
        hydro.gas_particles.add_particles(wind)

    channel.update({"from_wind": wind.new_channel_to(hydro.gas_particles)})
    channel.update({"to_wind": hydro.gas_particles.new_channel_to(wind)})
    channel.update({"from_stars": inner_binary.new_channel_to(hydro.dm_particles)})
    channel.update({"to_stars": hydro.dm_particles.new_channel_to(inner_binary)})

    # Setting up the bridge between gravity and hydrodynamics
    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,))
    gravhydro.timestep = min(dt, 2 * hydro.parameters.timestep)

    # Evolving the system
    moving_bodies = ParticlesSuperset([apep, wind])

    istep = 0
    save_every = 1

    while model_time < 150 | units.yr:
        model_time += dt

        inner_binary.Mwind += inner_binary.dmdt * dt
        new_sph = new_sph_particles_from_stellar_wind(inner_binary, mgas)
        if len(new_sph) > 0:
            wind.add_particles(new_sph)
            wind.synchronize_to(hydro.gas_particles)
        
        gravhydro.evolve_model(model_time)
        channel["to_apep"].copy()
        channel["to_wind"].copy()
        channel["to_stars"].copy()
        channel["to_wind"].copy_attributes(["u"])

        if istep % 1 / save_every == 0:
            filename = f"apep_{int(istep/save_every)}.hdf5"
            if os.path.exists(filename):
                os.remove(filename)
            write_set_to_file(moving_bodies, filename, "hdf5")
        istep += 1

    gravity.stop()
    hydro.stop()

if __name__ in ("__main__", "__plot__"):
    gravity_hydro_bridge()
