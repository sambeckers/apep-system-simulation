"""
hydro_sph
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
import numpy
from amuse.lab import *
from amuse import datamodel
from amuse.units import units
from amuse.lab import Particles
from amuse.units.constants import G
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.community.seba.interface import SeBa

# Own modules
from initialize_apep import Initialize_inner_binary


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
        Mgas = mgas * Ngas
        si.Mwind += Mgas
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
    inner_binary = (
        Initialize_inner_binary()
    )  # Carbon star first. Moved to center of mass.

    a = inner_binary.position.length().amax()
    vc = G * inner_binary.mass.sum() / a
    stellar = SeBa()
    stellar.particles.add_particles(inner_binary)
    stellar_to_framework = stellar.particles.new_channel_to(inner_binary)
    stellar.evolve_model(1 | units.Myr)
    stellar_to_framework.copy_attributes(["mass", "radius", "temperature"])
    dt = 0.1 | units.Myr
    stellar.evolve_model(
        (1 | units.Myr) + dt
    )  # evolving for a very short time just to see if it works (will evolve both the WR stars and supergiant separately in future)
    inner_binary[0].dmdt = (
        10 ** (-4.3) | units.MSun / units.yr
    )  # mass loss rate took from the ppt we made
    inner_binary[1].dmdt = (
        10 ** (-4.5) | units.MSun / units.yr
    )  # mass loss rate took from the ppt we made
    # stars.dmdt = (stellar.particles.mass-stars.mass)/dt
    inner_binary.Mwind = 0 | units.MSun
    inner_binary[0].terminal_wind_velocity = (
        3500 | units.kms
    )  # wind_velocity took from the ppt
    inner_binary[1].terminal_wind_velocity = (
        2100 | units.kms
    )  # wind_velocity took from the ppt
    stellar.stop()
    dt = 0.1 | units.day
    mgas = 0.1 * abs(
        inner_binary.dmdt.sum() * dt
    )  # mass of gas lost through stellar wind

    converter = nbody_system.nbody_to_si(1 | units.MSun, a)
    bodies = Particles(0)
    bodies.mass = mgas
    bodies.position = (0, 0, 0) | units.AU
    bodies.velocity = (0, 0, 0) | units.kms
    bodies.u = 0 | units.m**2 * units.s**-2
    bodies.h_smooth = 0.01 * a

    hydro = Fi(converter, redirection="none")
    if len(bodies) > 0:
        hydro.gas_particles.add_particles(bodies)
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.timestep = dt
    hydro.parameters.periodic_box_size = 1000 * a
    hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)

    moving_bodies = ParticlesSuperset([inner_binary, bodies])
    filename = "hydro_outflow.hdf5"
    istep = 0
    while (
        hydro.model_time < 2 | units.day
    ):  # evolving for 2 days just to see if this works
        inner_binary.Mwind += inner_binary.dmdt * dt
        new_sph = new_sph_particles_from_stellar_wind(inner_binary, mgas)

        if len(new_sph) > 0:
            bodies.add_particles(new_sph)
            bodies.synchronize_to(hydro.gas_particles)
        print("time=", hydro.model_time, "Ngas=", len(bodies), mgas * len(bodies))
        if len(bodies) > 100:
            hydro.evolve_model(hydro.model_time + dt)
            hydro_to_framework.copy()
            if istep % 1 == 0:
                filename = f"hydro_outflow_step_{istep}.hdf5"  # saving the system as new hdf5 file at each step
                write_set_to_file(
                    hydro.gas_particles, filename, "hdf5", append_to_file=False
                )

            istep += 1
    hydro.stop()


if __name__ in ("__main__", "__plot__"):
    main()
