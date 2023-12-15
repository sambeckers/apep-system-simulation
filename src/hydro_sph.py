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

# Own modules
from initialize_apep import Initialize_inner_binary
from initialize_apep import M_loss_WN, M_loss_WC, v_inf_wind_WN, v_inf_wind_WC, 
from plotting_routine import plot_sph_particles


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
    inner_binary = (
        Initialize_inner_binary()
    )  # Carbon star first. Moved to center of mass.

    a = inner_binary.position.length().amax()

    dt = 3 | units.day
    mgas = 100 * abs(
        inner_binary.dmdt.sum() * dt  # 0.01 |units.MEarth
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

    # moving_bodies = ParticlesSuperset([inner_binary, bodies])
    filename = "hydro_outflow.hdf5"
    istep = 0
    while (
        hydro.model_time < 10 | units.yr
    ):  # evolving for 2 days just to see if this works
        inner_binary.Mwind += inner_binary.dmdt * dt
        new_sph = new_sph_particles_from_stellar_wind(inner_binary, mgas)

        if len(new_sph) > 0:
            bodies.add_particles(new_sph)
            bodies.synchronize_to(hydro.gas_particles)
        print(
            "time=",
            hydro.model_time.in_(units.yr),
            "Ngas=",
            len(bodies),
            mgas * len(bodies),
        )
        if len(bodies) > 100:
            hydro.evolve_model(hydro.model_time + dt)
            hydro_to_framework.copy()
            savestep = 5
            if istep % savestep == 0:
                filename = f"hydro_outflow_step_{int(istep/savestep)}.hdf5"  # saving the system as new hdf5 file at each step
                if os.path.exists(filename):
                    os.remove(filename)
                write_set_to_file(
                    hydro.gas_particles, filename, "hdf5", append_to_file=False
                )
                plot_sph_particles(filename)
                if os.path.exists(filename):
                    os.remove(filename)

            istep += 1
    hydro.stop()


if __name__ in ("__main__", "__plot__"):
    main()
