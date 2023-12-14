"""
hydro_sph
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""

from __future__ import print_function
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.ext.evrard_test import uniform_unit_sphere
from initialize_apep import Initialize_inner_binary, Initialize_apep
from amuse.units.constants import G
from amuse.io import write_set_to_file

def new_sph_particles_from_stellar_wind(stars, mgas):
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

def gravity_hydro_bridge(a, ecc, t_end, n_steps, Rgas, Mgas, Ngas):
    stars = Initialize_inner_binary()
    a = stars.position.length().amax()
    dt = 0.1 | units.day
    mgas = 0.1 * abs(stars.dmdt.sum() * dt)

    stars.h_smooth = 0.0 * a
    stars.u = 0 | units.kms ** 2

    converter = nbody_system.nbody_to_si(stars.mass.sum(), a)
    gravity = ph4(converter, redirection="none")
    gravity.particles.add_particles(stars)
    gravity.parameters.epsilon_squared = (10 | units.RSun) ** 2

    channel_from_gravity = gravity.particles.new_channel_to(stars)
    channel_from_to_gravity = stars.new_channel_to(gravity.particles)

    ism = Particles(0)
    ism.mass = mgas
    ism.position = (0, 0, 0) | units.AU
    ism.velocity = (0, 0, 0) | units.kms
    ism.u = 0 | units.m ** 2 * units.s ** -2
    ism.h_smooth = 0.01 * a

    hydro = Fi(converter, redirection="none")
    hydro.parameters.timestep = dt / 8.
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.radiation_flag = False
    hydro.parameters.self_gravity_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.gamma = 1.
    hydro.parameters.isothermal_flag = True
    hydro.parameters.epsilon_squared = (10 | units.RSun) ** 2
    if len(ism) > 0:
        hydro.gas_particles.add_particles(ism)
    hydro.parameters.periodic_box_size = 10000 * a

    channel_from_hydro = hydro.gas_particles.new_channel_to(ism)
    channel_from_to_hydro = ism.new_channel_to(hydro.gas_particles)

    moving_bodies = ParticlesSuperset([stars, ism])
    model_time = 0 | units.day
    filename = "newstellargravhydro.hdf5"
    if len(ism) > 0:
        write_set_to_file(moving_bodies, filename, 'hdf5')

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,))
    gravhydro.timestep = min(dt, 2 * hydro.parameters.timestep)

    istep = 0
    while model_time < t_end:
        model_time += dt
        stars.Mwind += stars.dmdt * dt
        new_sph = new_sph_particles_from_stellar_wind(stars, mgas)
        if len(new_sph) > 0:
            ism.add_particles(new_sph)
            ism.synchronize_to(hydro.gas_particles)
        gravhydro.evolve_model(model_time)
        channel_from_gravity.copy()
        channel_from_hydro.copy()
        channel_from_hydro.copy_attributes(["u"])

        if istep % 1 == 0:
            filename = f"newstellargravhydro_{istep}.hdf5"
            write_set_to_file(moving_bodies, filename, 'hdf5')
        istep += 1

    gravity.stop()
    hydro.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 1000,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-N", dest="Ngas", type="int", default = 1024,
                      help="number of gas particles [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mgas", type="float", default = 1|units.MSun,
                      help="Mass of the gas [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rgas", type="float", default = 1|units.AU,
                      help="Size of the gas distribution [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="a", type="float", default = 0.2|units.AU,
                      help="initial orbital separation [%default]")
    result.add_option("-e", dest="ecc", type="float", default = 0.0,
                      help="initial orbital eccentricity [%default]")
    result.add_option("-t", unit=units.yr,
                      dest="t_end", type="float", default = 20|units.day,
                      help="end time of the simulation [%default]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments = new_option_parser().parse_args()
    gravity_hydro_bridge(**o.__dict__)
