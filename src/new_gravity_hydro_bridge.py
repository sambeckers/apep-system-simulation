"""
gravity_hydro_bridge
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
from __future__ import print_function
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.community.bhtree.interface import Bhtree
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.ext.evrard_test import uniform_unit_sphere
from initialize_apep import Initialize_inner_binary, Initialize_apep
from amuse.units.constants import G
from amuse.io import write_set_to_file
from amuse.community.gadget2.interface import Gadget2

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

def gravity_hydro_bridge():
    stars = Initialize_apep()
    binary = stars[0:2]
    print(binary)
    a = stars.position.length().amax()
    dt = 1 | units.yr
    mgas = .1 * abs(binary.dmdt.sum() * dt)

    stars.h_smooth = 0.0 * a
    stars.u = 0 | units.kms ** 2

    converter = nbody_system.nbody_to_si(stars.mass.sum(), a)
    #gravity = Bhtree(converter)
    gravity = ph4(converter, redirection="none")
    gravity.particles.add_particles(stars)
    gravity.parameters.epsilon_squared = (10 | units.RSun) ** 2

    channel_from_gravity = gravity.particles.new_channel_to(stars)
    channel_from_to_gravity = stars.new_channel_to(gravity.particles)

    ism = Particles(0)
    ism.mass = mgas
    ism.position = (0, 0, 0) | units.AU
    ism.velocity = (0, 0, 0) | units.kms
    ism.u = 0 | units.kms ** 2
    ism.h_smooth = 0.01 * a
    #hydro = Gadget2(converter)
    hydro = Fi(converter, redirection="none")
    #hydro.parameters.timestep = dt
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.radiation_flag = False
    hydro.parameters.self_gravity_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.gamma = 1.
    hydro.parameters.isothermal_flag = True
    hydro.parameters.epsilon_squared = (10 | units.RSun) ** 2

    hydro.parameters.timestep = 1 | units.s #Steven's patch
    print(hydro.model_time)
    hydro.evolve_model(hydro.model_time )
    hydro.parameters.timestep = dt
    hydro.evolve_model(0 |units.yr)
    if len(ism) > 0:
        hydro.gas_particles.add_particles(ism)
    hydro.parameters.periodic_box_size = 10000 * a

    channel_from_hydro = hydro.gas_particles.new_channel_to(ism)
    channel_from_to_hydro = ism.new_channel_to(hydro.gas_particles)

    moving_bodies = ParticlesSuperset([stars, ism])
    model_time = 0 | units.yr
    filename = "snewstellargravhydro.hdf5"
    if len(ism) > 0:
        write_set_to_file(moving_bodies, filename, 'hdf5')

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,))
    gravhydro.timestep = dt #min(dt, 2 * hydro.parameters.timestep)

    istep = 0
    save_every = 1

    while (model_time < 150 | units.yr):
        model_time += dt
        stars.Mwind += stars.dmdt * dt
        new_sph = new_sph_particles_from_stellar_wind(binary, mgas)
        if len(new_sph) > 0:
            ism.add_particles(new_sph)
            ism.synchronize_to(hydro.gas_particles)
        gravhydro.evolve_model(model_time)
        channel_from_gravity.copy()
        channel_from_hydro.copy()
        channel_from_hydro.copy_attributes(["u"])

        if istep % 1/save_every == 0:
            filename = f"snewstellargravhydro_{int(istep/save_every)}.hdf5"
            write_set_to_file(moving_bodies, filename, 'hdf5')
        istep += 1

    gravity.stop()
    hydro.stop()



if __name__ in ('__main__', '__plot__'):
    gravity_hydro_bridge()



# if __name__ in ('__main__', '__plot__'):
#     gravity_hydro_bridge()
