"""
gravity_hydro_bridge
Created on 02-12-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
from __future__ import print_function
import numpy
import os
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.community.bhtree.interface import Bhtree
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.ext.evrard_test import uniform_unit_sphere
from initialize_apep import Initialize_inner_binary, Initialize_apep
from amuse.units.constants import G
from amuse.io import write_set_to_file, read_set_from_file
from amuse.community.gadget2.interface import Gadget2
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH

def new_sph_particles_from_stellar_wind(stars, mgas):
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

        if si.name == "WC8":
            vrot = 530 | units.kms
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
            # print(add[ri].velocity.length().in_(units.kms))
        new_sph.add_particles(add)
    return new_sph


def inject_supernova_energy(gas_particles,
                            explosion_energy=1.0e+51|units.erg,
                            exploding_region=10|units.RSun):
    inner = gas_particles.select(
        lambda pos: pos.length_squared() < exploding_region**2,
        ["position"])
    print(len(inner), "innermost particles selected.")
    print("Adding", explosion_energy / inner.total_mass(), "of supernova " \
        "(specific internal) energy to each of the n=", len(inner), "SPH particles.")
    inner.u += explosion_energy / inner.total_mass()
    return inner


def gravity_hydro_bridge():
    stars = Initialize_apep()
    binary = stars[0:2]
    #print(binary)
    a = stars.position.length().amax()

    dt = (
        0.01 | units.yr
    )  # changed from 1 yr because the SPH particles were shooting away in that time
    binary_mass_loss_rate = abs(binary.dmdt.sum() * dt)
    mgas = 0.01 * binary_mass_loss_rate

    stars.h_smooth = 0.0 * a
    stars.u = 0 | units.kms**2

    converter = nbody_system.nbody_to_si(stars.mass.sum(), a)
    # gravity = Bhtree(converter)
    gravity = ph4(converter, redirection="none")
    gravity.particles.add_particles(stars)
    gravity.parameters.epsilon_squared = (10 | units.RSun) ** 2

    channel = {
        "from apep:": stars.new_channel_to(gravity.particles),
        "to_apep": gravity.particles.new_channel_to(stars),
    }
    # channel_from_gravity = gravity.particles.new_channel_to(stars)
    # channel_from_to_gravity = stars.new_channel_to(gravity.particles)

    ism = Particles(0)
    ism.mass = mgas
    ism.position = (0, 0, 0) | units.AU
    ism.velocity = (0, 0, 0) | units.kms
    ism.u = 0 | units.kms**2
    ism.h_smooth = 0.01 * a
    # hydro = Gadget2(converter)
    hydro = Fi(converter, redirection="none")
    hydro.parameters.timestep = dt / 8.0  # smaller_dt for hydro is better
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.radiation_flag = False
    hydro.parameters.self_gravity_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.gamma = 1.0
    hydro.parameters.isothermal_flag = True
    hydro.parameters.epsilon_squared = (10 | units.RSun) ** 2

    hydro.parameters.timestep = 1 | units.s  # Steven's patch
    print(hydro.model_time)
    hydro.evolve_model(hydro.model_time)
    hydro.parameters.timestep = dt
    hydro.evolve_model(0 | units.yr)
    if len(ism) > 0:
        hydro.gas_particles.add_particles(ism)
    # hydro.parameters.periodic_box_size = 10000 * a

    channel.update({"from_wind": ism.new_channel_to(hydro.gas_particles)})
    channel.update({"to_wind": hydro.gas_particles.new_channel_to(ism)})
    channel.update({"from_stars": binary.new_channel_to(hydro.dm_particles)})
    channel.update({"to_stars": hydro.dm_particles.new_channel_to(binary)})

    # channel_from_hydro = hydro.gas_particles.new_channel_to(ism)
    # channel_from_to_hydro = ism.new_channel_to(hydro.gas_particles)

    moving_bodies = ParticlesSuperset([stars, ism])
    #print(moving_bodies[1])
    model_time = 0 | units.yr
    filename = "snewstellargravhydroTEST.hdf5"
    if len(ism) > 0:
        write_set_to_file(moving_bodies, filename, "hdf5")

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,), False)
    gravhydro.timestep = min(dt, 2 * hydro.parameters.timestep)

    istep = 0
    save_every = 1
    first_time = True
    once_supernova = True
    time_for_supernova = .8 | units.yr
    pickle_file = './WN_structure.pkl'

    while model_time < 150 | units.yr:
        model_time += gravhydro.timestep
        binary.Mwind += binary.dmdt * dt
        print("Wind mass loss: ", binary.Mwind)

        # First time we must guarantee some sph particles otherwise it crashes
        if first_time:
            mass_sph = 0.1 * binary_mass_loss_rate

        if once_supernova:
            new_sph = new_sph_particles_from_stellar_wind(binary, mass_sph) #this line needs improvement as one of the stars of the binary is gone now and a non wind emitting core is there
        #else:
         #   new_sph = new_sph_particles_from_stellar_wind(binary[0], mass_sph)
        if first_time:
            mass_sph = mgas
            first_time = False

        if model_time>=time_for_supernova and once_supernova:
            once_supernova = False
            supernova_model = convert_stellar_model_to_SPH(None, #assuming WN goes supernova
        1000, #these sph particles will make the supernova
        seed=12345,pickle_file=pickle_file,
        with_core_particle=True,
        target_core_mass = 1.2|units.MSun # WN was 11 solar mass so 1.2 solar mass core should be a good guess
    )
            core, supernova_gas, core_radius = supernova_model.core_particle, supernova_model.gas_particles, supernova_model.core_radius
            print('*'*100)
            print("SUPERNOVA")
            print('*' * 100)
            print('Core Radius',core_radius)
            supernova_gas= inject_supernova_energy(supernova_gas, exploding_region=1 | units.RSun)
            if len(supernova_gas) > 0:
                ism.add_particles(supernova_gas)
                moving_bodies.remove_particle(particle=moving_bodies[1])
                stars.add_particle(particle=core)

        print("Total number of particles", len(ism))
        if len(new_sph) > 0:
            ism.add_particles(new_sph)
            print(ism[3].velocity.length().in_(units.kms))
            ism.synchronize_to(hydro.gas_particles)
        gravhydro.evolve_model(model_time)
        channel["to_apep"].copy()  # channel_from_gravity.copy()
        channel["to_wind"].copy()  # channel_from_hydro.copy()
        channel["to_stars"].copy()
        channel["to_wind"].copy_attributes(["u"])  # channel_from_hydro.copy_attributes(["u"])
        print(ism[3].velocity.length().in_(units.kms))

        if istep % 1 / save_every == 0:
            filename = f"stellargravhydro_supernova_{int(istep/save_every)}.hdf5"
            if os.path.exists(filename):
                os.remove(filename)
            write_set_to_file(moving_bodies, filename, "hdf5")
        istep += 1

    gravity.stop()
    hydro.stop()


if __name__ in ("__main__", "__plot__"):
    gravity_hydro_bridge()


# if __name__ in ('__main__', '__plot__'):
#     gravity_hydro_bridge()
