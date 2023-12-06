from amuse.ext.stellar_wind import new_stellar_wind
from amuse.units import units, nbody_system
from amuse.lab import Particles, ParticlesSuperset
from amuse.community.fi.interface import Fi
import os
from amuse.io import write_set_to_file

# Own modules
from initialize_apep import Initialize_inner_binary
from initialize_apep import M_loss_WN, M_loss_WC, v_inf_wind_WN, v_inf_wind_WC
from plotting_routine import plot_sph_particles


def main():
    inner_binary = (
        Initialize_inner_binary()
    )  # Carbon star first. Moved to center of mass.

    a = inner_binary.position.length().amax()

    dt = 0.1 | units.day
    mgas = 0.1 * abs(
        inner_binary.dmdt.sum() * dt
    )  # mass of each gas lost through stellar wind

    converter = nbody_system.nbody_to_si(1 | units.MSun, a)
    bodies = Particles(0)
    bodies.mass = mgas
    bodies.position = (0, 0, 0) | units.AU
    bodies.velocity = (0, 0, 0) | units.kms
    bodies.u = 0 | units.m**2 * units.s**-2
    bodies.h_smooth = 0.01 * a

    # wind = new_stellar_wind(mgas, target_gas=bodies, timestep=dt, derived_forces=True)
    # wind.particles.add_particles(inner_binary)
    # channel_to_wind = inner_binary.new_channel_to(wind.particles)

    hydro = Fi(converter, redirection="none")
    if len(bodies) > 0:
        hydro.gas_particles.add_particles(bodies)
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.timestep = dt
    hydro.parameters.periodic_box_size = 1000 * a
    hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)

    filename = "hydro_outflow.hdf5"
    istep = 0
    while (
        hydro.model_time < 20 | units.day
    ):  # evolving for 2 days just to see if this works
        # inner_binary.Mwind += inner_binary.dmdt * dt
        new_sph = new_stellar_wind(
            mgas, target_gas=bodies, timestep=dt, derive_from_evolution=True
        )
        # print(new_sph.internal_energy_from_velocity)
        # print(new_sph.particles)
        # new_sph.particles.add_particles(inner_binary)
        # setattr(
        #     new_sph.particles,
        #     "u",
        #     [0 | units.m**2 * units.s**-2 for i in range(len(new_sph.particles))],
        # )
        # print(new_sph.particles.u)

        if len(new_sph.particles) > 0:
            inner_binary.u = 0 | units.m**2 * units.s**-2
            new_sph.particles.add_particles(inner_binary)
            bodies.add_particles(new_sph.particles)
            bodies.synchronize_to(hydro.gas_particles)
        print("time=", hydro.model_time, "Ngas=", len(bodies), mgas * len(bodies))
        if len(bodies) > 100:
            hydro.evolve_model(hydro.model_time + dt)
            hydro_to_framework.copy()
            if istep % 1 == 0:
                filename = f"hydro_outflow_step_{istep}.hdf5"  # saving the system as new hdf5 file at each step
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


if __name__ == "__main__":
    main()
