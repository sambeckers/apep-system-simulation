"""
initialize_apep
Created on 26-11-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
# Importing modules
from amuse.units import units, constants
from amuse.lab import Particles
from amuse.community.seba.interface import SeBa
from amuse.ext.orbital_elements import (
    new_binary_from_orbital_elements,
    orbital_period_to_semimajor_axis,
)

# Wolf Rayet (WR) binary
a_binary = 67 | units.au
P_binary = 125 | units.yr
ecc_binary = 0.7
inc_binary = 0 | units.deg
omega_binary = 0 | units.deg

M_enclosed_binary = 30 | units.MSun
# A guess for the individual WN46b-type WR (WN) and WC8-type WR (WC) masses:
M_WN = 11 | units.MSun
M_WC = 19 | units.MSun

R_WN = 6 | units.RSun
R_WC = 6.3 | units.RSun

T_eff_WN = 65000 | units.K
T_eff_WC = 60000 | units.K

M_loss_WN = (
    4 * 10 ** (-5) | units.MSun / units.yr
)  # mass loss rate of the WN46b-tyoe WR
M_loss_WC = (
    2.9 * 10 ** (-5) | units.MSun / units.yr
)  # mass loss rate of the WC8-type WR

# O8-Iaf Supergiant (SG)
d_WR_binary_to_SG = 1700 | units.au
P_binary_SG = 10**4 | units.yr
M_SG = 40 | units.MSun
R_SG = 35.4 | units.RSun
T_eff_SG = 30000 | units.K
M_loss_SG = 10 ** (-5) | units.MSun / units.yr

# Stellar winds
v_inf_wind_WN = 3500 | units.kms
v_inf_wind_WC = 2100 | units.kms
T_wind_WN = 0.3 * T_eff_WN
T_wind_WC = 0.3 * T_eff_WC
ecc_SG = 0  # guess
inc_SG = 25 | units.deg  # guess


def Initialize_inner_binary():
    """
    Create a Particleset of the inner binary system.
    The carbon Wolf-Rayet star is the primary,
    the nitrogen Wolf-Rayet star is the secondary.
    """
    inner_binary = new_binary_from_orbital_elements(
        mass1=M_WC,
        mass2=M_WN,
        semimajor_axis=a_binary,
        eccentricity=ecc_binary,
        inclination=inc_binary,
        true_anomaly=90 | units.deg,
    )
    # inner_binary.move_to_center()
    setattr(inner_binary, "name", ["WC8", "WN46b"])

    inner_binary[inner_binary.name == "WC8"].radius = R_WC
    inner_binary[inner_binary.name == "WN46b"].radius = R_WN

    inner_binary[inner_binary.name == "WC8"].temperature = T_eff_WC
    inner_binary[inner_binary.name == "WN46b"].temperature = T_eff_WN

    inner_binary[inner_binary.name == "WC8"].terminal_wind_velocity = v_inf_wind_WC
    inner_binary[inner_binary.name == "WN46b"].terminal_wind_velocity = v_inf_wind_WN

    inner_binary[inner_binary.name == "WC8"].wind_temperature = T_wind_WC
    inner_binary[inner_binary.name == "WN46b"].wind_temperature = T_wind_WN

    inner_binary[inner_binary.name == "WC8"].dmdt = M_loss_WC
    inner_binary[inner_binary.name == "WN46b"].dmdt = M_loss_WN
    inner_binary.Mwind = 0 | units.MSun

    # stellar = SeBa()
    # stellar.particles.add_particles(inner_binary)
    # stellar_to_framework = stellar.particles.new_channel_to(inner_binary)
    # stellar.evolve_model(1 | units.Myr)
    # stellar_to_framework.copy_attributes(["mass", "radius", "temperature"])
    # dt = 0.1 | units.Myr
    # stellar.evolve_model(
    #     (1 | units.Myr) + dt
    # )  # evolving for a very short time just to see if it works (will evolve both the WR stars and supergiant separately in future)
    # inner_binary[0].dmdt = M_loss_WC  # mass loss rate took from the ppt we made
    # inner_binary[1].dmdt = M_loss_WN  # mass loss rate took from the ppt we made
    # # stars.dmdt = (stellar.particles.mass-stars.mass)/dt
    # inner_binary.Mwind = 0 | units.MSun
    # inner_binary[
    #     0
    # ].terminal_wind_velocity = v_inf_wind_WC  # wind_velocity took from the ppt
    # inner_binary[
    #     1
    # ].terminal_wind_velocity = v_inf_wind_WN  # wind_velocity took from the ppt
    # stellar.stop()

    return inner_binary


def Initialize_apep():
    """
    Create a Particleset of the entire apep system.
    The carbon Wolf-Rayet star is the primary,
    the nitrogen Wolf-Rayet star is the secondary.
    The O8-Iaf supergiant is the tertiary.
    """
    inner_binary = Initialize_inner_binary()

    companion_binary = new_binary_from_orbital_elements(
        mass1=M_WN + M_WC,
        mass2=M_SG,
        semimajor_axis=d_WR_binary_to_SG,
        eccentricity=ecc_SG,
        inclination=inc_SG,
        longitude_of_the_ascending_node=45 | units.deg,
    )
    companion_binary.position -= companion_binary[0].position
    companion_binary.velocity -= companion_binary[0].velocity

    inner_binary.add_particle(companion_binary[1])
    # inner_binary.move_to_center()
    ternary_system = inner_binary
    setattr(ternary_system, "name", ["WC8", "WN46b", "O8"])
    return ternary_system


# print(Initialize_inner_binary())
# print(Initialize_apep())
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_stars_3d(stars):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x_coords = stars.x.value_in(units.au)
    y_coords = stars.y.value_in(units.au)
    z_coords = stars.z.value_in(units.au)
    ax.scatter(x_coords, y_coords, z_coords, c='blue', marker='o')
    ax.set_xlabel('X (AU)')
    ax.set_ylabel('Y (AU)')
    ax.set_zlabel('Z (AU)')
    plt.show()

stars = Initialize_apep()
plot_stars_3d(stars)
"""
