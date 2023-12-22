"""
initialize_apep
Created on 26-11-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

- Define the initial conditions of the Apep system, according to the parameters 
from Callingham et al. 2019, 2020; Han et al. 2020; del Palacio et al. 2022.
- Create the initial binary system of the WR stars, and add the O8-Iaf supergiant.
"""
# Importing modules
from amuse.units import units, constants
from amuse.lab import Particles
from amuse.community.seba.interface import SeBa
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

# Wolf Rayet (WR) binary
a_binary = 67 | units.au
P_binary = 125 | units.yr
ecc_binary = 0.7
inc_binary = 0 | units.deg # 25 degrees according to Han et al. 2020, but we keep WR binary in the x-y plane for rotation
omega_binary = 0 | units.deg

M_enclosed_binary = 30 | units.MSun
# A guess for the individual WN46b-type WR (WN) and WC8-type WR (WC) masses:
M_WN = 11 | units.MSun
M_WC = 19 | units.MSun

R_WN = 6 | units.RSun
R_WC = 6.3 | units.RSun

T_eff_WN = 65000 | units.K
T_eff_WC = 60000 | units.K

M_loss_WN = 4 * 10 ** (-5) | units.MSun / units.yr # mass loss rate of the WN46b-tyoe WR
M_loss_WC = 2.9 * 10 ** (-5) | units.MSun / units.yr # mass loss rate of the WC8-type WR

# O8-Iaf Supergiant (SG)
d_WR_binary_to_SG = 1700 | units.au
P_binary_SG = 10**4 | units.yr
M_SG = 40 | units.MSun
R_SG = 35.4 | units.RSun
T_eff_SG = 30000 | units.K
M_loss_SG = 10 ** (-5) | units.MSun / units.yr
ecc_SG = 0  # guess, circular
inc_SG = 25 | units.deg  # guess

# Stellar winds
v_inf_wind_WN = 3500 | units.kms
v_inf_wind_WC = 2100 | units.kms
T_wind_WN = 0.3 * T_eff_WN
T_wind_WC = 0.3 * T_eff_WC

def Initialize_inner_binary():
    """
    Create a Particle set of the inner binary system (the "Central Engine").
    - Carbon burning WR star is the primary
    - Nitrogen burning WR star is the secondary
    """
    inner_binary = new_binary_from_orbital_elements(
        mass1=M_WC,
        mass2=M_WN,
        semimajor_axis=a_binary,
        eccentricity=ecc_binary,
        inclination=inc_binary,
        true_anomaly=90 | units.deg,
    )
    setattr(inner_binary, "name", ["WC8", "WN46b"])

    # Assigining properties to the stars:
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

    return inner_binary


def Initialize_apep():
    """
    Create a Particleset of the entire apep system.
    - Carbon burning WR star is the primary,
    - Nitrogen burning WR star is the secondary,
    - O8-Iaf supergiant is the tertiary.
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
    ternary_system = inner_binary
    setattr(ternary_system, "name", ["WC8", "WN46b", "O8"])

    ternary_system[ternary_system.name == "O8"].radius = R_SG
    ternary_system[ternary_system.name == "O8"].temperature = T_eff_SG
    ternary_system[ternary_system.name == "O8"].dmdt = M_loss_SG

    return ternary_system


# Plot initial system in 3D (for testing purposes)
def plot_initial_apep():
    bodies = Initialize_apep()
    x = bodies.x.value_in(units.AU)
    y = bodies.y.value_in(units.AU)
    z = bodies.z.value_in(units.AU)

    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d

    with plt.rc_context({'axes.edgecolor': 'white',
                            'xtick.color': 'white',
                            'ytick.color': 'white',
                            'figure.facecolor': 'black',
                            'axes.facecolor': 'black',
                            'axes.labelcolor': 'white',
                            'axes.titlecolor': 'white'}):
        fig = plt.figure(dpi=450)
        ax = plt.axes(projection='3d')
        ax.dist = 13
        ax.scatter3D(x, y, z, s=100, c=['blue', 'cyan', 'red'], marker='*')
        ax.set_xlabel('x [AU]')
        ax.set_ylabel('y [AU]')
        ax.set_zlabel('z [AU]')
        for i in plt.legend(handles=[ax.scatter3D([], [], [],marker="*", color='blue', label='WC8'),
                                ax.scatter3D([], [], [], marker="*", color='cyan', label='WN46b'),
                                ax.scatter3D([], [], [], marker="*", color='red', label='O8')],
                        loc='upper center', bbox_to_anchor=(0.5, 1),
                        ncol=2, fancybox=True, shadow=True).get_texts():
            i.set_color("white")
        plt.show()
# plot_initial_apep()