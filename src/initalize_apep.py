"""
initialize_apep
Created on 26-11-23

@author(s): Sam Beckers, Divyansh Srivastava, Vincent van Rie 

{Outline of code}
"""
#Importing modules
from amuse.units import units, constants
from amuse.lab import Particles

#Wolf Rayet (WR) binary
a_binary = 67|units.au
P_binary = 125|units.yr
ecc_binary = 0.7
inc_binary = 25|units.deg
omega_binary = 0|units.deg

M_enclosed_binary = 30 |units.MSun
# A guess for the individual WN46b-type WR (WN) and WC8-type WR (WC) masses:
M_WN = 11|units.MSun
M_WC = 19|units.MSun

R_WN = 6|units.RSun
R_WC = 6.3|units.RSun

T_eff_WN = 65000|units.K
T_eff_WC = 60000|units.K

M_loss_WN = 4*10**(-5)|units.MSun/units.yr #mass loss rate of the WN46b-tyoe WR
M_loss_WC = 2.9*10**(-5)|units.MSun/units.yr #mass loss rate of the WC8-type WR

#O8-Iaf Supergiant (SG)
d_WR_binary_to_SG = 1700|units.au
P_binary_SG = 10**4|units.yr
M_SG = 40|units.MSun
R_SG = 35.4|units.RSun
T_eff_SG = 30000|units.K
M_loss_SG = 10**(-5)|units.MSun/units.yr

#Stellar winds
v_inf_wind_WN = 3500|units.kms
v_inf_wind_WC = 2100|units.kms
T_wind_WN = 0.3*T_eff_WN
T_wind_WC = 0.3*T_eff_WC
