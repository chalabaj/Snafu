import numpy as np
# CONVERSION CONSTANTS
AU_FS = 0.02418884326505e0   # atomic units to femtosecs
AU_EV = 27.21138386
AMU = 1822.888484264545e0    # atomic mass unit me = 1 AMU*atomic weight
ANG_BOHR = 1.889726132873e0  # agstroms to bohrs
HBAR_EV = 6.58211889e-16
HBAR_AU = 1.00000e0
BOHR_ANG = 1/ANG_BOHR
np.set_printoptions(precision=5)  # for print and stability testing

#  DEFAULT VALUES 
#  might be overwriten by input.in values
ener_thresh = 1.000
hop_thresh = 0.5
dt = 4.00
hop = False
step = 0
dE = 0.0  # energy change since initial energies
Etot_init = 0.0  # setting variable , total energy at the beginning
Etot_prev = 0.0
sim_time = 0.0
prob = 0.0
tera_mpi = 0
write_freq = 10
restart_freq = 100