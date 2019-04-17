#  DEFAULT VALUES - do not change unless you know what you do!!!!
#  might be overwriten by input.in values
ener_thresh = 1.000   # (in eV)    
hop_thresh = 0.50      # (in eV)    
write_freq = 10
restart = 0
restart_freq = 100
vel_adj = 0
timestep = 4
method = "lz"
hop = False
step = 0
dE = 0.000         # energy change since initial (t=0) energy
Etot_init = 0.000  # setting variable , total energy at the beginning
Etot_prev = 0.000
prob = 0.000
sim_time = 0.000
max_terachem_time = 10000 # wait for 1max_terachem_time s until terachem send data, otherwise too long step
liner = ("_") * 140

# tera parameters for calc_forces, will be changed only if tera_mpi <>0
comm = 0
MO = 0
CiVecs = 0
NAC = 0
blob = 0
SMatrix = 0
civec_size = 0
nbf_size = 0 
blob_size = 0
qmcharges = 0
TDip = 0
Dip = 0
