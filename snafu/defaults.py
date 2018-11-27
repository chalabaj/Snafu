#  DEFAULT VALUES 
#  might be overwriten by input.in values
ener_thresh = 1.000   # (in eV)    
hop_thresh = 0.5      # (in eV)    
write_freq = 10
restart = 0
restart_freq = 100
vel_adj = 0
timestep = 4
method = "lz"
hop = False
step = 0
dE = 0.0  # energy change since initial energies
Etot_init = 0.0  # setting variable , total energy at the beginning
Etot_prev = 0.0
prob = 0.0
tera_mpi = 0
sim_time = 0.0
max_terachem_time = 500 # wait for 1max_terachem_time s until terachem send data, otherwise too long step
liner = ("_") * 100
