
try:
    import numpy as np  
    import random 
    import math
    from constants import *
    from errors import  error_exit
except ImportError as ime:
    if ime.name is None: 
        print("Import in some of the modules ({})",
              " in snafu dir failed. Exiting...".format(ime.name))
        exit(1)
 
def calc_lz_hopp(method, state, pot_eners,
                 pot_eners_array, Ekin, dt):
   
    hop = False             
    instate = state

    # add last calculated pot energies to evaluate hop 
    # vstack along row , axis = 0
    pot_eners_array = np.vstack((pot_eners_array, pot_eners))  
    print("Pot_eners_array: {}".format(pot_eners_array))

    prob = [calc_prob(instate, outstate, pot_eners_array, dt) 
            for outstate in range(0, len(pot_eners)) if outstate != instate ]

    max_prob_row = (np.argmax(prob, axis = 0)[3])  # row with max probaBILITY 
    max_prob = prob[max_prob_row][3]
    theta = random.random()
    print("Max probability: {}".format(prob[max_prob_row][3]))
    
    if prob[max_prob_row][3] > theta:
        outstate = prob[max_prob_row][1]

        dEpot = pot_eners_array[2][outstate] - pot_eners_array[2][instate]

        if dEpot < Ekin: 
            hop = True
            print("dEpot {:.4f}  < Ekin: {:.4f}".format(float(dEpot), Ekin),
                  "Hop = {},".format(hop),
                  "Outstate: {}, Instate: {}".format(instate, outstate),       
                  "Probability: {}".format(prob[max_prob_row][3]),
                  "\nRandon number: {}".format(theta))
            pot_eners_array = np.delete(pot_eners_array, obj = -1, axis = 0) 
            # hop to another PES, remove appended energies from last propagation, comment for test
            
    if not hop:
        pot_eners_array = np.delete(pot_eners_array, 0, axis = 0)
        outstate = instate
    # uncoment for test     pot_eners_array = np.delete(pot_eners_array, 0, axis = 0)    
    return(hop, outstate, pot_eners_array, max_prob)
   
def calc_prob(instate, outstate, pot_eners_array, dt):
    """
    PHYSICAL REVIEW A 84, 014701 (2011)
    Nonadiabatic nuclear dynamics of atomic 
    collisions based on branching classical trajectories
    Andrey K. Belyaev1,2 and Oleg V. Lebedev1
    """
    
    #print("Instate {}, Outstate {}".format(instate, outstate))
    
    # energy gap in time T-DT [0], T [1], T+DT [2]
    Z = [ (abs(pot_eners_array[step][instate] - pot_eners_array[step][outstate])) for step in range(0,3)]  
    print("Z[t-dt]: {:.4f}  Z[t]: {:.4f}  Z[t+dt]: {:.4f}".format(Z[0],Z[1],Z[2]))
    
    # three point minima of adiaabatic splitting Zjk
    if Z[0] > Z[1] and Z[1] < Z[2]: 
        sec_der = ((Z[2] - 2 * Z[1] + Z[0]) / (dt ** 2))
        prob = math.exp(-math.pi/(2 * HBAR_AU) 
                        * (math.sqrt(Z[1] ** 3 / sec_der)))
        if prob > 1.00:
            error_exit(6)
        print("Z({}->{}) minimum  with".format(instate, outstate), 
              "dE/a.u. = {}. ".format(Z[1]),
              "Second derivative at minima: {},".format(sec_der),
              "Z**3: {}".format( Z[1]**3),
              "Probability: {}".format(prob))
    else:
        prob = 0.0000
    return([instate, outstate, Z[1], prob])
