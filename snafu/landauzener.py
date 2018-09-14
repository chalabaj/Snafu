
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
 
def calc_hopp(method, state, pot_eners,
              pot_eners_array, Ekin, dt):
    
    """
    Calculate hopping probability according to the following paper
    Landau Zener type surface hopping algorithms
    Andrey K. Belyaev, Caroline Lasser, and Giulio Trigila
    The Journal of Chemical Physics 140, 224108 (2014) 
    Also factor for velocity rescaling is determined by energy conservation
    """

    hop = False             
    instate = state

    # add last calculated pot energies to evaluate hop 
    # vstack along row , axis = 0
    pot_eners_array = np.vstack((pot_eners_array, pot_eners))  
    #print("Pot_eners_array HOP:\n",
    #      "{}".format(pot_eners_array))

    probs = [calc_prob(instate, outstate, pot_eners_array, dt) 
             for outstate in range(0, len(pot_eners)) if outstate != instate]

    max_prob_row = (np.argmax(probs, axis = 0)[3])  # row with max probaBILITY 
    max_prob = probs[max_prob_row][3]
    theta = random.random()
    #print("Max probability: {}".format(probs[max_prob_row][3]))
    
    if max_prob > theta:
        outstate = probs[max_prob_row][1]
        dEpot = pot_eners_array[1][outstate] - pot_eners_array[1][instate]
        # energy conservation criteria
        if dEpot < Ekin and dEpot < 0.02: 
            hop = True
            print("Hop {} --> {}".format(instate, outstate),
                  "dEpot {:.4f}  < Ekin: {:.4f}".format(float(dEpot), Ekin),     
                  "\nProbability: {} ".format(probs[max_prob_row][3]),
                  "Randon number: {}".format(theta))
            v_scaling_fac = math.sqrt(1-(dEpot / Ekin))  # lower upper
    else:
        v_scaling_fac = -1
    if not hop:
        outstate = instate
        v_scaling_fac = 1
    # uncoment for test pot_eners_array = np.delete(pot_eners_array, 0, axis = 0)    
    return(hop, outstate, v_scaling_fac, max_prob)

def calc_prob(instate, outstate, pot_eners_array, dt):
    
    """
    PHYSICAL REVIEW A 84, 014701 (2011)
    Nonadiabatic nuclear dynamics of atomic 
    collisions based on branching classical trajectories
    Andrey K. Belyaev1,2 and Oleg V. Lebedev1
    """

    #print("Instate {}, Outstate {}".format(instate, outstate))
    # energy gap in time T-DT [0], T [1], T+DT [2]
    prob = 0.0000
    Z = [(abs(pot_eners_array[step][instate] 
              - pot_eners_array[step][outstate])) for step in range(0,3)]
    #print("Z[-dt]: {:.4f}  Z[t]: {:.4f} Z[dt]: {:.4f}".format(Z[0],Z[1],Z[2]))

    # three point minima of adiaabatic splitting Zjk
    if Z[0] > Z[1] and Z[1] < Z[2]: 
        sec_der = ((Z[2] - 2 * Z[1] + Z[0]) / (dt ** 2))
        prob = math.exp(-math.pi/(2 * HBAR_AU) 
                        * (math.sqrt(Z[1] ** 3 / sec_der)))
        if prob > 1.00:
            error_exit(6)
        print("Z({}->{}) minimum  with".format(instate, outstate),
              "dE/a.u. = {}. ".format(Z[1]),
              #"Second derivative at minima: {},".format(sec_der),
              #"Z**3: {}".format( Z[1]**3),
              "Probability: {}".format(prob))

    return([instate, outstate, Z[1], prob])
