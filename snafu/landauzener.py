
try:
    import numpy as np  
    import random 
    import math
    from constants import *
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
    pot_eners_array = np.vstack((pot_eners_array, pot_eners))  # vstack along row , axis = 0
    print("Pot_eners_array: {}".format(pot_eners_array))
    
    prob = [calc_prob(instate, outstate, pot_eners_array, dt) for outstate in range(0,len(pot_eners)) if outstate != instate ]
    print("Probabilities: {}".format(prob))

    max_prob_row = (np.argmax(prob, axis = 0)[3]) # find row with max probaBILITY 
    theta = random.random()

    print("Maximal probability: {} | randon number: {}".format(prob[max_prob_row][3], theta))

    if prob[max_prob_row][3] > theta:
        outstate = prob[max_prob_row][1]
        instate = prob[max_prob_row][0]
        dEpot = pot_eners_array[2][outstate] - pot_eners_array[2][instate]
        if dEpot < Ekin: 
            hop = True
            print("dEpot {}  < Ekin: {} , Hop = {}, Instate: {}, Outstate: {}".format(dEpot, Ekin, hop, instate, outstate))
            pot_eners_array = np.delete(pot_eners_array, obj = -1, axis = 0) # hop to another PES, remove appended energies from last propagation
            
    if not hop:
        pot_eners_array = np.delete(pot_eners_array, 0, axis = 0)
        outstate = instate
        
    return(hop, outstate, pot_eners_array)
   
def calc_prob(instate, outstate, pot_eners_array, dt):

    minima = False
    # three point minima
    Z = [ (abs(pot_eners_array[step][instate] - pot_eners_array[step][outstate])) for step in range(0,3)]  # 1: T-DT, 2: T, 3: T+DT
    
    print("Z[t-dt]: {}\nZ[t]: {}\nZ[t+dt]: {}".format(Z[0],Z[1],Z[2]))
    if Z[0] > Z[1] and Z[1] < Z[2]:
        minima = True      
        print("Calculating hop prob. between {} - {} state with dE/a.u. = {} ".format(instate, outstate, Z[1]))
        sec_der = ((Z[2] - 2 * Z[1] + Z[0]) / (dt ** 2))
        print("sec_derivative: {}".format(sec_der))
        prob = math.exp(-math.pi/(2 * HBAR_AU) * (math.sqrt(Z[1]**3 / sec_der) ))
        if prob > 1:
            print("prob larger than 1, somthing wrong")
            exit(1)
        print(prob)
    else:
        prob = 0.0
    return([instate, outstate, Z[1], prob])
