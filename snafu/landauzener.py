
try:
    import numpy as np  
    import random  
except ImportError as ime:
    if ime.name is None:  # module could have been removed or module file renamed
        print("Import in some of the modules ({}) in snafu dir failed. Exiting...".format(ime.name))
        exit(1)

      
def calc_lz_hopp(method, state, pot_eners,
                 pot_eners_array, Ekin, dt):
    
                
    hop = False             
    instate = state
    # add last calculated pot energies to evaluate hop 
    # based on three point minima
    pot_eners_array = np.append(pot_eners_array, pot_eners, axis = 0)
    print(pot_eners_array)
    
    prob = [calc_prob(instate, outstate, pot_eners_array) for outstate in range(0,len(pot_eners) if outstate <> instate ]
    print(minima)
    
    
    max_prob = (np.argmax(prob, axis = 0)[3]) # find row with max probaBILITY 
    if prob[max_prob] > random.random()
        if prob[2] < Ekin   # TODO depend if instateou <> outstate???
            hop = True
     
    if conservation
    
    if more than 1 possibilities
    highest prob wins 
      compare with random number 
      
    
    
    # hop to another PES, remove appended energies
    if hop:
        pot_eners_array = np.delete(pot_eners_array, 0, axis = 0)  
    else:
        pot_eners_array = np.delete(pot_eners_array, -1, axis = 0)
    
    return(hop, pot_eners_array)
   
def calc_prob(instate, outstate, pot_eners_array, dt):

    hbar_ev = 6.58211889e-16
    minima = False

    Z = [ (abs(pot_eners_array[step][instate] - pot_eners_array[step][oustate]) * au_eV) for step in range(1,4)]  # 1: T-DT, 2: T, 3: T+DT
    
    print(Z)
    if Z[0] > Z[1] and Z[1] < Z[2]:
        minima = True      
        print("Calculating hop prob. between {} - {} state with dE = {} ".format(minima[0], minima[1], minima[2]))
        sec_der = ((Z[2] - 2 * Z[1] + Z[0]) / (dt ** 2))
        prob = math.exp(-math.pi/(2 * hbar_ev)*(math.sqrt( Z[1]**3 / sec_der ) ))
        print(prob)
    return([instate, outstat, Z[1], prob])
