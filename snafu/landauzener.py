

try:
    import numpy as np    
except ImportError as ime:
    if ime.name is None:  # module could have been removed or module file renamed
        print("Import in some of the modules ({}) in snafu dir failed. Exiting...".format(ime.name))
        exit(1)

        
def calc_lz_hopp(method, state, pot_eners,
                 pot_eners_array, Ekin):
    
    # if hop remove last pot ener
                
    hop = False             
    print("lal") 
    # add last calculated pot energies to evaluate hop 
    # based on three point minima
    pot_eners_array = np.append(pot_eners_array, pot_eners, axis = 0)
    
    minima = [minima_search(state, st, pot_eners_array) for st in range(0,len(pot_eners) if st <> state ]
    
    calc_prob(minima[inout][3]) for inout in range(0, len(minima)) if minima[inout][3]     
       
    if minima
     check energy conservation
     
    if conservation
    
    if more than 1 possibilities
    highest prob wins 
      compare with random number 
      
    
    
    # hop to another PES, remove appended energies
    if hop:
        pot_eners_array = np.delete(pot_eners_array, 0, axis = 0 )  
    else:
        pot_eners_array = np.delete(pot_eners_array, -1, axis = 0 )
    
    return(hop, pot_eners_array)
   
def minima_search(instate, outstate, pot_eners_array):

    minima = False

    Z = [ abs(pot_eners_array[step][instate] - pot_eners_array[step][oustate]) for step in range(1,4)]  # T-DT, T, T+DT
    
    if Z[0] > Z[1] and Z[1] < Z[2]:
        minima = True      

    return([statein, stateout, minima]) 