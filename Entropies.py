# import libraries 
import numpy as np 

def entropy(Q_lower_block_triang, quasi_stat_dist):
    '''Calculates individual and total state entropies of the system. 
    Inputs are Q_lower_block_triang: the lower block
    triangular matrix form of Q and quasi_stat_dist: the QSD.'''
    H = 0 # intialise total entropy
    n_LBT_states = len(Q_lower_block_triang)  # define number of states 
    H_vec = np.zeros(n_LBT_states)    # initialise array for individual state entropies

    for i in range(n_LBT_states):   # looping over all the states
        Hi = 0       # initialise current individual state entropy
        for j in range(n_LBT_states):   # looping over all the states
            if (i == j):   # if the same state
                p = 0
            else:             # calculating the individual entropy, see eq (27)-(29)
                p = -Q_lower_block_triang[i,j]/Q_lower_block_triang[i,i]
            if (p == 0):
                plogp = 0
            else:
                plogp = p*np.log(p)
            Hi += plogp    # reassigning the individual entropy value
        H_vec[i] = -Hi     # recording the individual entropy in the array
        H -= quasi_stat_dist[i]*Hi   # updating the overall entropy
    return H_vec, H    # returns array of individual state entropies and the overall entropy