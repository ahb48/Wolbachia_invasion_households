# importing libraries
import numpy as np
from scipy.linalg import expm
from Finding_dictionary_keys import find_keys

def Pget(t_start,t_range,Q,steps,initial_state,state_dict):
    '''Returns the probability distribution for the range of time points entered given by the solution of the ME.'''
    n_states = len(Q)                   # number of states is equal to one of the dimensions of the transition matrix Q
    P_vec = np.zeros((steps,n_states))  # initialise P vector to store probabilities
    t = np.linspace(t_start,t_start+t_range,steps)    # range of time points calculating probabilities over
    P0 = np.zeros(n_states)             # initialising array for initial probability distribution
    indx = find_keys(state_dict,initial_state)   # finding the keys/ index of the initial conditon
    P0[indx] = 1   # setting the probability of being in the initial state defined by the above keys to 1
    
    for i in range(steps):  # looping over each time point
        P_vec[i,:] = P0@expm(Q*t[i])   # calculating and storing the probability distribution at time t
    
    return P_vec,t   # return probability vector and time range