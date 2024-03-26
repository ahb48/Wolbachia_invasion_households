# import libraries
import itertools
from Rate_transitions import get_transition, get_transition_Hughes
import numpy as np

def getQ(state_dict,params_dict):
    '''Constructs the full Q matrix. This is not the Q matrix ordered into lower trianguler form.'''
    n_states = len(state_dict)           # number of states is equivalent to the length of the state dictionary.
    Q = np.zeros((n_states, n_states))   # initialising Q matrix
    # identify and calculate the transition rate for every state pair and store in Q matrix
    for (s1, s2) in itertools.product(range(n_states), repeat=2):  # looping over each state pairing
        Q[s1,s2] = get_transition(state_dict[s1], state_dict[s2], params_dict)

    # the diagonal elements of Q are the negative row sums
    rowsumQ = np.sum(Q,1)  # find the row sums
    for s in range(n_states):
        Q[s,s] = -rowsumQ[s]

    return(Q)   # return the full Q matrix