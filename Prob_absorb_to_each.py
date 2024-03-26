# import libraries
import numpy as np
import itertools
from Rate_transitions import get_transition, get_transition_Hughes
from Finding_dictionary_keys import find_keys

def prob_reach_absorb(state_dict, trans_dict, absorb_state, params_dict):
    '''Returns the probabilities of reaching the inputted absorbing state from each possible transient state.'''
    # initialise a matrix Q for the transition rate q_ij betweeen states i and j
    n_states = len(state_dict)               # total number of states
    Q = np.zeros((n_states, n_states))

    # identify and calculate the transition rate for every state pair and store in Q matrix
    for (s1, s2) in itertools.product(range(n_states), repeat=2):
        Q[s1,s2] = get_transition(state_dict[s1], state_dict[s2], params_dict)

    # the diagonal elements of Q are the negative row sums
    rowsumQ = np.sum(Q,1)
    for s in range(n_states):
        Q[s,s] = -rowsumQ[s]

    # initialise a matrix Qcc for the transition rate q_ij betweeen transient states i and j
    n_transient = len(trans_dict) 
    Qcc = np.zeros((n_transient, n_transient))
    qc = np.zeros((n_transient,1))
    
    # identify and calculate the transition rate for every state pair and store in Qcc matrix
    for (s1, s2) in itertools.product(range(n_transient), repeat=2):
        Qcc[s1,s2] = get_transition(trans_dict[s1], trans_dict[s2], params_dict)

    # the diagonal elements of Qcc are equivalent to the diagonal elements of the respective rows of Q
    for s in range(n_transient):     # looping over all the transients states
        diag = find_keys(state_dict,trans_dict[s])[0]    # find key for current transient state
        Qcc[s,s] = Q[diag,diag]
        # transition rate from current transient state to absorbing 
        qc[s,0] = get_transition(trans_dict[s],absorb_state, params_dict)   
    
    soln = np.linalg.solve(-Qcc, qc)  # solving the full linear system of equations
    return soln, Qcc     # returns the probabilities and the sub-q matrix of transient states


def prob_reach_absorb_Hughes(state_dict, trans_dict, absorb_state, params_dict):
    '''Returns the probabilities of reaching the inputted absorbing state from each possible transient state.'''
    # initialise a matrix Q for the transition rate q_ij betweeen states i and j
    n_states = len(state_dict)              # total number of states
    Q = np.zeros((n_states, n_states))

    # identify and calculate the transition rate for every state pair and store in Q matrix
    for (s1, s2) in itertools.product(range(n_states), repeat=2):
        Q[s1,s2] = get_transition_Hughes(state_dict[s1], state_dict[s2], params_dict)

    # the diagonal elements of Q are the negative row sums
    rowsumQ = np.sum(Q,1)
    for s in range(n_states):
        Q[s,s] = -rowsumQ[s]

    # initialise a matrix Qcc for the transition rate q_ij betweeen transient states i and j
    n_transient = len(trans_dict)                 
    Qcc = np.zeros((n_transient, n_transient))
    qc = np.zeros((n_transient,1))
    
    # identify and calculate the transition rate for every state pair and store in Qcc matrix
    for (s1, s2) in itertools.product(range(n_transient), repeat=2):
        Qcc[s1,s2] = get_transition_Hughes(trans_dict[s1], trans_dict[s2], params_dict)

    # the diagonal elements of Qcc are equivalent to the diagonal elements of the respective rows of Q
    for s in range(n_transient):        # looping over all the transients states
        diag = find_keys(state_dict,trans_dict[s])[0]    # find key for current transient state
        Qcc[s,s] = Q[diag,diag]
        # transition rate from current transient state to absorbing 
        qc[s,0] = get_transition_Hughes(trans_dict[s], absorb_state, params_dict)
    
    soln = np.linalg.solve(-Qcc, qc)     # solving the full linear system of equations
    return soln,Qcc                      # returns the probabilities and the sub-q matrix of transient states