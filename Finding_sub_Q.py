# import libraries
import numpy as np
import itertools
from Finding_dictionary_keys import find_keys
from Rate_transitions import get_transition, get_transition_rev, get_transition_Hughes
from Tridiagonalisation import tridiagonal

def getQk(state_dict_k, state_dict, Q, params_dict):
    '''Identify and calculate the transition rates for every state pair and store in Q_k matrix. This is for the 3 mosquito model rates'''
    Q_k = np.zeros((len(state_dict_k), len(state_dict_k))) # initialise sub q matrix
    for (s1, s2) in itertools.product(range(len(state_dict_k)), repeat=2): # loop over states
        # the diagonal elements of Qk are the negative row sums of full Q
        if s1 == s2:
            diag = find_keys(state_dict,state_dict_k[s1])[0]
            Q_k[s1,s1] = Q[diag,diag]
        else:   # find the other elements using get_transition function 
            Q_k[s1,s2] = get_transition(state_dict_k[s1], state_dict_k[s2], params_dict)
    # covert the sub q matrix to tridiagonal    
    Q_k,key_list = tridiagonal(Q_k,state_dict_k)
    
    return Q_k,key_list  # return Q_k and ordered list of states


def getQk_Hughes(state_dict_k,state_dict,Q,params_dict):
    '''Identify and calculate the transition rate for every state pair and store in Q_k matrix. This is for the 30 mosquito model with the rates comparable to Hughes mean-field model.'''
    Q_k = np.zeros((len(state_dict_k), len(state_dict_k))) # initialise sub q matrix
    for (s1, s2) in itertools.product(range(len(state_dict_k)), repeat=2): # loop over states
        # the diagonal elements of Qk are the negative row sums of full Q
        if s1 == s2:
            diag = find_keys(state_dict,state_dict_k[s1])[0]
            Q_k[s1,s1] = Q[diag,diag]
        else:   # find the other elements using get_transition function 
            Q_k[s1,s2] = get_transition_Hughes(state_dict_k[s1], state_dict_k[s2], params_dict)
    # covert the sub q matrix to tridiagonal    
    Q_k,key_list = tridiagonal(Q_k,state_dict_k)
    
    return Q_k,key_list  # return Q_k and ordered list of states