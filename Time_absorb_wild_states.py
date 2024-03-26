# import libraries
import numpy as np
from Prob_absorb_to_each import prob_reach_absorb, prob_reach_absorb_Hughes
from scipy.optimize import fsolve

def absorb_time_wolb(max_pop,initial_guess,params_dict):
    '''Returns the expected time reach the Wolbachia-only state space i.e. invasion is successfull for each possible transient state. For the 3 mosquito model, no reversion.'''
    # full state space dictionary
    state_dict = {index: np.array((i, j)) for index, (i, j) in enumerate([(i, j) for i in range(max_pop + 1) for j in   range(max_pop + 1) if i + j <= max_pop])}
    # transient state space dictionary, under no reversion this is the mixed state space
    trans_dict = {index: np.array((i, j)) for index, (i, j) in enumerate([(i, j) for i in range(1,max_pop + 1) for j in range(1,max_pop + 1) if i + j <= max_pop])}
    n_transient = len(trans_dict)    # the number of transient states
    
    # initialising an array for the probabilities of reaching a given absorbing state
    prob_reach_wolb = np.zeros(n_transient) 
    for i in range(max_pop):               # looping over all the possible absorbing states
        absorb_state = np.array([0,i+1])   # current absorbing state  
        # extracting probability of reaching absorbing state and sub-q matrix of transient states
        ac, Qcc = prob_reach_absorb(state_dict,trans_dict,absorb_state,params_dict)
        prob_reach_wolb[:] += np.transpose(ac)[0]  # adding up all the probabilities of reaching the Wolbachia-only states
        
    def equations(u_values):
        '''Returns set of equations need to solve for expected time to absorption to the Wolbachia-only state space.'''
        eqns = []  # initialising list will append equations to

        for i in range(n_transient):    # looping over all the transient states
            # initialising an array to contain the product of the probability of reaching 
            # the Wolbachia-only state space from each transient state and the corresponding expected time to get there
            au_prod = np.zeros((n_transient + max_pop, 1))   
            for j in range(n_transient):
                au_prod[j] = prob_reach_wolb[j] * u_values[j]
            # we multiply each entry of au_prod by the transition rate of moving from each Wolbachia-only 
            # state to the current transient state and take the sum 
            # the last term accounts for the -a*_(m,w) term in eq (14)
            eqn_i = (Qcc[i,:] @ au_prod[:n_transient]) + prob_reach_wolb[i] 
            eqns.append(eqn_i)       # append equation for ith transient state to list
        return np.concatenate(eqns)    # return set of equations to be solved for expected time to absorption

    # Use fsolve to solve for u_values i.e. the expected times to reach the Wolbachia-only state space
    u_solutions = fsolve(equations, initial_guess)

    return u_solutions  # return solutions


def absorb_time_wolb_Hughes(max_pop,initial_guess, params_dict):
    '''Returns the expected time reach the Wolbachia-only state space i.e. invasion is successfull for each possible transient state. For the 30 mosquito model.'''
    # full state space dictionary
    state_dict = {index: np.array((i, j)) for index, (i, j) in enumerate([(i, j) for i in range(max_pop + 1) for j in   range(max_pop + 1) if i + j <= max_pop])}
    # transient state space dictionary, under no reversion this is the mixed state space
    trans_dict = {index: np.array((i, j)) for index, (i, j) in enumerate([(i, j) for i in range(1,max_pop + 1) for j in range(1,max_pop + 1) if i + j <= max_pop])}
    n_transient = len(trans_dict)    # the number of transient states
    
    # initialising an array for the probabilities of reaching a given absorbing state
    prob_reach_wolb = np.zeros(n_transient)
    for i in range(max_pop):                 # looping over all the possible absorbing states
        absorb_state = np.array([0, i+1])    # current absorbing state  
        # extracting probability of reaching absorbing state and sub-q matrix of transient states
        ac, Qcc = prob_reach_absorb_Hughes(state_dict, trans_dict, absorb_state, params_dict)
        prob_reach_wolb[:] += np.transpose(ac)[0]   # adding up all the probabilities of reaching the Wolbachia-only states
        
    def equations(u_values):
        '''Returns set of equations need to solve for expected time to absorption to the Wolbachia-only state space.'''
        eqns = []     # initialising list will append equations to

        for i in range(n_transient):   # looping over all the transient states
            # initialising an array to contain the product of the probability of reaching 
            # the Wolbachia-only state space from each transient state and the corresponding expected time to get there
            au_prod = np.zeros((n_transient + max_pop, 1))
            for j in range(n_transient):
                au_prod[j] = prob_reach_wolb[j] * u_values[j]
            # we multiply each entry of au_prod by the transition rate of moving from each Wolbachia-only 
            # state to the current transient state and take the sum 
            # the last term accounts for the -a*_(m,w) term in eq (14)    
            eqn_i = (Qcc[i,:] @ au_prod[:n_transient]) + prob_reach_wolb[i]
            eqns.append(eqn_i)          # append equation for ith transient state to list
        return np.concatenate(eqns)     # return set of equations to be solved for expected time to absorption

    # Use fsolve to solve for u_values i.e. the expected times to reach the Wolbachia-only state space
    u_solutions = fsolve(equations, initial_guess)

    return u_solutions   # return solutions