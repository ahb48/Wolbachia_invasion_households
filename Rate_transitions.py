
def get_transition(state1, state2, p):
    '''Define a function to identify the transition (if there is one) that connects two states.
    The transition is from state 1 to state 2. These are the transition rates corresponding to the 3 mosquito model.'''

    # the change between state 1 and state 2, as an np.array
    state_diff = state2 - state1 
    
    # the only admissible changes for our model are when one population increaes/decreases by 1 and other does not change
    if state_diff[0] == 1 and state_diff[1] == 0:
        # population 1 increases by 1 (reproduction)
        rate = max(0, p['b1']*state1[0]*(1 - sum(state1)/p['K']))  # full wild-type birth rate
    elif state_diff[0] == -1 and state_diff[1] == 0:
        # population 1 decreases by 1 (mortality)
        rate = p['d1']*state1[0]                                     # wild-type death rate
    elif state_diff[0] == 0 and state_diff[1] == 1:
        # population 2 increases by 1 (reproduction)
        rate =  max(0, p['b2']*state1[1]*(1 - sum(state1)/p['K']))  # full Wolbachia birth rate
    elif state_diff[0] == 0 and state_diff[1] == -1:
        # population 2 decreases by 1 (mortality)
        rate = p['d2']*state1[1]                                      # Wolbachia death rate 
    else:
        # all other transition are invalid
        rate = 0
        
    return rate  # returns required rate


from Hughes_model import F_hughes as F # importing Dye's larval density function
# will need to alter for alternative larval density function

def get_transition_Hughes(state1, state2, p):
    '''Define a function to identify the transition (if there is one) that connects two states.
    The transition is from state 1 to state 2. These are the transition rates corresponding to the 30 mosquito model
    With rates comparable to the mean-field model.'''

    # the change between state 1 and state 2, as an np.array
    state_diff = state2 - state1 
    m = state1[0]  # relabelling initial wild-type value, for compactness
    w = state1[1]  # Wolbachia value
    
    # the only admissible changes for our model are when one population increaes/decreases by 1 and other does not change
    if state_diff[0] == 1 and state_diff[1] == 0:
        # population 1 increases by 1 (reproduction)
        if m + w == 0:  # if no mosquitoes
            rate = 0
        else:
            rate = p['b1']*((m*(m + (1-p['v'])*p['phi']*w) + w*((1-p['u'])*m + \
                        (1-p['v'])*p['phi']*w))/(m + w))*F(m+w,p['h'],p['k']) # full wild-type birth rate
    elif state_diff[0] == -1 and state_diff[1] == 0:
        # population 1 decreases by 1 (mortality)
        rate = p['d1']*state1[0]                        # wild-type death rate
    elif state_diff[0] == 0 and state_diff[1] == 1:
        # population 2 increases by 1 (reproduction)
        rate = p['b1']*p['v']*p['phi']*w*F(m+w,p['h'],p['k'])  # full Wolbachia birth rate
    elif state_diff[0] == 0 and state_diff[1] == -1:
        # population 2 decreases by 1 (mortality)
        rate = p['d2']*state1[1]                        # Wolbachia death rate 
    else:
        # all other transition are invalid
        rate = 0
        
    return rate  # return required transition rate