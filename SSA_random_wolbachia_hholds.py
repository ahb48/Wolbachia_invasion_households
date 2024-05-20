# import required libraries
import random
import numpy as np

''' Household class'''           
class household: 
    '''Constructor method used to initialise household attributes'''
    def __init__(self, m, w, b, u, v, K, d, h, k, delta, phi, tau, num):
        self.m = m    # no. of wild-type mosquitoes in the household
        self.w = w    # no. of wolbachia-infected mosquitoes in the household
        
        self.num = num       # no. of households of this type
        
        # propensities - rates per household are fixed for each object in the class, but total rate depends
        # on number of households
        
        # total birth rates, just used for calculation, not stored
        # if no mosquitoes, no births
        # do I need to add a check for atleast two mosquitoes for a birth?
        
        if (self.m + self.w == 0) or (self.m+self.w >= K):  # and if more than K mosquitoes in household, density dependence
            self.birth_m = 0                                # is set to zero
            self.birth_w = 0
        # otherwise density dependent birth rate, 
        else:
            # this birth rate formulation taken from Hughes
            #self.birth_m = b*((self.m*(self.m + (1-v)*phi*self.w) + self.w*((1-u)*self.m + \
            #                        (1-v)*phi*self.w))/(self.m + self.w))*(1-(self.m+self.w)/K)
            #self.birth_w = b*v*phi*self.w*(1-(self.m+self.w)/K)
            self.birth_m = b*((self.m*(self.m + (1-v)*phi*self.w) + self.w*((1-u)*self.m + \
                                    (1-v)*phi*self.w))/(self.m + self.w))*(np.exp(-h*(self.m + self.w)**k))
            self.birth_w = b*v*phi*self.w*(np.exp(-h*(self.m + self.w)**k))

        # wild-type and wolbachia-infected death rates
        self.death_m = d*self.m
        self.death_w = d*delta*self.w
        
        # wild-type/ wolbachia infected leave rates
        self.leave_m = tau*self.m
        self.leave_w = tau*self.w
        
        # construct the event list and weights for event determination in the get_household_event function, 
        # these never change
        # the weights are equal to the propensities since everything scales with the number of households
        self.events = ['birth_m', 'birth_w', 'death_m', 'death_w', 'leave_m', 'leave_w']
        self.event_weights = [self.birth_m, self.birth_w, self.death_m, self.death_w, self.leave_m, self.leave_w]
        
        # sum the propensities to get total rate event for one household with the given configuration, which never changes
        self.prop_ind = sum(self.event_weights)
      
        # multiply the propensity for each household by the number of households with this configuration  
        # the total rate event for households of this configuration, which changes when the number
        # of households with the configuration changes
        
        self.prop = self.prop_ind*self.num

    '''class function that updates the total number of households in the configuration and total propensity'''
    # call this function whenever the number of households with that configuration changes    
    def update_number(self, increment):
        self.num += increment 
        self.prop = self.prop_ind*self.num

    '''class function that determines which event occurs'''
    # the SSA returns that a broad event has occurred for a household of a given configuration, this function 
    # splits that broad event into a more specific event
    def get_household_event(self):
        # the list of events and their weights are a fixed attribute of the object, just need to select one at random
        event = random.choices(self.events, weights = self.event_weights)[0] 
        # return the event to be processed
        return event        
        
'''Free class - define attributes of the free-living mosquitoes'''           
class free:
    '''Constructor method used to initialise attributes'''    
    def __init__(self, m, w, d, delta, rho):    
        self.m = m        # free wild-type mosquitoes
        self.w = w        # free Wolbachia-infected mosquitoes
        self.num = 0      # number of households in this state - always equal to 0 since not a household state
        
        # rate parameters that never change after initialisation
        
        # we have not included deaths in the disperser pool but we could do, probably evidence for it
        #self.mort_rate_m = d         # mortality rate of wildtype male 
        #self.mort_rate_w = d*delta   # mortality rate of infected male
        
        self.leave_rate = rho        # leave rate of any mosquito
        
        # rates that do change, initialise at 0, then set using the update function below
        #self.death_m = 0
        #self.death_w = 0
        self.leave_m = 0
        self.leave_w = 0
        self.prop = 0    # sum of all event rates for the object
        
        self.update_rates()
     
    '''rate update method'''    
    # updates rate parameters that change when population size changes
    def update_rates(self):
       # self.death_m = self.mort_rate_m*self.m 
       # self.death_w = self.mort_rate_w*self.w
        self.leave_m = self.leave_rate*self.m 
        self.leave_w = self.leave_rate*self.w  
        self.prop = sum([self.leave_m, self.leave_w]) 
    
    ''' number update method '''
    # updates the total number of free wild-type or Wolbachia-infected and propensities
    def update_number(self, mtype, increment):
        # update population numbers
        if mtype == 'm':
            self.m += increment
        elif mtype == 'w':
            self.w += increment
        
        self.update_rates()
      
    '''class function the determines which event occurs'''
    # the SSA returns that a broad event has occurred for a household of a given configuration, this function 
    # splits that broad event into a more specific event
     
    def get_free_event(self):
        # choose an event at random with weights corresponding to propensities (which depend on population size)
        event = random.choices(['leave_fm', 'leave_fw'], weights = [self.leave_m, self.leave_w])[0] 
        return event

    
'''Other functions - probably should eventually be moved to a separate file to the classes'''        

'''Function to create dictionary keys that summarise household attributes.'''
def make_key_house(m, w):
    key = str(m) + "-" + str(w)
    return key  


'''function to generate event details and update the state dictionary accordingly'''

def do_event(states_dict, event_key, event, params_dict): 
    # Gillespie produces an event_key - the state in which an event occurs, and an event - the name of the event that occurs
    # this function determines which other states are affected by the event, and then updates the state
    # in which the event occurs, and any other states affected
    
    # an event involves two states - the starting state decreases by 1, the end state increases by 1
    # event_key gives the starting state, need to construct the end state
    
    # if the event occurred in the free state
    if event_key == 'free':
        # get the event    
        if (event == 'leave_fm') or (event == 'leave_fw'):
            # leaving the free state removes one individual from the free state, and adds it to a random household
            # get the key of that household randomly, in proportion to number of households in each state
            key_list = [key for key in states_dict]     # get the keys of all extant states
            # get the number of households in that state, note that the 'free' state is included,
            # but always has num = 0 so will never be chosen
            weights_list = [value.num for value in states_dict.values()] 
            house_key_out = random.choices(key_list, weights = weights_list)[0]
            
            # extract household state before event occurs      
            m = states_dict[house_key_out].m
            w = states_dict[house_key_out].w
            
            if (event == 'leave_fm'):    
                # update the free state immediately
                states_dict['free'].update_number('m', -1)
                # construct the key for the household state for which the number of households will be increased by 1
                house_key_in = make_key_house(m+1, w)
            elif (event == 'leave_fw'):
                states_dict['free'].update_number('w', -1)
                house_key_in = make_key_house(m, w+1)
    # if event did not occur in free population, must have been in a household                    
    else: 
        # get the household state in which the event ocurred
        # extract household states before event occurs      
        m = states_dict[event_key].m
        w = states_dict[event_key].w 
 
        # set the key of the household state in which the event occurs, the number of households in this state will
        # be reduced by 1
        house_key_out = event_key
        
        # update household states according to event that occurs.
        if event == 'birth_m':
            house_key_in = make_key_house(m+1, w)
        elif (event == 'birth_w'):
            house_key_in = make_key_house(m, w+1)
        elif (event == 'death_m'):
            house_key_in = make_key_house(m-1, w)
        elif (event == 'death_w'):
            house_key_in = make_key_house(m, w-1)  
        elif (event == 'leave_m'): 
            house_key_in = make_key_house(m-1, w)
            states_dict['free'].update_number('m', 1) 
        elif (event == 'leave_w'):
            house_key_in = make_key_house(m, w-1)
            states_dict['free'].update_number('w', 1)
      #  elif (event == 'unviable'):
      #      house_key_in = None
      #      house_key_out = None
         
    # at this point, the free state has been updated if necessary
    # house_key_out is the key of the household state for which the number of households has been reduced by 1
    # house_key_in is the key of the household state for which the number of households has been increased by 1
    # if either of these keys is None, it means the event did not involve a reducing/increasing households in
    # a given state because that part of the event occurred in the free population
    
    # update the number of households in the affected states
    update_household_number(states_dict, house_key_out, -1, params_dict)
    update_household_number(states_dict, house_key_in, 1, params_dict)
        
    return
    
'''function to update the number of households in a given state'''
        
def update_household_number(states_dict, house_key, inc, params_dict):
    # update the number of households in state house_key by amount inc
    # if state does not exist, create it; if number in state drops to 0, remove it
    
    if house_key is None:
        # no household to update, event must have affected free population only
        return
    else: 
        # if this household type already exists, its key will be in the dictionary, so call its update_number
        # function, which will also recalculate its propensity
        if house_key in states_dict:            
            states_dict[house_key].update_number(inc)
            # if the number of households in the state is 0, remove it
            if (states_dict[house_key].num == 0):
                states_dict.pop(house_key) 
        # if the household does not already exist, add it to the dictionary with an initial number inc
        else:
            m = get_value(house_key, 'm')
            w = get_value(house_key, 'w')
            states_dict[house_key] =  household(m, w, params_dict['b'], params_dict['u'], params_dict['v'], \
                                                params_dict['K'], params_dict['d'], params_dict['h'], \
                                                params_dict['k'],params_dict['delta'], params_dict['phi'], \
                                                params_dict['tau'], inc)
    return 
        
'''function to get the number of mosquitoes of given type from household key'''    
    
def get_value(house_key, component):       
    # get the number of mosquitoes of type component in a household from its key
    # key is of the form m-w so need to cut string at the -
    component_strings = house_key.split('-')  # gives a list of integers in string format
    if component == 'm':
        value = int(component_strings[0])
    elif component == 'w':
        value = int(component_strings[1])

    return value    
        
'''function to get the total population sizes from the state dictionary'''    
    
def get_population_size(states_dict):                      
    # sum up the total population size from the state dictionary
    m_pop, w_pop = 0, 0
                            
    for key, value in states_dict.items():
        if key == 'free':
            m_pop += value.m
            w_pop += value.w
        else:
            m_pop += value.m*value.num  # need to multiply by the number of households in this state
            w_pop += value.w*value.num
                            
    return m_pop, w_pop
                      
'''function to initialise the gillespie simulation'''    
                            
def initialise_gillespie(params_dict):
    '''Initialise dictionary for initial household configurations.'''
    
    wrange = np.arange(1,params_dict['K']+1)
    #uvec = np.load('wolb_enter_dist.npy')#('uvec.npy')
    wrand = random.choices(wrange, weights = params_dict['uvec']) # random.choices(range(1, params_dict['K']+1), k=params_dict['H'])
    wsize = np.zeros(params_dict['K'])
    for k in range(params_dict['K']):
        wsize[k] = wrand.count(k+1)
    
    # initially all households are in state (m0, w0); make a key corresponding to this household state
    key_house = make_key_house(params_dict['m0'], 1) 
    
    # initialise states dictionary with all households (total number H) in this state 
    states_dict = {key_house: household(params_dict['m0'], 1, params_dict['b'], params_dict['u'], \
                                        params_dict['v'], params_dict['K'], params_dict['d'], params_dict['h'], \
                                        params_dict['k'], params_dict['delta'], \
                                        params_dict['phi'], params_dict['tau'], wsize[0])}
    
    for k in range(params_dict['K']-1):
        key_house = make_key_house(params_dict['m0'], k+2)
        states_dict[key_house] = household(params_dict['m0'], k+2, params_dict['b'], params_dict['u'], \
                                        params_dict['v'], params_dict['K'], params_dict['d'], params_dict['h'], \
                                        params_dict['k'], params_dict['delta'], \
                                        params_dict['phi'], params_dict['tau'], wsize[k+1])
        
    # add the free state to the states dictionary, initially m_free0, w_free0 in the free state
    key_free = 'free' 
    states_dict[key_free] = free(params_dict['m_free0'], params_dict['w_free0'], params_dict['d'], \
                                 params_dict['delta'], params_dict['rho'])  

    # set up the time points
    t = params_dict['t_start']                       # set current time to start time
    t_vec = np.arange(params_dict['t_start'], params_dict['t_end']+1, step = 1)  # time points for Gillespie output, daily

    i_out = 1                                        # set time point counter to first point after start point
    t_out = t_vec[i_out]                             # set time for next output

    m_vec = np.zeros_like(t_vec)                     # to store no. of wild-type mosquitoes at each time point
    w_vec = np.zeros_like(t_vec)                     # to store no. of Wolbachia-infected mosquitoes at each time point
    #hhold_size_max = np.zeros_like(t_vec)            # to store maximum household size at each time_point
    #num_revert = np.zeros_like(t_vec)
    tOut = np.zeros_like(t_vec)                      # to store time exact time points used
                            
    m_vec[0] = params_dict['m0']*params_dict['H'] + params_dict['m_free0']   # initial no. of wild-type mosquitoes
    w_vec[0] = np.sum(wrand)*params_dict['H'] + params_dict['w_free0']   # initial no. of Wolbachia-infected mosquitoes
                            
    return states_dict, t, i_out, t_out, t_vec, m_vec, w_vec, tOut
        
'''Function to carry out the Gillespie SSA and return results. Many of the functions we have created above feed into this.'''
    
def gillespie(params_dict):
    # initialise states dictionary, time counters for simulation, vectors of m, w population at each time point
    # scalars of m, w populations now, parameters 
    states_dict, t, i_out, t_out, time_points, m_vec, w_vec, tOut = initialise_gillespie(params_dict)  
    num_revert = 0
    
    while t >= 0:  # keeps iterating until break when passes final output time
        # extract the keys, could use .keys and .values instead of .items() but not sure if order is preserved
        keys_list = list(states_dict.keys())
        prop_list = [value.prop for (key, value) in states_dict.items()]
        
        # find the time until the next reaction    
        prop_sum = sum(prop_list)            # total reaction rate
        u = np.random.random()               # uniform random number from (0, 1)
        delta_t = 1/prop_sum * np.log(1/u)   # time to next reaction
            
        ### find the next reaction    
        # prop_list contains the propensities of an event occuring in each extant household state, plus the free state
        # determine the state involved in the next reaction randomly by 
        # selecting state key from list randomly with weights given by the event propensities
        ## may be better to create a structure that associates keys and propensities to avoid cryptic buggs arising
        ## from unexpected reordering?
        # random.choices returns a list, here of length 1, we just want the first element
        event_key = random.choices(keys_list, weights = prop_list)[0]  

        # event_key is the household state, or free state, where an event has occured
        if event_key == 'free':                                       # if event occurs in free population
            event = states_dict[event_key].get_free_event()           # determine event that occurs
        else:                                                         # if event occurs in household population
            event = states_dict[event_key].get_household_event()      # determine event that occurs
        
        # event is a string indicating which broad event has occured
        # process that event
        do_event(states_dict, event_key, event, params_dict) 
                                       
        # update time 
        t += delta_t
        m_pop, w_pop = get_population_size(states_dict)
        
        if t >= t_out:               # check if past next time point
            tOut[i_out] = t          # if so, record time event occured (should be close to time point)
            m_vec[i_out] = m_pop     # record no. wild-type mosquitoes at time step
            w_vec[i_out] = w_pop
        
            # stop if we've reached or passed the final timepoint
            if (t >= params_dict['time_points'][-1]): break
         
            # update the time point for the next output 
            i_out += 1                 
            t_out = time_points[i_out]
            
            ### checking what is the maximum household size reached
            # make a list of the key-value pairs in the households dictionary
            households_list_out = list(states_dict.items())
            # extract the keys
            keys_list_out = [x[0] for x in households_list_out if x[0]!='free'] 
            num_revert = np.sum([1 for key in keys_list_out if (states_dict[key].w==0 and states_dict[key].m>0)])
            #hhold_size_max[i_out] = np.max([states_dict[key].m+states_dict[key].w for key in keys_list_out])
            
        if  ((m_pop == 0) and (w_pop == 0)) or num_revert > 0: break  # if population has become extinct end sim 
            
    # returns arrays that recorded no. of wild-type, Wolbachia-infected mosquitoes 
    # and the maximum household size at each time step        
    return(m_vec,w_vec,num_revert)   