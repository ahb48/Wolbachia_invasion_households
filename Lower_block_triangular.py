# import libraries
from Finding_sub_Q import getQk, getQk_Hughes
import numpy as np
from Rate_transitions import get_transition, get_transition_Hughes

def LBTQ(Q, state_dict, state_dict_S1, state_dict_S2, state_dict_S3, max_pop, params_dict):
    '''Finding the sub-q matrices and their respective ordered lists of states in the class. This is for the 3 mosquito model.'''
    # when no reversion is possible there are 3 communicating classes
    # (1) the wild-type-only, (2) the Wolbachia-only, (2) the mixed states
    Q1,key_list1 = getQk(state_dict_S1,state_dict,Q,params_dict) # finding the sub-q matrices for each communicating class
    Q2,key_list2 = getQk(state_dict_S2,state_dict,Q,params_dict)
    Q3,key_list3 = getQk(state_dict_S3,state_dict,Q,params_dict)
                                       
    Q_lower_block_triang = np.zeros_like(Q)  # initialing the Q matrix which we will write in lower block triangular form
    # it still has the same shape as Q in it's original form
    Q_lower_block_triang = np.delete(Q_lower_block_triang, 0, axis=0)   # but we delete the rows and columns corresponding
    Q_lower_block_triang = np.delete(Q_lower_block_triang, 0, axis=1)   # to the (0,0) state
    # initialising reordered state dictionary with class S1
    # this will match the order of the Q matrix in lower block triangular form
    state_dict_relabel = {index: np.array((i, 0)) for index, i in enumerate([i for i in range(1,max_pop + 1)])}
    
    S1_len = Q1.shape[0]  # row/column length of S1 class
    S2_len = Q2.shape[0]  # S2 class
    S3_len = Q3.shape[0]  # S3 class

    for i in range(S1_len):  # populating Q matrix with Q1 along the diagonal starting top left
        for j in range(S1_len):
            Q_lower_block_triang[i,j] = Q1[i,j]
        
    for i in range(S2_len):    # populating Q next with Q2 along the diagonal
        for j in range(S2_len):
            Q_lower_block_triang[S1_len+i,S1_len+j] = Q2[i,j] 
        state_dict_relabel[S1_len+i] = state_dict_S2[i]   # adding S2 states to state dictionary
        
    for i in range(S3_len):      # populating Q lastly along the diagonal with Q3
        for j in range(S3_len):
            Q_lower_block_triang[S1_len+S2_len+i,S1_len+S2_len+j] = Q3[i,j]
        
    for s1 in range(S2_len):   # adding in transition terms from S2 states to S1 (this in the lower triangle portion of Q)
        for s2 in range(S1_len):
            Q_lower_block_triang[S1_len+s1,s2] = get_transition(state_dict_S2[s1], state_dict_S1[s2], params_dict)
        
    for s1 in range(S3_len):
        s3_indx = key_list3[s1] # ensuring correct order of states in S3
        state_dict_relabel[S1_len+S2_len+s1] = state_dict_S3[s3_indx]  # adding S3 states to state dictionary
        for s2 in range(S1_len):  # adding transition terms from S3 to S1
            Q_lower_block_triang[S1_len+S2_len+s1,s2] = get_transition(state_dict_S3[s3_indx], state_dict_S1[s2], params_dict)
        for s2 in range(S2_len):  # adding transition terms from S3 to S2
            Q_lower_block_triang[S1_len+S2_len+s1,S1_len+s2] = get_transition(state_dict_S3[s3_indx], state_dict_S2[s2], params_dict)
        
    return(Q_lower_block_triang,state_dict_relabel)  # returning the lower block triangular matrix and the reordered
                                                     # full state dictionary

def LBTQ_Hughes_comp(Q,state_dict,state_dict_S1,state_dict_S2,state_dict_S3,max_pop,params_dict):
    '''Finding the sub-q matrices and their respective ordered lists of states in the class. This is for the 30 mosquito model without reversion.'''
    # when no reversion is possible there are 3 communicating classes
    # (1) the wild-type-only, (2) the Wolbachia-only, (2) the mixed states
    # finding the sub-q matrices for each communicating class
    Q1,key_list1 = getQk_Hughes(state_dict_S1,state_dict,Q,params_dict)
    Q2,key_list2 = getQk_Hughes(state_dict_S2,state_dict,Q,params_dict)
    Q3,key_list3 = getQk_Hughes(state_dict_S3,state_dict,Q,params_dict)
                                       
    Q_lower_block_triang = np.zeros_like(Q)  # initialing the Q matrix which we will write in lower block triangular form
    # it still has the same shape as Q in it's original form
    Q_lower_block_triang = np.delete(Q_lower_block_triang, 0, axis=0)  # but we delete the rows and columns corresponding
    Q_lower_block_triang = np.delete(Q_lower_block_triang, 0, axis=1)  # to the (0,0) state
    # initialising reordered state dictionary with class S1
    # this will match the order of the Q matrix in lower block triangular form
    state_dict_relabel = {index: np.array((i, 0)) for index, i in enumerate([i for i in range(1,max_pop + 1)])}
    
    S1_len = Q1.shape[0]  # row/column length of S1 class
    S2_len = Q2.shape[0]  # S2 class
    S3_len = Q3.shape[0]  # S3 class

    for i in range(S1_len):  # populating Q matrix with Q1 along the diagonal starting top left
        for j in range(S1_len):
            Q_lower_block_triang[i,j] = Q1[i,j]
        
    for i in range(S2_len):    # populating Q next with Q2 along the diagonal
        for j in range(S2_len):
            Q_lower_block_triang[S1_len+i,S1_len+j] = Q2[i,j] 
        state_dict_relabel[S1_len+i] = state_dict_S2[i]   # adding S2 states to state dictionary
        
    for i in range(S3_len):      # populating Q lastly along the diagonal with Q3
        for j in range(S3_len):
            Q_lower_block_triang[S1_len+S2_len+i,S1_len+S2_len+j] = Q3[i,j]
        
    for s1 in range(S2_len):   # adding in transition terms from S2 states to S1 (this in the lower triangle portion of Q)
        for s2 in range(S1_len):
            Q_lower_block_triang[S1_len+s1,s2] = get_transition_Hughes(state_dict_S2[s1], state_dict_S1[s2], params_dict)
        
    for s1 in range(S3_len):
        s3_indx = key_list3[s1] # ensuring correct order of states in S3
        state_dict_relabel[S1_len+S2_len+s1] = state_dict_S3[s3_indx]  # adding S3 states to state dictionary
        for s2 in range(S1_len):  # adding transition terms from S3 to S1
            Q_lower_block_triang[S1_len+S2_len+s1,s2] = get_transition_Hughes(state_dict_S3[s3_indx], state_dict_S1[s2], params_dict)
        for s2 in range(S2_len):  # adding transition terms from S3 to S2
            Q_lower_block_triang[S1_len+S2_len+s1,S1_len+s2] = get_transition_Hughes(state_dict_S3[s3_indx], state_dict_S2[s2], params_dict)
        
    return(Q_lower_block_triang,state_dict_relabel)  # returning the lower block triangular matrix and the reordered
                                                     # full state dictionary



def LBTQ_Hughes(Q,state_dict,state_dict_S1,state_dict_S2,max_pop,params_dict):
    '''Finding the sub-q matrices and their respective ordered lists of states in the class. This is for the 30 mosquito model with reversion.'''
    # when reversion is possible there are 2 communicating classes
    # (1) the wild-type-only, (2) the mixed states and the Wolbachia-only states
    # finding the sub-q matrices for each communicating class
    Q1,key_list1 = getQk_Hughes(state_dict_S1,state_dict,Q,params_dict)
    Q2,key_list2 = getQk_Hughes(state_dict_S2,state_dict,Q,params_dict)

    Q_lower_block_triang = np.zeros_like(Q)  # initialing the Q matrix which we will write in lower block triangular form
    # it still has the same shape as Q in it's original form
    Q_lower_block_triang = np.delete(Q_lower_block_triang, 0, axis=0)   # but we delete the rows and columns corresponding
    Q_lower_block_triang = np.delete(Q_lower_block_triang, 0, axis=1)   # to the (0,0) state
    # initialising reordered state dictionary with class S1
    # this will match the order of the Q matrix in lower block triangular form
    state_dict_relabel = {index: np.array((i, 0)) for index, i in enumerate([i for i in range(1,max_pop + 1)])}
    
    S1_len = Q1.shape[0]  # row/column length of S1 class
    S2_len = Q2.shape[0]  # S2 class

    for i in range(S1_len):  # populating Q matrix with Q1 along the diagonal starting top left
        for j in range(S1_len):
            Q_lower_block_triang[i,j] = Q1[i,j]
        
    for i in range(S2_len):    # populating Q next with Q2 along the diagonal
        for j in range(S2_len):
            Q_lower_block_triang[S1_len+i,S1_len+j] = Q2[i,j] 
        state_dict_relabel[S1_len+i] = state_dict_S2[i]   # adding S2 states to state dictionary
        
    for s1 in range(S2_len):   # adding in transition terms from S2 states to S1 (this in the lower triangle portion of Q)
        s2_indx = key_list2[s1]
        state_dict_relabel[S1_len+s1] = state_dict_S2[s2_indx]  # adding S2 states to state dictionary
        for s2 in range(S1_len):
            Q_lower_block_triang[S1_len+s1,s2] = get_transition_Hughes(state_dict_S2[s2_indx], state_dict_S1[s2], params_dict)
        
    return(Q_lower_block_triang,state_dict_relabel)  # returning the lower block triangular matrix and the reordered
                                                     # full state dictionary