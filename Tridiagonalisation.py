# import libararies
import numpy as np

def tridiagonal(Qk, state_dict):
    '''Takes a square matrix Qk and swaps rows and columns as neccessary to produce tridiagonal matrix.'''
    columns = Qk.shape[1]  # no. of columns
    rows = Qk.shape[0]     # no. of rows
    key_list = [key for key, value in state_dict.items()]    # list of states in class
    
    # column opperations
    for col in range(int(columns/2)):     # loop over first half of the columns 
        if all(Qk[:,col][2+col:] != 0):   # check if column in the correct position has the right no. of zeros
            for cols in range(columns):   # if not loop over the columns
                if np.all(Qk[:, cols][2+col:] == 0):   # search for the correct column to make a swap with
                    col_indx = cols                    # record which column to swap with
            Qk[:,[col,col_indx]] = Qk[:,[col_indx,col]]   # swap columns
            key_list[col], key_list[col_indx] = key_list[col_indx], key_list[col] # reorder states of class list 
            

    # row opperations    
    for row in range(int(rows/2)):    # loop over first half of rows    
        if all(Qk[row,:][2+row:] != 0):   # check if row in correct position
            for rs in range(rows):        # if not loop over rows
                if np.all(Qk[rs, :][2+row:] == 0):  # search for correct row to make a swap with
                    row_indx = rs                   # record which row to swap with
            Qk[[row,row_indx]] = Qk[[row_indx,row]]  # swap rows
            
    return Qk,key_list  # return tridiagonal matrix and reordered key_list