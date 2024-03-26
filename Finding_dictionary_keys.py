# importing libraries
import numpy as np

def find_keys(my_dict, target_value):
    '''Using a list comprehension to get keys with the specified value from the dictionary inputted.'''
    keys = [key for key, value in my_dict.items() if np.array_equal(value, target_value)]
    return keys # returns keys