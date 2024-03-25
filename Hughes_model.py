# import libraries
import numpy as np

def Hughes_ODEs(t, n, u, v, phi, delta, b, d, dw, Q, h, k):
    '''Returns ODE outputs for adapted mean-field model under parameter values given.'''
    nm,nw = n  # unpack number of wild-types and number of Wolbachia-infected mosquitoes
    if nm < 0:  # checking populations not going negative
        nm = 0
    if nw < 0:
        nw = 0
    if ((nm == 0) and (nw == 0)):  # if no mosquitoes 
        zm = 0
    else:
        zm = (nm*(nm + (1-v)*phi*nw) + nw*((1-u)*nm + (1-v)*phi*nw))/(nm + nw)  # altered zm
    zw = v*phi*nw
    
    # ODEs
    nm_dot = b*zm*F_hughes(nm+nw,h,k) - d*nm    # describes wild-type population
    nw_dot = b*zw*F_hughes(nm+nw,h,k) - dw*nw   # Wolbachia
    
    return(nm_dot, nw_dot)  # return ODE values

def F(x,Q):
    '''Returns alternative larval density function.'''
    if x >= Q:
        return 0
    else:
        return (1-x/Q)
    
def F_hughes(x,h,k):
    '''Returns Dye's larval density function.'''
    return(np.exp(-h*(x)**k))