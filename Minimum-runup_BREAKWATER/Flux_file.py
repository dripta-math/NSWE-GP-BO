##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np

def Flux(q, N, eps, g):
    
    F = np.zeros((N-1,2,1))
    F[:,0,:] = q[:,1,:];
    F[:,1,:] = np.divide(np.multiply(q[:,1,:], q[:,1,:]), q[:,0,:]+eps) + 0.5 * g * np.multiply(q[:,0,:], q[:,0,:]);
    
    return  F