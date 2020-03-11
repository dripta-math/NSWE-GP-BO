##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np

def time_step(W,CFL,dx, eps, g):
    
    v = np.divide(W[:,1,:], (W[:,0,:] + eps)) ;
    dt = CFL*dx/np.max( np.abs(v) + 2*np.sqrt(g*W[:,0,:]) );
    
    return dt