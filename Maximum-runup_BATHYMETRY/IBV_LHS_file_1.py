##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np

def IBV_LHS(Amp, T_per, t, h):
    
    omega = 2*np.pi/T_per;
    
    h_left = max(h[:,0]);  
    
#    h_LBC = h_left + Amp * np.sin( omega * t) * np.heaviside(4*T_per - t, 0.5);
    h_LBC = h_left + Amp * np.sin( omega * t);
    
    return  h_LBC