##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np

def IBV_LHS_Nwave(Amp1, Amp_ratio, t1, t2, t, h):
    
  
    
    h_left = max(h[:,0]);  # depth at the left boundary in meters
    

    T_per = 4*(t1-t2);
    
    t10 = 3*(t1-t2);
    t20 = 2*(t1-t2);
    
    arg1 = (2*np.pi/T_per) * (t - t10);
    arg2 = (2*np.pi/T_per) * (t - t20);
    Amp2 = Amp1*Amp_ratio;
    h_LBC = h_left  +  Amp1*np.power(1/np.cosh(arg1), 2) - Amp2*np.power(1/np.cosh(arg2), 2); 
    
    
    return  h_LBC