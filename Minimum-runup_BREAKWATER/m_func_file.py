##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np

def m_func(x,y):   
    m = 0.5*np.multiply( (np.sign(x) + np.sign(y)) , np.minimum(np.abs(x), np.abs(y)) ) ;
    return  m