##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
def runup_inund_calc(W, xc, h):
    
    
    for i in range(xc.shape[0]):
        
        j = xc.shape[0] - (i+1);
        
        if W[j,0] > 0.01:
            
            runup = - h[j,0];
            inund = xc[j,0];
            
            break
        
        
    return [runup, inund]