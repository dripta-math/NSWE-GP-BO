##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
from Flux_file import Flux
import numpy as np

def Roe_Flux(WL, WR, N, eps, g):
    
    hL = WL[:,0,:] ;
    hR = WR[:,0,:] ;
    
    hL_sq = np.sqrt(hL) ;
    hR_sq = np.sqrt(hR) ;
    
    vL = WL[:,1,:]/(WL[:,0,:] + eps) ;
    vR = WR[:,1,:]/(WR[:,0,:] + eps) ;
    
    mu_1 = 0.5*(hL + hR) ;
    mu_2 = np.divide( np.multiply(hL_sq, vL) + np.multiply(hR_sq, vR) , (hL_sq + hR_sq +eps) );
    
    c = np.sqrt(g * mu_1);
    
    s1 = np.sign(mu_2 + c) ;
    s2 = np.sign(mu_2 - c) ;
    
    # U matrix
    U_11 = np.divide( np.multiply(s2, mu_2 + c) - np.multiply(s1, mu_2 - c)  , 2*(c+eps) ) ;
    U_12 = np.divide(s1 - s2 , 2*(c+eps) );
    U_21 = np.divide( np.multiply((s2 -s1), (np.multiply(mu_2, mu_2) - np.multiply(c,c)) ) , 2*(c+eps) ); 
    U_22 = np.divide( np.multiply(s1, mu_2 + c) - np.multiply(s2, mu_2 - c), 2*(c+eps) );
    
    FR = Flux(WR, N, eps, g) ;
    FL = Flux(WL, N, eps, g) ;
    
    Fdiff = (FR - FL)/2;
    Fadd  = (FR + FL)/2;

    F_roe = np.zeros((N-1,2,1));
    F_roe[:,0,:] = Fadd[:,0,:] - (np.multiply(U_11, Fdiff[:,0,:]) + np.multiply(U_12, Fdiff[:,1,:]));  
    F_roe[:,1,:] = Fadd[:,1,:] - (np.multiply(U_21, Fdiff[:,0,:]) + np.multiply(U_22, Fdiff[:,1,:]));      
        
    return F_roe