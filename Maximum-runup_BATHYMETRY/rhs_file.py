##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np
from m_func_file import m_func
from IBV_LHS_file_1 import IBV_LHS
from Roe_flux_file import Roe_Flux

def rhs(W, t, N, dx, eps, h, g, Amp, T_per, xc):
    
    h_hu = np.zeros((N,2,1));
    h_hu[:,0,:] = W[:,0,:] ;
    h_hu[:,1,:] = np.divide(W[:,1,:], (W[:,0,:] + eps)) ; 
   
    dW = h_hu[1:, :, :] - h_hu[0:-1, :, :] ;
    DW = h_hu[2:, :, :] - 2*h_hu[1:-1, :, :] + h_hu[0:-2, :, :] ;
    
    D12W  = m_func(DW[0:-1, : , : ], DW[1:, : , : ]) ;

    Splus  = dW[2:-1 , : , :]  - 0.5 * D12W[1:, : , :] ; 
    Sminus = dW[1:-2, : , :] + 0.5 * D12W[0:-1, :, :] ;             

    S = m_func(Splus, Sminus); 

    WL = np.zeros((N-1,2,1));
    WR = np.zeros((N-1,2,1));
                                     
    WL[2:(N-2),0,:] = h_hu[2:(N-2), 0, :] + 0.5 * S[:,0,:]; 
    WR[1:(N-3),0,:] = h_hu[2:(N-2), 0, :] - 0.5 * S[:,0,:];                  
                  
    div1 = np.divide(WR[1:N-3,0,:], (W[2:N-2,0,:] + eps));
    WL[2:(N-2),1,:] = h_hu[2:(N-2), 1, :] + 0.5 * np.multiply(S[:,1,:], div1);  

    div2 = np.divide(WL[2:N-2,0,:], (W[2:N-2,0,:] + eps));
    WR[1:(N-3),1,:] = h_hu[2:(N-2), 1, :] - 0.5 * np.multiply(S[:,1,:], div2);
                  
    # Edge conditions for cells near the left-side boundary
    WL[0,:,:] = h_hu[0,:,:] ;
    WR[0,:,:] =  0.375*h_hu[0,:,:] + 0.75*h_hu[1,:,:] - 0.125*h_hu[2,:,:] ;
    WL[1,:,:] = -0.125*h_hu[0,:,:] + 0.75*h_hu[1,:,:] + 0.375*h_hu[2,:,:] ;

    # Edge conditions for cells near the right-side boundary
    WR[N-2,:,:] = h_hu[N-1,:,:] ;
    WL[N-2,:,:] = -0.125*h_hu[N-3,:,:] + 0.75*h_hu[N-2,:,:] + 0.375*h_hu[N-1,:,:] ;
    WR[N-3,:,:] =  0.375*h_hu[N-3,:,:] + 0.75*h_hu[N-2,:,:] - 0.125*h_hu[N-1,:,:] ;


    hL = h[0:-1,:] ;
    hR = h[1:,:];

    hi = np.minimum(hL, hR) ;         # Element wise minimum
    HLs = np.maximum( W[0:-1,0,:] - hL + hi, np.zeros((N-1,1)) ) ;   # Element wise maximum
    HRs = np.maximum( W[1:,0,:]   - hR + hi, np.zeros((N-1,1)) ) ;

    WL[:,0,:] = HLs ;
    WL[:,1,:] = np.multiply(HLs, WL[:,1,:]) ;

    WR[:,0,:] = HRs ;
    WR[:,1,:] = np.multiply(HRs, WR[:,1,:]) ;

    #####################################################
    F_num = Roe_Flux(WL, WR, N, eps, g) ;

    rhs = np.zeros((N,2,1)); 

    rhs[0:-1,:,:] = rhs[0:-1,:,:] - F_num ;
    rhs[1:,:,:]   = rhs[1:,:,:]   + F_num ;

    rhs[0:-1,1,:] = rhs[0:-1,1,:] - 0.5*g*( np.power(W[0:-1,0,:],2) - np.power(HLs,2) ) ;
    rhs[1:,1,:]   = rhs[1:,1,:]   + 0.5*g*( np.power(W[1:,0,:],2) - np.power(HRs,2) ) ;

    # Boundary conditions at left-side
    x0 = min(xc[:,0]);
    h_LBC  = IBV_LHS(Amp, T_per, t, h) ;
    c = np.sqrt(g*WL[0,0,0]);         # speed
    u_LBC = WL[0,1,0]/(WL[0,0,0]+eps) + (1 - h_hu[0,0,0]/h_LBC)*c ;
    Flux_LBC_1 = h_LBC * u_LBC ;
    Flux_LBC_2 = h_LBC * np.power(u_LBC,2) + 0.5 * g * np.power(h_LBC,2) ;
    rhs[0,:,:] = rhs[0,:,:] + np.array([ [Flux_LBC_1], [Flux_LBC_2] ]); 

    # Boundary conditions at right-side
    h_RBC = h_hu[N-1,0,0] ;
    u_RBC = 0 ;
    Flux_RBC_1 = h_RBC * u_RBC ;
    Flux_RBC_2 = h_RBC * np.power(u_RBC,2) + 0.5 * g * np.power(h_RBC,2) ;
    rhs[N-1,:,:] = rhs[N-1,:,:] - np.array([ [Flux_RBC_1], [Flux_RBC_2] ]) ;   
    
    rhs = rhs/dx; 
    
    return rhs 