##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import GPy
import numpy as np
import scipy.io 
import os
import matplotlib.pyplot as plt
import cmath
from numpy import linalg as LA


def bathy_GP(h_bathy, k_lensc, sigma, xc, x_bathy, slope_upper, slope_lower, slope):
    
    slope_1 = slope_lower; 
    slope_2 = slope_upper; 
    
    h_LBC = h_bathy[0,0];    # Depth at the left-side boundary

    x_inp_GP = x_bathy; 
    y_inp_GP = h_bathy;      
    
    x_inp_GP = np.append(x_inp_GP, [[0]], axis=0);
    y_inp_GP = np.append(y_inp_GP, [[0]], axis=0);  

    fit = np.polyfit([xc[0,0], 0], [h_LBC, 0] , 1);
    fit_fn = np.poly1d(fit) 
    
    mf = GPy.core.Mapping(1,1)
    mf.f = lambda x:fit_fn(x)
    mf.update_gradients = lambda a,b: None
    
    x_out_GP = np.zeros([1,1]);
    iter1 = 0;
    # Select all points (in x) that lie within (in the range of) the breakwater width
    for i in range(xc.shape[0]):
        if (xc[i,0] <= 0):
            iter1 = iter1 + 1;
            if (iter1==1):
                x_out_GP[0,0] = xc[i,0];

            elif (iter1>1):
                x_out_GP = np.append(x_out_GP, [[xc[i,0]]], axis=0);

    kernel2 = GPy.kern.RBF(input_dim=1, ARD=True, variance=sigma, lengthscale=k_lensc)

    gauss = GPy.likelihoods.Gaussian(variance=0.1)

    exact = GPy.inference.latent_function_inference.ExactGaussianInference()

    m1 = GPy.core.GP(X=x_inp_GP, Y=y_inp_GP, mean_function = mf, kernel=kernel2, likelihood=gauss, inference_method=exact)

    [gp_m, gp_var] = m1.predict(x_out_GP)
    
    h_bathy_only = np.zeros((xc.shape)); 
    j=0;
    for i in range(xc.shape[0]):
        if (xc[i,0] <= 0): 
            h_bathy_only[i,0] = gp_m[j,0];
            j = j + 1;
                        
    if (h_bathy_only[0,0]>-slope_1*xc[0,0]):
        h_bathy_only[0,0] = -slope_1*xc[0,0];
    
    if (h_bathy_only[0,0]<-slope_2*xc[0,0]):
        h_bathy_only[0,0] = -slope_2*xc[0,0]; 
    
    for i in range(iter1-1): 
        if (h_bathy_only[i+1,0]>=h_bathy_only[0,0]):
            h_bathy_only[i+1,0] = h_bathy_only[0,0];
        if (h_bathy_only[i+1,0]>-slope_1*xc[i+1,0]):
            h_bathy_only[i+1,0] = -slope_1*xc[i+1,0];
        if (h_bathy_only[i+1,0]<-slope_2*xc[i+1,0]):
            h_bathy_only[i+1,0] = -slope_2*xc[i+1,0];
            
    
    for i in range(xc.shape[0]):
        if (xc[i,0] <= 0):
            h_bathy_only[i,0] = h_bathy_only[i,0]/(1+np.exp(0.8*(xc[i,0]-(-4)))) +  (-slope*xc[i,0])/(1+np.exp(-0.8*(xc[i,0]-(-4))) );
            
    
    return [gp_m, h_bathy_only]           