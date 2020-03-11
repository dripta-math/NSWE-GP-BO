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


def brk_GP(h_brk, k_lensc, sigma, xc_brk, brk_width,  xc, x_brk):
    
    

    x_inp_GP = xc_brk + x_brk;
    y_inp_GP = h_brk;      
    
    
    x_out_GP = np.zeros([1,1]);
    iter1 = 0;
    # Select all points (in x) that lie within (in the range of) the breakwater width
    for i in range(xc.shape[0]):
        if ((xc_brk - brk_width/2) <= xc[i,0] <= (xc_brk + brk_width/2)):
            iter1 = iter1 + 1;
            if (iter1==1):
                x_out_GP[0,0] = xc[i,0];
                brk_index_start = i;
            elif (iter1>1):
                x_out_GP = np.append(x_out_GP, [[xc[i,0]]], axis=0);
                brk_index_end = i;

    kernel2 = GPy.kern.RBF(input_dim=1, ARD=True, variance=sigma, lengthscale=k_lensc)

    gauss = GPy.likelihoods.Gaussian(variance=0.1)

    exact = GPy.inference.latent_function_inference.ExactGaussianInference()

    m1 = GPy.core.GP(X=x_inp_GP, Y=y_inp_GP, kernel=kernel2, likelihood=gauss, inference_method=exact)

    [gp_m, gp_var] = m1.predict(x_out_GP)
    
 
    for i in range(gp_m.shape[0]):
        if gp_m[i,0]<2:
            gp_m[i,0]=2;


    h_brk_only = np.zeros((xc.shape)); 
    j=0;
    for i in range(xc.shape[0]):
        if (brk_index_start <= i <= brk_index_end):

            h_brk_only[i,0] = gp_m[j,0];
            j = j + 1;
    
    return [gp_m, h_brk_only]        