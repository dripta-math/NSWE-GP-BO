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

from IBV_LHS_file_1 import IBV_LHS
from m_func_file import m_func
from rhs_file import rhs
from Roe_flux_file import Roe_Flux
from time_step_file import time_step
from runup_inund_calc_file import runup_inund_calc
from breakwater_interpolate_file import brk_GP 
    
import matplotlib.pyplot as plt
import time

import variables_file

def FV_2D(x1):
    
    sigma = x1[0][0]
    k_lensc = x1[0][1]
    xc_brk = x1[0][2]
    h_brk = np.array([  [x1[0][3]], [x1[0][4]], [x1[0][5]], [x1[0][6]], [x1[0][7]], [x1[0][8]], [x1[0][9]] ]);
    
    dx = variables_file.dx;  
    Amp = variables_file.Amp; 
    T_per = variables_file.T_per; 
    omega = variables_file.omega; 
    CFL = variables_file.CFL; 
    g = variables_file.g; 
    eps = variables_file.eps; 
    lb = variables_file.lb; 
    ub = variables_file.ub; 
    N = variables_file.N; 
    h0 = variables_file.h0; 
    T = variables_file.T; 
    xc = variables_file.xc;
    brk_width = variables_file.brk_width;
    x_brk = variables_file.x_brk;
    
    

     #from the breakwater_interpolate file
    [gp_brk, h_brk_only]  = brk_GP(h_brk, k_lensc, sigma, xc_brk, brk_width,  xc, x_brk);

    h = np.zeros(xc.shape);
    for i in range(xc.shape[0]):
        if (h_brk_only[i,0] != 0):
            h[i,0]  = h_brk_only[i,0];
        else:
            h[i,0]  = h0[i,0];
    

    
    W = np.zeros((N,2,1)) ;
    W[:,0,:] =  np.maximum(h, eps);

    t = 0;
    iteration = 0;
    

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)

    runup_all = np.zeros([1,1]);
    inund_all = np.zeros([1,1]);
    time_all = np.zeros([1,1]);

    while (t<T):
        
        dt = time_step(W, CFL, dx, eps, g) ;
        # Runge-kutta 
        k1 = dt * rhs(W,          t,          N, dx, eps, h, g, Amp, T_per, xc) ;
        k2 = dt * rhs(W + 0.5*k1, t + 0.5*dt, N, dx, eps, h, g, Amp, T_per, xc) ;
        k3 = dt * rhs(W + 0.5*k2, t + 0.5*dt, N, dx, eps, h, g, Amp, T_per, xc) ;               
        k4 = dt * rhs(W + k3,     t + dt,     N, dx, eps, h, g, Amp, T_per, xc) ;
        W  = W + (k1 + 2*k2 + 2*k3 + k4)/6 ;
    
    
        [runup, inund] = runup_inund_calc(W, xc, h);

        t = t + dt;
        iteration = iteration +1; 

    
        if (iteration==1):
            runup_all[0,0] = runup;
            inund_all[0,0] = inund;

        elif (iteration>1):
            runup_all = np.append(runup_all, [[runup]], axis=0);
            inund_all = np.append(inund_all, [[inund]], axis=0);


        ax1.clear()
        ax1.plot(xc, W[:,0,:]-h, 'b-')
        ax1.plot(xc, -h, 'r-')
        ax1.set_xlabel('x (m)', fontsize=13)
        ax1.set_ylabel('y (m)', fontsize=13)
        plt.pause(1e-17)
        time.sleep(0.000001)

    plt.show()
    
    runup_max = max(runup_all[:,0]);
    inund_max = max(inund_all[:,0]);
    
    return runup_max