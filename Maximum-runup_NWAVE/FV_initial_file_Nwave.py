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

from m_func_file import m_func
from rhs_file import rhs
from Roe_flux_file import Roe_Flux
from time_step_file import time_step
from runup_inund_calc_file import runup_inund_calc
from bathy_interpolate_file import bathy_GP 
   
import matplotlib.pyplot as plt
import time

import variables_file

import scipy
from scipy.signal import find_peaks


def FV_2D_Nwave(x1):
    
    sigma = x1[0][0]
    k_lensc = x1[0][1]
    h_bathy_pts = np.array([  [x1[0][2]], [x1[0][3]], [x1[0][4]], [x1[0][5]], [x1[0][6]], [x1[0][7]], [x1[0][8]] ]);  # [x1[0][7]], [x1[0][8]] ]);
    amp_ratio = x1[0][9];    # Amplitude ratio 
    t1 = x1[0][10];           # Time of crest peak 
    t2 = x1[0][11];           # Time of triugh peak   
        
    dx = variables_file.dx;  
    Amp = variables_file.Amp; 
    CFL = variables_file.CFL; 
    g = variables_file.g; 
    eps = variables_file.eps; 
    lb = variables_file.lb; 
    ub = variables_file.ub; 
    N = variables_file.N; 
    h0 = variables_file.h0; 
    xc = variables_file.xc;
    x_bathy = variables_file.x_bathy;
    slope_upper = variables_file.slope_upper;
    slope_lower = variables_file.slope_lower;
    slope = variables_file.slope;


    T = 5*(t1 - t2) + 40;
    
     #from the breakwater_interpolate file
    [gp_bathy, h_bathy_only]  = bathy_GP(h_bathy_pts, k_lensc, sigma,  xc, x_bathy, slope_upper, slope_lower, slope);

    h = np.zeros(xc.shape); 
    for i in range(xc.shape[0]):
        if (xc[i,0] <= 0):
            h[i,0]  = h_bathy_only[i,0];
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
        k1 = dt * rhs(W,          t,          N, dx, eps, h, g, xc, Amp, amp_ratio, t1, t2) ;
        k2 = dt * rhs(W + 0.5*k1, t + 0.5*dt, N, dx, eps, h, g, xc, Amp, amp_ratio, t1, t2) ;
        k3 = dt * rhs(W + 0.5*k2, t + 0.5*dt, N, dx, eps, h, g, xc, Amp, amp_ratio, t1, t2) ;               
        k4 = dt * rhs(W + k3,     t + dt,     N, dx, eps, h, g, xc, Amp, amp_ratio, t1, t2) ;
        W  = W + (k1 + 2*k2 + 2*k3 + k4)/6 ;

    
        [runup, inund] = runup_inund_calc(W, xc, h);
        
        iteration = iteration +1; 
    
        if (iteration==1):
            runup_all[0,0] = runup;
            inund_all[0,0] = inund;
            time_all[0,0] = t;
        elif (iteration>1):
            runup_all = np.append(runup_all, [[runup]], axis=0);
            inund_all = np.append(inund_all, [[inund]], axis=0);
            time_all = np.append(time_all, [[t]], axis=0);
        
        t = t + dt;
        
    
        ax1.clear()
        ax1.plot(xc, W[:,0,:]-h, 'b-')
        ax1.plot(xc, -h, 'r-')
        ax1.set_xlabel('x (m)', fontsize=13)
        ax1.set_ylabel('y (m)', fontsize=13)
        plt.pause(1e-17)
        time.sleep(0.000001)

    plt.show()
    

    index_peaks = find_peaks(runup_all[:,0], height=0);
    runup_max = runup_all[index_peaks[0][0],0];
    

    return [runup_max, time_all, runup_all]