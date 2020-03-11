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
from FV_initial_file import FV_2D

from GPyOpt.methods import BayesianOptimization

import variables_file
xc = variables_file.xc;
h0 = variables_file.h0; 
brk_width = variables_file.brk_width;
N = variables_file.N;
x_brk = variables_file.x_brk;

#-------------------------------------------------
import pickle    
f1 = open('FV_input_brk_variable.pickle', 'rb');
[inputs_all, runup_mat, time_all_list, runup_all_list] = pickle.load(f1);

for i in range(xc.shape[0]):
    if h0[N-1-i,0] > 2:
        index_brk_lower = i;
        x_brk_shallower = xc[N-1-i,0] - brk_width/2;
        break
        
x_brk_deeper = xc[0,0] + brk_width/2;  


x_brk_deeper = variables_file.x_brk_deeper;
x_brk_shallower = variables_file.x_brk_shallower;

domain =  [{'name': 'sigma', 'type': 'continuous', 'domain':(0.4,4), 'dimensionality':1}, 
           {'name': 'lengthscale', 'type': 'continuous', 'domain':(2.5,5), 'dimensionality':1}, 
           {'name': 'x_centre_brk', 'type': 'continuous', 'domain':(x_brk_deeper,x_brk_shallower), 'dimensionality':1}, 
           {'name': 'h_brk_1', 'type': 'continuous', 'domain':(2,6), 'dimensionality':1},
           {'name': 'h_brk_2', 'type': 'continuous', 'domain':(2,6), 'dimensionality':1},
           {'name': 'h_brk_3', 'type': 'continuous', 'domain':(2,6), 'dimensionality':1},
           {'name': 'h_brk_4', 'type': 'continuous', 'domain':(2,6), 'dimensionality':1}, 
           {'name': 'h_brk_5', 'type': 'continuous', 'domain':(2,6), 'dimensionality':1}, # ];
           {'name': 'h_brk_6', 'type': 'continuous', 'domain':(2,6), 'dimensionality':1},
           {'name': 'h_brk_7', 'type': 'continuous', 'domain':(2,6), 'dimensionality':1}];

                 

xc_list = xc[:,0].tolist(); 
h0_list = h0[:,0].tolist();  


constrains =  [ {'name': 'constr_1', 'constraint': 'x[:,3] - np.interp(x[:,2] + {},  {}, {})'.format(x_brk[0,0], xc_list, h0_list) },
                {'name': 'constr_2', 'constraint': 'x[:,4] - np.interp(x[:,2] + {},  {}, {})'.format(x_brk[1,0], xc_list, h0_list) },
                {'name': 'constr_3', 'constraint': 'x[:,5] - np.interp(x[:,2] + {},  {}, {})'.format(x_brk[2,0], xc_list, h0_list) },
                {'name': 'constr_4', 'constraint': 'x[:,6] - np.interp(x[:,2] + {},  {}, {})'.format(x_brk[3,0], xc_list, h0_list) },
                {'name': 'constr_5', 'constraint': 'x[:,7] - np.interp(x[:,2] + {},  {}, {})'.format(x_brk[4,0], xc_list, h0_list) },
                {'name': 'constr_5', 'constraint': 'x[:,8] - np.interp(x[:,2] + {},  {}, {})'.format(x_brk[5,0], xc_list, h0_list) },
                {'name': 'constr_5', 'constraint': 'x[:,9] - np.interp(x[:,2] + {},  {}, {})'.format(x_brk[6,0], xc_list, h0_list) },
                {'name': 'constr_11', 'constraint': '2 - x[:,3]'},
                {'name': 'constr_12', 'constraint': '2 - x[:,4]'},
                {'name': 'constr_13', 'constraint': '2 - x[:,5]'},
                {'name': 'constr_14', 'constraint': '2 - x[:,6]'},
                {'name': 'constr_15', 'constraint': '2 - x[:,7]'},
                {'name': 'constr_16', 'constraint': '2 - x[:,8]'},
                {'name': 'constr_17', 'constraint': '2 - x[:,9]'}];



kernel = GPy.kern.Matern52(input_dim=10, ARD=True, variance=1, lengthscale=[1,1,1,1,1,1,1,1,1,1]);
myBopt = BayesianOptimization(f=FV_2D, X= inputs_all[:,0,:], Y=runup_mat, domain=domain, constraints=constrains, 
                              kernel=kernel,
                              acquisition_type ='EI', initial_design_numdata = 0, 
                              exact_feval=True,
                              verbosity=True,
                              cost_withGradients=None,
                              model_type='GP', 
                              acquisition_optimizer_type='lbfgs') 
myBopt.run_optimization(max_iter=100)
myBopt.plot_convergence()

Bopt_out = open('breakwater_variable_position.pickle', 'wb');
pickle.dump([myBopt.X, myBopt.Y, kernel.variance[0], list(kernel.lengthscale)], Bopt_out);
Bopt_out.close();

min_runup = np.minimum.accumulate(myBopt.Y).ravel()
iteration = np.arange(1,min_runup.shape[0]+1)
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.plot(iteration, min_runup, 'r-')
ax2.set_xlabel('Iteration no.', fontsize=13)
ax2.set_ylabel('minimum runup (m)', fontsize=13)
plt.show()