##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np
import cmath
import variables_file
from FV_initial_file_new import FV_2D
xc = variables_file.xc;
h0 = variables_file.h0; 

N = variables_file.N;

xc_list = xc[:,0].tolist(); 
h0_list = h0[:,0].tolist();

slope_upper = variables_file.slope_upper;   
slope_lower = variables_file.slope_lower;   
slope = variables_file.slope;
dep_limit = variables_file.dep_limit;
x_bathy = variables_file.x_bathy;


import random

num_cases = 100;

sigma_initial = np.zeros((num_cases,1));
k_lensc_initial = np.zeros((num_cases,1)); 
h_bathy_1_initial = np.zeros((num_cases,1));
h_bathy_2_initial = np.zeros((num_cases,1));
h_bathy_3_initial = np.zeros((num_cases,1));
h_bathy_4_initial = np.zeros((num_cases,1));
h_bathy_5_initial = np.zeros((num_cases,1));
h_bathy_6_initial = np.zeros((num_cases,1));
h_bathy_7_initial = np.zeros((num_cases,1));

inputs_all = np.zeros((num_cases,1,9));

for i in range(num_cases):
    sigma_initial[i,0]   = random.uniform(0.4, 2);
    k_lensc_initial[i,0] = random.uniform(5, 15);
    h_bathy_1_initial[i,0] = random.uniform(-slope_upper*x_bathy[0,0], np.min((dep_limit, -slope_lower*x_bathy[0,0])));
    h_bathy_2_initial[i,0] = random.uniform(-slope_upper*x_bathy[1,0], np.min((dep_limit, -slope_lower*x_bathy[1,0])));
    h_bathy_3_initial[i,0] = random.uniform(-slope_upper*x_bathy[2,0], np.min((dep_limit, -slope_lower*x_bathy[2,0])));
    h_bathy_4_initial[i,0] = random.uniform(-slope_upper*x_bathy[3,0], np.min((dep_limit, -slope_lower*x_bathy[3,0])));
    h_bathy_5_initial[i,0] = random.uniform(-slope_upper*x_bathy[4,0], np.min((dep_limit, -slope_lower*x_bathy[4,0])));
    h_bathy_6_initial[i,0] = random.uniform(-slope_upper*x_bathy[5,0], np.min((dep_limit, -slope_lower*x_bathy[5,0])));
    h_bathy_7_initial[i,0] = random.uniform(-slope_upper*x_bathy[6,0], np.min((dep_limit, -slope_lower*x_bathy[6,0])));
    
    inputs_all[i,:,:] = np.array([[ sigma_initial[i,0], 
                                  k_lensc_initial[i,0], 
                                  h_bathy_1_initial[i,0], 
                                  h_bathy_2_initial[i,0], 
                                  h_bathy_3_initial[i,0], 
                                  h_bathy_4_initial[i,0], 
                                  h_bathy_5_initial[i,0], 
                                  h_bathy_6_initial[i,0], 
                                  h_bathy_7_initial[i,0] ]]);
    


runup_mat = np.zeros((num_cases,1));

import pickle    

for i in range(num_cases):
    sim_no = i;
    data_simul = inputs_all[sim_no, :, :].tolist(); 
    [runup_mat[i,0], time_all, runup_all] = FV_2D(data_simul);
    
    if(i==0):
        runup_all_list = [ runup_all ];
        time_all_list = [ time_all ];
        f_FV_input = open('FV_input_bathy_max_new.pickle', 'wb');
        pickle.dump([inputs_all, runup_mat, time_all_list, runup_all_list], f_FV_input);      
        f_FV_input.close()
    else:
        f_FV_output = open('FV_input_bathy_max_new.pickle', 'rb');
        [a_arbit, b_arbit, time_all_list, runup_all_list] = pickle.load(f_FV_output);
        time_all_list.append( time_all );
        runup_all_list.append( runup_all );
        f_FV_output.close()
        
        f_FV_input = open('FV_input_bathy_max_new.pickle', 'wb');
        pickle.dump([inputs_all, runup_mat, time_all_list, runup_all_list], f_FV_input);      
        f_FV_input.close()

