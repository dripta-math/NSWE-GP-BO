##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np
import cmath
import variables_file
from FV_initial_file_new import FV_2D
import pickle 
import random

xc = variables_file.xc;
h0 = variables_file.h0; 
brk_width = variables_file.brk_width;
N = variables_file.N;
x_brk = variables_file.x_brk;
x_brk_deeper = variables_file.x_brk_deeper;
x_brk_shallower = variables_file.x_brk_shallower;

xc_list = xc[:,0].tolist(); 
h0_list = h0[:,0].tolist();

num_cases = 100;

sigma_initial = np.zeros((num_cases,1));
k_lensc_initial = np.zeros((num_cases,1)); 
xc_brk_initial  = np.zeros((num_cases,1));
h_brk_1_initial = np.zeros((num_cases,1));
h_brk_2_initial = np.zeros((num_cases,1));
h_brk_3_initial = np.zeros((num_cases,1));
h_brk_4_initial = np.zeros((num_cases,1));
h_brk_5_initial = np.zeros((num_cases,1));
h_brk_6_initial = np.zeros((num_cases,1));
h_brk_7_initial = np.zeros((num_cases,1));

inputs_all = np.zeros((num_cases,1,10));

for i in range(num_cases):
    sigma_initial[i,0]   = random.uniform(0.4, 4);
    k_lensc_initial[i,0] = random.uniform(2.5, 10);
    xc_brk_initial[i,0]  = random.uniform(x_brk_deeper, x_brk_shallower);
    h_brk_1_initial[i,0] = random.uniform(2, np.interp(xc_brk_initial[i,0] + x_brk[0,0],  xc_list, h0_list));
    h_brk_2_initial[i,0] = random.uniform(2, np.interp(xc_brk_initial[i,0] + x_brk[1,0],  xc_list, h0_list));
    h_brk_3_initial[i,0] = random.uniform(2, np.interp(xc_brk_initial[i,0] + x_brk[2,0],  xc_list, h0_list));
    h_brk_4_initial[i,0] = random.uniform(2, np.interp(xc_brk_initial[i,0] + x_brk[3,0],  xc_list, h0_list));
    h_brk_5_initial[i,0] = random.uniform(2, np.interp(xc_brk_initial[i,0] + x_brk[4,0],  xc_list, h0_list));
    h_brk_6_initial[i,0] = random.uniform(2, np.interp(xc_brk_initial[i,0] + x_brk[5,0],  xc_list, h0_list));
    h_brk_7_initial[i,0] = random.uniform(2, np.interp(xc_brk_initial[i,0] + x_brk[6,0],  xc_list, h0_list));
    
    inputs_all[i,:,:] = np.array([[ sigma_initial[i,0], 
                                  k_lensc_initial[i,0], 
                                  xc_brk_initial[i,0], 
                                  h_brk_1_initial[i,0], 
                                  h_brk_2_initial[i,0], 
                                  h_brk_3_initial[i,0], 
                                  h_brk_4_initial[i,0], 
                                  h_brk_5_initial[i,0], 
                                  h_brk_6_initial[i,0], 
                                  h_brk_7_initial[i,0] ]]);
    
runup_mat = np.zeros((num_cases,1));

for i in range(num_cases):
    sim_no = i;
    data_simul = inputs_all[sim_no, :, :].tolist(); 
    [runup_mat[i,0], time_all, runup_all] = FV_2D(data_simul);
    
    if(i==0):
        runup_all_list = [ runup_all ];
        time_all_list = [ time_all ];
        f_FV_input = open('FV_input_brk_variable.pickle', 'wb');
        pickle.dump([inputs_all, runup_mat, time_all_list, runup_all_list], f_FV_input);      
        f_FV_input.close()
    else:
        f_FV_output = open('FV_input_brk_variable.pickle', 'rb');
        [a_arbit, b_arbit, time_all_list, runup_all_list] = pickle.load(f_FV_output);
        time_all_list.append( time_all );
        runup_all_list.append( runup_all );
        f_FV_output.close()
        
        f_FV_input = open('FV_input_brk_variable.pickle', 'wb');
        pickle.dump([inputs_all, runup_mat, time_all_list, runup_all_list], f_FV_input);       
        f_FV_input.close()
