##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np

eps = np.finfo(float).eps;
g = 9.81;
CFL = 1.5;

Amp = 0.25;       # amplitude of the crest

lb = -40;                
ub = 60;    
N = 1300;      # total number of spatial discretizations

dx = (ub-lb)/N; 
x = np.transpose(np.array([np.arange(lb, ub+0.000000000001, dx)]));  
xc = 0.5*(x[0:-1] + x[1:]);

dep_limit = 7;

slope = 0.25;
y=slope*xc; 
h0=-y;


slope_upper = 0.05;
y_upper = slope_upper * xc;
h_upper = -y_upper;

slope_lower = 0.4;
y_lower = slope_lower * xc;
h_lower = -y_lower;

interval = ((0)-(lb+0))/7;
x_bathy = np.array([[lb], [lb+interval],[lb+2*interval], [lb+3*interval], [lb+4*interval],[lb+5*interval], [lb+6*interval] ] );