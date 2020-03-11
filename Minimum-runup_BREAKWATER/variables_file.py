##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np

eps = np.finfo(float).eps;
g = 9.81;
CFL = 1.5;
# Wave characteristics 
omega = 0.5;
T_per = 2*np.pi/omega;    # time period of the wave 

Amp = 0.25;       # amplitude of the wave

lb = -40;   
ub = 35;
N = 950;     # total number of spatial discretizations

dx = (ub-lb)/N; 
x = np.transpose(np.array([np.arange(lb, ub+0.000000000001, dx)])); 
xc = 0.5*(x[0:-1] + x[1:]);


slope = 0.25;
y=slope*xc; 
h0=-y;


brk_width = 5;              #Width of the breakwater
x_brk_deeper = -35;         #Deeper limit of the center of the breakwater
x_brk_shallower = -25;      #Shallower limit of the center of the breakwater

T = 4*T_per;  #49;       # total time of the simulation

x_brk = np.array([[-brk_width/2 + 0.01],[-brk_width/3], [-brk_width/6], [0],[brk_width/6], [brk_width/3], [brk_width/2 - 0.01]]);