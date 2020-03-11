##############################################################################
# Authors: Dripta Mj,    RKMVERI, Howrah -711202, India
#          Denys Dutykh, Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, 
#                        CNRS, LAMA, 73000 Chambéry, France
##############################################################################
import numpy as np
import matplotlib.pyplot as plt

t=np.transpose([np.arange(1,151)])

t1=85;
t2=55;

Amp1=4;
Amp2=2;

T_per1=4*(t1-t2); 
T_per2=4*(t1-t2); 

h_left = 0;


h_LBC = np.zeros(t.shape);

for i in range(t.shape[0]):
    
    arg1 = (2*np.pi/T_per1) * (t[i,0] - t1);
    arg2 = (2*np.pi/T_per2) * (t[i,0] - t2);
    
    h_LBC[i,0] = h_left  +  Amp1*np.power(1/np.cosh(arg1), 2) - Amp2*np.power(1/np.cosh(arg2), 2);
    
fig = plt.figure();
ax1 = fig.add_subplot(1,1,1);
ax1.clear()
ax1.plot(t, h_LBC, 'b-')

ax1.set_xlabel('time (s)', fontsize=13)
ax1.set_ylabel('elevation (m)', fontsize=13)