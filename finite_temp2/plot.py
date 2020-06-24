# Plotting commands in python

import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Manual Change to make plots nicer.
plt.rcParams['figure.figsize'] = 7,5
mpl.rcParams['font.size'] = 12

# Read Files	
vel={}
eps={}
temp_list = [0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]
cmap = mpl.cm.get_cmap('coolwarm')
for temp in temp_list:
	vel['momentum'],vel['velocity'] = np.loadtxt('velon0_'+str(temp)+'.dat',unpack=True)
	plt.plot(vel['momentum'],vel['velocity'],color=cmap(temp*10),label='t = '+str(temp))
plt.xlabel(r'$k/\Lambda_0$')
plt.ylabel(r'$\dfrac{v_{\Lambda\to 0}(k)}{v_F}$')
plt.title('Temperature Dependence of Renormalised Velocity')
plt.xlim(0,1)
plt.legend()
plt.savefig("temperature_dependence_velocity.pdf")
plt.close()

for temp in temp_list:
	eps['momentum'],eps['dielectric'] = np.loadtxt('epson0_'+str(temp)+'.dat',unpack=True)
	plt.plot(eps['momentum'],eps['dielectric'],color=cmap(temp*10),label='t = '+str(temp))
plt.xlabel(r'$q/\Lambda_0$')
plt.ylabel(r'$\epsilon_{\Lambda\to 0}(q)$')
plt.title('Temperature Dependence of Renormalised Dielectric Function')
plt.xlim(0,1)
plt.ylim(1,2.5)
plt.legend()
plt.savefig("temperature_dependence_dielectric.pdf")
plt.close()
