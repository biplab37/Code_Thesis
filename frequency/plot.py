# Plotting commands in python

import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Manual Change to make plots nicer.
plt.rcParams['figure.figsize'] = 9,6
mpl.rcParams['font.size'] = 13

# Read Files	
vel={}
eps={}
freq_list = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
cmap = mpl.cm.get_cmap('coolwarm')

for freq in freq_list:
	eps['momentum'],eps['dielectric'] = np.loadtxt('epson0_'+str(freq)+'.dat',unpack=True)
	plt.plot(eps['momentum'],eps['dielectric'],color=cmap(freq),label=r'$\omega =$ '+str(freq))
plt.xlabel(r'$q/\Lambda_0$')
plt.ylabel(r'$\epsilon_{\Lambda\to 0}(q)$')
plt.title('Frequency Dependence of Renormalised Dielectric Function')
plt.xlim(0,1)
plt.legend()
plt.savefig("frequency_dependence.pdf")
plt.close()