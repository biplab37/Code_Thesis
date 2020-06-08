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
vel['momentum'],vel['velocity'] = np.loadtxt('velon0.dat',unpack=True)
eps['momentum'],eps['dielectric'] = np.loadtxt('epson0.dat',unpack=True)

# vel = pd.read_csv("/home/biplab/Desktop/Code/fortran/bauer/velon0.dat",names=('momentum','velocity'),delim_whitespace=True,dtype='float64')
# eps = pd.read_csv("/home/biplab/Desktop/Code/fortran/bauer/epson0.dat",names=('momentum','dielectric'),delim_whitespace=True,dtype='float64')
def fitting_func(x,a,b):
	return a+b*np.log(x)

popt,pcov = curve_fit(fitting_func,vel['momentum'],vel['velocity'])
y_fit = popt[0] + popt[1]*np.log(vel['momentum'])

plt.plot(vel['momentum'],y_fit,'black',label='a + b*log(x)\na='+str(np.round(popt[0],2))+' b='+str(np.round(popt[1],2)))
plt.plot(vel['momentum'],vel['velocity'],linewidth='3',label='velocity')
plt.xlabel(r'$k/\Lambda_0$')
plt.ylabel(r'$v_{\Lambda\to 0}(k)$')
plt.title('Renormalised Velocity')
plt.xlim(0,1)
plt.legend()
plt.savefig("renormalised_velocity.pdf")
plt.close()

plt.plot(eps['momentum'],eps['dielectric'],linewidth='2',label='dielectric')
plt.xlabel(r'$q/\Lambda_0$')
plt.ylabel(r'$\epsilon_{\Lambda\to 0}(q)$')
plt.title('Renormalised Dielectric Function')
plt.xlim(0,1)
plt.legend()
plt.savefig("renormalised_dielectric.pdf")
plt.close()