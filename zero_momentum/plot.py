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
vel['temp'],vel['velocity'] = np.loadtxt('zero_momentum_temp_dependence.dat',unpack=True)

# vel = pd.read_csv("/home/biplab/Desktop/Code/fortran/bauer/velon0.dat",names=('momentum','velocity'),delim_whitespace=True,dtype='float64')
# eps = pd.read_csv("/home/biplab/Desktop/Code/fortran/bauer/epson0.dat",names=('momentum','dielectric'),delim_whitespace=True,dtype='float64')
def fitting_func1(x,a,b):
	return a+b*np.log(x)

def fitting_func2(x,a,b,c):
	return a+b*np.log(x) + c*np.log(x)*np.log(x)

popt1,pcov1 = curve_fit(fitting_func1,vel['temp'],vel['velocity'])
popt2,pcov2 = curve_fit(fitting_func2,vel['temp'],vel['velocity'])
y_fit1 = popt1[0] + popt1[1]*np.log(vel['temp'])
y_fit2 = popt2[0] + popt2[1]*np.log(vel['temp']) + popt2[2]*np.log(vel['temp'])*np.log(vel['temp'])

plt.plot(vel['temp'],vel['velocity'],linewidth='5',label='velocity')
plt.plot(vel['temp'],y_fit1,'black',label='a + b*log(x)\na='+str(np.round(popt1[0],2))+' b='+str(np.round(popt1[1],2)))
plt.plot(vel['temp'],y_fit2,label='a + b*log(x) + c*log(x)^2\na='+str(np.round(popt2[0],2))+' b='+str(np.round(popt2[1],2))+' c='+str(np.round(popt2[2],2)))
plt.xlabel(r'$t$')
plt.ylabel(r'$\dfrac{v_{\Lambda\to 0}(k\to0)}{v_F}$')
plt.title('Temperature Dependence of Velocity at zero Momentum')
plt.xlim(0,1)
plt.legend()
plt.savefig("zero_momentum_temp_dependence.pdf")
plt.close()


def fitting_func3(x,a,b,c,d):
	return a+b*(np.exp(-c/x) - np.exp(-d/x))

popt,pcov = curve_fit(fitting_func3,vel['temp'],vel['velocity'])
y_fit3 = popt[0] + popt[1]*(np.exp(-popt[2]/vel['temp']) - np.exp(-popt[3]/vel['temp']))
print(popt,pcov)

plt.plot(vel['temp'],y_fit3,'black',label=r'a + b$(e^{-c/t} - e^{-d/t})$'+'\na='+str(np.round(popt[0],2))+' b='+str(np.round(popt[1],2))+' c='+str(np.round(popt[2],2))+' d='+str(np.round(popt[3],2)))
plt.plot(vel['temp'],vel['velocity'],linewidth='3',label='velocity')
plt.xlabel(r'$t$')
plt.ylabel(r'$\dfrac{v_{\Lambda\to 0}(k\to0)}{v_F}$')
plt.title('Temperature Dependence of Velocity at zero Momentum')
plt.xlim(0,1)
plt.legend()
plt.savefig("zero_momentum_temp_dependence1.pdf")
plt.close()
