#!/usr/bin/env python
# coding: utf-8

# # Plots
# ---

import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.figsize'] = 7,5
mpl.rcParams['font.size'] = 12

vel = {}
eps = {}
velpole = {}
epspole = {}
velpole['momentum'],velpole['velocity'] = np.loadtxt('velcon.dat',unpack=True)
epspole['momentum'],epspole['dielectric'] = np.loadtxt('epscon.dat',unpack=True)
vel['momentum'],vel['velocity'] = np.loadtxt('../coupled2/velon0.dat',unpack=True)
eps['momentum'],eps['dielectric'] = np.loadtxt('../coupled2/epson0.dat',unpack=True)

plt.plot(vel['momentum'],vel['velocity'],label='without pole')
plt.plot(velpole['momentum'],velpole['velocity'],label='with pole')
plt.xlabel(r'$k/\Lambda_0$')
plt.ylabel(r'$v_{\Lambda\to 0}(k)$')
plt.title('Renormalised Velocity')
plt.xlim(0,1)
plt.legend()
plt.savefig("renormalised_velocity.pdf")
plt.close()

plt.plot(eps['momentum'],eps['dielectric'],linewidth='2',label='without pole')
plt.plot(epspole['momentum'],epspole['dielectric'],linewidth='2',label='with pole')
plt.xlabel(r'$q/\Lambda_0$')
plt.ylabel(r'$\epsilon_{\Lambda\to 0}(q)$')
plt.title('Renormalised Dielectric Function')
plt.xlim(0,1)
plt.legend()
plt.savefig("renormalised_dielectric.pdf")
plt.close()

## Wavefunction Renormalisation Factor
wavefunc = {}
wavefunc['momentum'],wavefunc['wavefunc'] = np.loadtxt('wavefunc1.dat',unpack=True)

plt.plot(wavefunc['momentum'],1./(1-wavefunc['wavefunc']),label='wavefunction renormalisation factor')
plt.xlabel(r'$k/\Lambda_0$')
plt.ylabel(r'$Z_{\Lambda\to 0}(k)$')
plt.ylim(-20,30)
plt.title('Wavefunction Renormalisation factor')
plt.xlim(0,1)
plt.legend()
plt.savefig("wavefunction_renormalisation_factor.pdf")
plt.close
