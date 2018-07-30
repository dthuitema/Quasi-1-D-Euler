print "Script for Basic Plotting Functions"
import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import csv
import sys
#format_input='PDF';
fname = raw_input("Enter unique filename for figures:")
fname1 = "_pressure.pdf"
fname2 = "_mach.pdf"
fname3 = "_res_dens.pdf"
fname_pressure = fname + fname1
fname_mach = fname + fname2
fname_dens = fname + fname3

mpl.rcParams['font.family'] = ['serif']
rc('text', usetex=True)
m_in_flow = open('flow.csv', 'rb')
reader = csv.reader(m_in_flow)
Pt = []
M = []
x = []
for row in reader:
    x.append(float(row[0]))
    Pt.append(float(row[1]))
    M.append(float(row[4]))
m_in_flow.close()
m_in_err = open('error.csv')
reader = csv.reader(m_in_err)
Res_Dens = []
iter_num = []
for row in reader:
    Res_Dens.append(float(row[1]))
	iter_num.append(int(row[0]))
m_in_err.close()

#Plotting Section
fig, ax = plt.subplots()
ax.plot(x,Pt,label='Pressure Distribution')
plt.xlabel(r'$x$')
plt.ylabel(r'$p/p_t$')
legend = ax.legend(loc='upper center', shadow=False, fontsize='medium')
plt.grid(True)
plt.savefig(fname_pressure,dpi=None, facecolor='w', edgecolor=None,
            orientation='landscape',papertype=None, format=None,
            transparent=True, bbox_inches=None, pad_inches=0.1,
            frameon=False)
plt.show()

fig, ax = plt.subplots()
ax.plot(x,M,label='Mach Number Distribution')
plt.xlabel(r'$x$')
plt.ylabel(r'$M$')
legend = ax.legend(loc='upper center', shadow=False, fontsize='medium')
plt.grid(True)
plt.savefig(fname_mach,dpi=None, facecolor='w', edgecolor=None,
            orientation='landscape',papertype=None, format=None,
            transparent=True, bbox_inches=None, pad_inches=0.1,
            frameon=False)
plt.show()

fig, ax = plt.subplots()
ax.semilogy(iter_num,Res_Dens,label='Density Residual')
plt.xlabel(r'iterations')
plt.ylabel(r'Residual ($\rho$)')
legend = ax.legend(loc='upper center', shadow=False, fontsize='medium')
plt.grid(True)
plt.savefig(fname_dens,dpi=None, facecolor='w', edgecolor=None,
            orientation='landscape',papertype=None, format=None,
            transparent=True, bbox_inches=None, pad_inches=0.1,
            frameon=False)
plt.show()
