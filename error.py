# -*- coding: utf-8 -*-
###############################################################################
#
#                 Error of Boris algorithm vs the timestep
#
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import constants as cst
from particle import Particle
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.figsize'] = (6,6)

def compute(dt, E0, B0):
    Np = 5             # Number of cyclotronic periods

    Tf = Np * 2 * np.pi / w0
    Nt = int(Tf // dt)  # Number of timesteps
    t = np.arange(0, Nt)*dt
    x, y, z = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))    # positions taken by the particle along 
    vx, vy, vz = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))   # velocities

    part = Particle(mass, charge)
    part.initPos(0, 0, 0)
    part.initSpeed(200, 0, 0)

    for i in range(Nt):
        part.push(dt, E0, B0)
        x[i], y[i], z[i] = part.r
        vx[i], vy[i], vz[i] = part.v
    
    # Diagnostics
    vx_th = 200 * np.cos(w0*t)
    vx_error = np.abs(vx - vx_th) / np.max(vx_th)
    return(np.average(vx_error), np.std(vx_error))
    
# Parameters and fields
charge = cst.elemCharge
mass = cst.Mp
E0 = np.array((0, 0, 0))
B0 = np.array((0, 1, 0))
w0 = np.abs(charge) * np.sqrt(np.sum(B0*B0, axis=0)) / mass

Nsteps = 50
steps = np.logspace(-4, 1, num=Nsteps)/w0
errorsAvg, errorsStd = np.zeros(Nsteps), np.zeros(Nsteps)
dt = 0.001 / w0       # timestep

for i in range(Nsteps):
    print(i, steps[i])
    errorsAvg[i], errorsStd[i] = compute(steps[i], E0, B0)
    
np.savetxt('errorsAvg.txt', errorsAvg, delimiter=' ')
np.savetxt('errorsStd.txt', errorsStd, delimiter=' ')
    
plt.loglog(steps * w0, errorsAvg, '-o')
plt.xlabel(r'$w_0 dt$')
plt.ylabel(r'$\eta$')
plt.tight_layout()
plt.show()

    


