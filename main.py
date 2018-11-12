# -*- coding: utf-8 -*-
###############################################################################
#
#                           The Boris method
#
###############################################################################

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import constants as cst
from particle import Particle
mpl.rcParams['legend.fontsize'] = 10

N = 1       # Number of particles

# Parameters and fields
charge = cst.elemCharge
mass = cst.Mp
E0 = np.array((0, 0, 0))
B = np.array((0, 1, 0))
w0 = np.abs(charge) * np.sqrt(np.sum(B*B, axis=0)) / mass
dt = 0.01 / w0       # timestep
Np = 10             # Number of cyclotronic periods

Tf = Np * 2 * np.pi / w0
Nt = int(Tf // dt)  # Number of timesteps
t = np.arange(0, Nt)*dt
x, y, z = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))    # positions taken by the particle along 
vx, vy, vz = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))   # velocities

part = Particle(mass, charge)
part.initPos(0, 0, 0)
part.initSpeed(200, 0, 0)

for i in range(Nt):
    E = E0*np.cos(w0*t[i])
    part.push(dt, E0, B)
    x[i], y[i], z[i] = part.r
    vx[i], vy[i], vz[i] = part.v
    
# Diagnostics
E = 0.5 * mass * (vx**2 + vy**2 + vz**2)
vx_th = 200 * np.cos(w0*t)
vx_error = np.abs(vx - vx_th)
print(np.max(vx_error))

# Outputs
plt.plot(t, vx, label='numeric')
plt.plot(t, vx_th, label='analytic')
plt.plot(t, vx_error, label='erreur absolue')
plt.legend()
plt.xlim((t[0], t[-1]))
plt.show()

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot(x, y, z)
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show()
plt.plot(t,E)
plt.xlim((t[0], t[-1]))
plt.show()
    
    





