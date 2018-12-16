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
plt.rcParams["figure.figsize"] = (9,8)
plt.rcParams['font.size'] = 18

N = 3       # Number of particles

# Parameters and fields
charge = cst.elemCharge
mass = cst.Mp
E0 = np.array((0, 25, 0))
B0 = np.array((0, 0, 1))
w0 = np.abs(charge) * np.sqrt(np.sum(B0*B0, axis=0)) / mass
dt = 0.01 / w0       # timestep
Np = 10             # Number of cyclotronic periods

Tf = Np * 2 * np.pi / w0
Nt = int(Tf // dt)  # Number of timesteps
t = np.arange(0, Nt)*dt
x, y, z = np.zeros((N,Nt)), np.zeros((N,Nt)), np.zeros((N,Nt))    # positions taken by the particle along
vx, vy, vz = np.zeros((N,Nt)), np.zeros((N,Nt)), np.zeros((N,Nt))   # velocities

type = np.dtype(Particle)
part = np.zeros((N), dtype=type)
for n in range(N):
    part[n] = Particle(mass, charge)
    part[n].initPos(0, 0, 0)
    part[n].initSpeed(100 + n*100, 0, 0)

    for i in range(Nt):
        part[n].push(dt, E0, B0)
        x[n,i], y[n,i], z[n,i] = part[n].r
        vx[n,i], vy[n,i], vz[n,i] = part[n].v

# Diagnostics
E = np.zeros((N,Nt))
for n in range(N):
    E[n,:] = 0.5 * mass * (vx[n,:]**2 + vy[n,:]**2 + vz[n,:]**2)

# Outputs

## 3D trajectory
fig = plt.figure()
ax = fig.gca(projection='3d')
for n in range(N):
    ax.plot(x[n,:]*1e6, y[n,:]*1e6, z[n,:], label='part ' + str(n))
ax.set_xlabel('$x$ [µm]')
ax.set_ylabel('$y$ [µm]')
ax.set_zlabel('$z$ [m]')
plt.legend()
# plt.savefig('ExB_3d.png')
plt.show()

## Energy vs. time
for n in range(N):
    plt.plot(t, E[n], label='part ' + str(n))
plt.xlabel('$t$ [s]')
plt.ylabel('$E_c$ [J]')
plt.xlim((t[0], t[-1]))
plt.legend()
plt.tight_layout()
# plt.savefig('ExB_E.png')
plt.show()
