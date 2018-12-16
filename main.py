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

N = 1       # Number of particles

# Parameters and fields
charge = cst.elemCharge
mass = cst.Mp
E0 = np.array((0, 0, 0))
B0 = np.array((0, 0, 1))
w0 = np.abs(charge) * np.sqrt(np.sum(B0*B0, axis=0)) / mass
dt = 0.001 / w0       # timestep
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
    part.push(dt, E0, B0)
    x[i], y[i], z[i] = part.r
    vx[i], vy[i], vz[i] = part.v

# Diagnostics
E = 0.5 * mass * (vx**2 + vy**2 + vz**2)
vx_th = 200 * np.cos(w0*t)
vx_error = np.abs(vx - vx_th)
print(vx_error.max() / vx.max())
vy_th = 200 * np.sin(w0*t + np.pi)
vy_error = np.abs(vy - vy_th)
print(vy_error.max() / vy.max())

# Outputs
## Speed along x
plt.plot(t, vx, label='numeric')
plt.plot(t, vx_th, label='analytic')
plt.plot(t, vx_error, label='erreur absolue')
plt.xlabel('$t$ [s]')
plt.ylabel('$v_x$ [m/s]')
plt.legend()
plt.xlim((t[0], t[-1]))
plt.tight_layout()
# plt.savefig('test_vx.png')
plt.show()

## Speed along y
plt.plot(t, vy, label='numeric')
plt.plot(t, vy_th, label='analytic')
plt.plot(t, vy_error, label='erreur absolue')
plt.xlabel('$t$ [s]')
plt.ylabel('$v_y$ [m/s]')
plt.legend()
plt.xlim((t[0], t[-1]))
plt.tight_layout()
# plt.savefig('test_vy.png')
plt.show()

## 3D trajectory
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x*1e6, y*1e6, z)
ax.set_xlabel('$x$ [µm]')
ax.set_ylabel('$y$ [µm]')
ax.set_zlabel('$z$ [m]')
# plt.savefig('test_3d.png')
plt.show()

## Energy vs. time
plt.plot(t, E)
plt.xlabel('$t$ [s]')
plt.ylabel('$E_c$ [J]')
plt.xlim((t[0], t[-1]))
plt.tight_layout()
# plt.savefig('test_E.png')
plt.show()
