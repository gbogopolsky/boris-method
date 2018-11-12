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
mpl.rcParams['legend.fontsize'] = 18

N = 1       # Number of particles
part = np.zeros((N))

# Parameters and fields
charge = cst.elemCharge
mass = cst.Mp
E = np.array((0, 0, 0))
B = np.array((0, 1, 0))
w0 = np.abs(charge) * np.sqrt(np.sum(B*B, axis=0)) / mass
dt = 0.1 / w0       # timestep
Np = 10             # Number of cyclotronic periods

Tf = Np * 2 * np.pi / w0
Nt = int(Tf // dt)  # Number of timesteps
t = np.arange(0, Nt)*dt
x, y, z = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))    # positions taken by the particle along 
vx, vy, vz = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))   # velocities

part1 = Particle(mass, charge)
part1.initPos(0, 0, 0)
part1.initSpeed(200, 0, 100)

for i in range(Nt):
    part1.push(dt, E, B)
    x[i], y[i], z[i] = part1.r
    vx[i], vy[i], vz[i] = part1.v
    
# Diagnostics
E = 0.5 * mass * (vx**2 + vy**2 + vz**2)

# Outputs
plt.plot(t,E)
plt.show()
    
    





