# -*- coding: utf-8 -*-
###############################################################################
#
#                           Test of magnetic field
#
###############################################################################

import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import constants as cst
import particle as part
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.figsize'] = (6,6)

#N = 50
#z = 0
#x, y = np.arange(N)+1, np.arange(N)+1
#Bx, By, Bz = np.zeros((N,N)), np.zeros((N,N)),  np.zeros((N,N)) 

tok = part.Tokamak(4,2,1000,20,10000)
#for i in range(N):
#    Bx[i], By[i], Bz[i] = tok.getB(x[i], y[i], z)

#plt.contourf(x, y, Bx, 64)
#plt.colorbar()
#plt.show()
#plt.contourf(x, y, By, 64)
#plt.colorbar()
#plt.show()

N = 1       # Number of particles

# Parameters and fields
charge = cst.elemCharge
mass = cst.Mp
E0 = np.array((0, 0, 0))
B0 = np.array((0, 0, 1))
w0 = np.abs(charge) * np.sqrt(np.sum(B0*B0, axis=0)) / mass
dt = 0.01 / w0       # timestep
Np = 300             # Number of cyclotronic periods

Tf = Np * 2 * np.pi / w0
Nt = int(Tf // dt)  # Number of timesteps
t = np.arange(0, Nt)*dt
x, y, z = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))    # positions taken by the particle along 
vx, vy, vz = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))   # velocities

part1 = part.Particle(mass, charge)
part1.initPos(tok.R, 0, 0)
part1.initSpeed(1e7, 1e7, 1)

for i in range(Nt):
    E = E0
    B = tok.getB(*part1.r)
    part1.push(dt, E, B)
    x[i], y[i], z[i] = part1.r
    vx[i], vy[i], vz[i] = part1.v
    

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
plt.plot(t,E)
plt.xlim((t[0], t[-1]))
plt.show()
