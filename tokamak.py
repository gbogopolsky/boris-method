# -*- coding: utf-8 -*-
###############################################################################
#
#                           Test of magnetic field
#
###############################################################################

import sys

#path = "/home/lucas_fuster/Documents/boris-method/"
#sys.path.append(path)

import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import constants as cst
import particle as part
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.figsize'] = (6,6)


########################### Paramètres Tore Supra ##############################
# Grand rayon: R = 2.40 m
# Petit rayon: r = 0.72 m
# Nombre de bobines: N = 18
# Nombre de spires par bobine: n = 2028
# Courant nominal dans les bobines: I_T = 1400 A
# Courant plasma: I_P = 1.5 MA
#
# Source : https://inis.iaea.org/collection/NCLCollectionStore/_Public/18/075/18075102.pdf   p II.14
################################################################################


tok = part.Tokamak(2.40, 0.72, 2028, 18, 1400, 1.5e6, 'T+P')

N = 1       # Number of particles

# Parameters and fields
charge = cst.elemCharge
mass = cst.Mp
E0 = np.array((0, 0, 0))
B0 = np.array((0, 0, 1))
w0 = np.abs(charge) * np.sqrt(np.sum(B0*B0, axis=0)) / mass
dt = 0.01 / w0       # timestep
Np = 10          # Number of cyclotronic periods

Tf = Np * 2 * np.pi / w0
Nt = int(Tf // dt)  # Number of timesteps
t = np.arange(0, Nt)*dt
x, y, z = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))    # positions taken by the particle along 
vx, vy, vz = np.zeros((Nt)), np.zeros((Nt)), np.zeros((Nt))   # velocities

part1 = part.Particle(mass, charge)
part1.initPos(tok.R+0.1*tok.a, 0, 0)
part1.initSpeed(1e7, 1e7, 0)

for i in range(Nt):
    E = E0
    B = tok.getB(*part1.r)
    part1.push(dt, E, B)
    x[i], y[i], z[i] = part1.r
    vx[i], vy[i], vz[i] = part1.v
    
x_line, y_line, z_line = tok.get_field_lines()


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z)
for i in range(16):
    ax.plot(x_line[i,:], y_line[i,:], z_line[i,:], 'r')
# ax.plot(x_line, y_line, z_line, 'r')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
