# -*- coding: utf-8 -*-
###############################################################################
#
#                 Particle class for the Boris method
#
###############################################################################

import numpy as np
import constants as cst

# We define the particle class
class Particle:
    def __init__(self, mass, charge):
        self.mass = mass
        self.charge = charge
        self.r = np.zeros((3), dtype='float64')
        self.v = np.zeros((3), dtype='float64')
        
    def initPos(self, x, y, z):
        self.r[:] = (x, y, z)
    
    def initSpeed(self, vx, vy, vz):
        self.v[:] = (vx, vy, vz)

    def push(self, dt, E, B):
        """
        Push the particles using Boris' method: updates position and speed of
        the particle.
        Input:  Timestep (float32)
                Electric field at the particle's position (Numpy((3), dtype='float32'))
                Magnetic field at the particle's position (Numpy((3), dtype='float32'))
        """
        # Updating speed
        
        ## Adding half of the magnetic impulse
        self.v += self.charge * E * dt / (2 * self.mass)
        
        ## Rotation of speed according to B
        t = self.charge * dt * B / (2 * self.mass)
        v1 = self.v + np.cross(self.v, t)
        self.v += np.cross(v1, 2 /(1 + np.abs(t)**2) * t)
        
        ## Adding second half of the magnetic impulse
        self.v += self.charge * E * dt / (2 * self.mass)
        
        # Updating position
        self.r += dt * self.v
        
        
class Tokamak:
    def __init__(self, R, r, n, N, I):
        self.R = R
        self.r = r
        self.n = n
        self.N = N        
        self.I = I
        
    def getB(self,x,y,z):
        r_x = x - self.R*x/np.sqrt(x**2+y**2)
        r_y = y - self.R*y/np.sqrt(x**2+y**2)
        r_z = z
        if (np.sqrt(x**2 + y**2) < (self.R-self.r)) or (np.sqrt(x**2 + y**2) > (self.R+self.r)):
            B = np.array((0, 0, 0), dtype='float64')
        else:
            if np.sqrt(r_x**2 + r_y**2 + r_z**2) > self.r:
                B = np.array((0, 0, 0), dtype='float64')
            else:
                B = np.array((0, 0, 0), dtype='float64')
                B[0] = - ((cst.mu0*self.n*self.N*self.I) / (2*np.pi*(x**2+y**2)))*y
                B[1] = ((cst.mu0*self.n*self.N*self.I) / (2*np.pi*(x**2+y**2)))*x
                B[2] = 0
        return(B)