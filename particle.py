# -*- coding: utf-8 -*-
###############################################################################
#
#                 Particle class for the Boris method
#
###############################################################################

import numpy as np

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
