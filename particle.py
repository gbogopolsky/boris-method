# -*- coding: utf-8 -*-
###############################################################################
#
#                 Particle & tokamak class for the Boris method
#
#                          L. Fuster & G. Bogopolsky
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
    def __init__(self, R, a, n, N, I_T, I_P, config):
        self.R = R
        self.a = a
        self.n = n
        self.N = N
        self.I_T = I_T
        self.I_P = I_P
        self.config = config

    def getB(self,x,y,z):
        r_x = x - self.R*x/np.sqrt(x**2+y**2)
        r_y = y - self.R*y/np.sqrt(x**2+y**2)
        r_z = z
        r = np.sqrt(r_x**2 + r_y**2 + r_z**2)

        Bt = cst.mu0*self.n*self.N*self.I_T
        Bp = cst.mu0*self.I_P/(2*np.pi*(0.4*self.a)**2)

        if self.config == 'T':
            if r > self.a:
                BT = np.array((0, 0, 0), dtype='float64')
            else:
                BT = np.array((0, 0, 0), dtype='float64')
                BT[0] = - Bt * y / (2*np.pi*(x**2+y**2))
                BT[1] =  Bt * x / (2*np.pi*(x**2+y**2))
                BT[2] = 0
            return(BT)

        elif self.config == 'T+P':
            if r > self.a:
                BT = np.array((0, 0, 0), dtype='float64')
                BP = np.array((0, 0, 0), dtype='float64')
            else:
                BT = np.array((0, 0, 0), dtype='float64')
                BT[0] = - Bt * y / (2*np.pi*(x**2+y**2))
                BT[1] =  Bt * x / (2*np.pi*(x**2+y**2))
                BT[2] = 0

                BP = np.array((0, 0, 0), dtype='float64')
                BP[0] = -Bp*z*x/np.sqrt(x**2+y**2)
                BP[1] = -Bp*z*y/np.sqrt(x**2+y**2)
                BP[2] = Bp * (x**2 + y**2 - self.R**2)
            return(BT+BP)
        return()

    def get_field_lines(self):
        """
        Returns the field lines for the magnetic field of the tokamak.
        """
        thetas = np.array([2*np.pi*i/12 for i in range(16)])
        ds = 0.004
        N = int(2*np.pi*(self.R+self.a)/ds)

        p = np.zeros((16,3,N), dtype='float64')

        p[:,0,0] = self.R + (0.95*self.a)*np.cos(thetas)
        p[:,1,0] = 0
        p[:,2,0] = (0.95*self.a)*np.sin(thetas)

        for j in range(16):
            for i in range(N-1):
                p_interm = p[j,:,i]

                p1 = self.getB(*p_interm)#/(np.sum(self.getB(*p_interm)**2)+1)
                p_interm = p[j,:,i] + (ds/2)*p1
                p2 = self.getB(*p_interm)#/(np.sum(self.getB(*p_interm)**2)+1)
                p_interm = p[j,:,i] + (ds/2)*p2
                p3 = self.getB(*p_interm)#/(np.sum(self.getB(*p_interm)**2)+1)
                p_interm = p[j,:,i] + ds*p3
                p4 = self.getB(*p_interm)#/(np.sum(self.getB(*p_interm)**2)+1)

                p[j,:,i+1] = p[j,:,i] + (ds/6)*( p1 + 2*p2 + 2*p3 + p4 )
        return(p[:,0,:], p[:,1,:], p[:,2,:])
