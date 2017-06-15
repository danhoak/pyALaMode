#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A crude python implementation of Nic Smith's a la Mode
# https://github.com/nicolassmith/alm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import constants

class beamComponent:

    # In all of the following, 'pos' is the z-position of the first part of the optical element that the beam encounters (i.e. front surface)

    def Lens(self,pos,f,label):
        # a thin lens, f>0 is converging (i.e. PLCX), f<0 is diverging (i.e. PLCC)
        self.ABCD = array(((1.0,0.0),(-1.0/f,1.0)))
        self.z = pos
        self.label = label
        self.type = 'Lens'

    def Dist(self,pos,d,label):
        # propagation through free space
        self.ABCD = array(((1.0,d),(0.0,1.0)))
        self.z = pos
        self.label = label
        self.type = 'Dist'

    def Mirror(self,pos,R,label):
        # a curved mirror, R>0 is for concave surface (converging)
        self.ABCD = array(((1.0,0.0),(-2.0/R,1.0)))
        self.z = pos
        self.label = label
        self.type = 'Mirror'

    def Slab(self,pos,d,n,label):
        # a slab of dielectric material, thickness d, index of refraction n
        self.ABCD = array(((1.0,d/n),(0.0,1.0)))
        self.z = pos
        self.label = label
        self.type = 'Slab'


class beamPath:

    def __init__(self):
        self.elementList = []

    def clearElements(self):
        self.elementList = []

    def addElement(self,element):
        self.elementList.append(element)

    def sortElements(self):
        # sorts beam path elements in ascending order of z

        elementList_raw = self.elementList
        z_raw = []
        for i in range(len(elementList_raw)):
            z_raw.append(elementList_raw[i].z)

        index_new = argsort(z_raw)

        elementList_new = []
        for i in range(len(index_new)):
            elementList_new.append(elementList_raw[index_new[i]])
        
        self.elementList = elementList_new


    def printElementLabels(self):
        # prints the labels of the optical elements in the order they are stored in elementList
        for i in range(len(self.elementList)):
            print self.elementList[i].label


    def getTransferMatrix(self):
        # returns the ABCD matrix for the beam path
        # have to sort the elementList first to make sure things are in the proper order (matrix multiplication is not commutative!)

        self.sortElements()

        M = eye(2,2)  # start with the identity and multiply from there

        for i in range(len(self.elementList)):
            M = dot(self.elementList[i].ABCD,M)

        self.transferMatrix = M
        return self.transferMatrix


    def setWavelength(self,lam):
        # set the wavelength of the laser beam, in meters
        self.wavelength = lam


    def setInputBeam(self,input_waist,input_waistpos):
        # define an input beam incident to the elements in the beam path; parameters are for a beam at z=0
        # input_waist is the virtual waist position of the input beam, 1/e^2 radius that would be observed in a power measurement, in meters
        # input_waistpos is z position of the input waist, in meters

        self.input_w0 = input_waist
        self.input_z0 = input_waistpos

    def setInputBeamPointing(self,input_displacement,input_angle):
        # define an input beam incident to the elements in the beam path; parameters are for a beam at z=0

        self.input_x = input_displacement
        self.input_theta = input_angle


    def q(self,z,w0,z0):
        # define a complex beam parameter at a given location, from a waist and a waist position (both in meters)
        # q(z) = q0 + (z - z0)
        return (constants.pi*(w0**2)/self.wavelength)*1.0j + (z-z0)

    
    def unpackComplexBeamParameter(self,z,q):
        # given a complex beam parameter at a given z, returns the waist size & waist position
        # note that the q provided must be defined at z!

        w0 = sqrt(self.wavelength*q.imag/constants.pi)
        z0 = z - q.real

        return w0, z0


    def propagateGaussianBeam(self):
        # propagate the beam through the entire beam path
        # first, calculate the complex beam parameter at the start of the path, sort the beam path, get the transfer matrix
        # then apply the transfer matrix to the beam parameter and return the complex number

        self.input_q = self.q(0.0,self.input_w0,self.input_z0)

        self.getTransferMatrix()

        vector = array(((self.input_q),(1.0)))
        M = dot(self.transferMatrix,vector)
        q_new = M[0]/M[1]

        return q_new


    def traceRay(self):
        # perform ray tracing through the entire beam path (using ABCD matrix)
        # first, sort the beam path, get the transfer matrix
        # then apply the transfer matrix to the beam parameters and store the final position & angle

        self.getTransferMatrix()

        vector = array(((self.input_x),(self.input_theta)))
        M = dot(self.transferMatrix,vector)
        self.final_x = M[0]
        self.final_theta = M[1]


    def getBeamSize(self,z1,q,z2):
        # Returns the 1/e^2 radius (in power) of the beam at position z1 for the complex beam parameter q defined at z2
        # w(z)**2 = w0 * [1 + (lambda * z / (pi * w0**2))**2]  for a waist at z=0
        
        w0,z0 = self.unpackComplexBeamParameter(z2,q)

        return w0 * sqrt(1 + (self.wavelength * (z1-z0) / (constants.pi * w0**2))**2)


    def getLastElementPos(self):
        # returns the z position of the last element in the path
        # remembers to add the distance if the last element is a 'Dist'
        if (self.elementList[-1].type == 'Dist') or (self.elementList[-1].type == 'Slab'):
            return self.elementList[-1].z + self.elementList[-1].ABCD[0,1]
        else:
            return self.elementList[-1].z

    def getLastElementLabel(self):
        return self.elementList[-1].label


    def getAngleSensitivity(self):

        # returns the sensitivity of an element (QPD) at the end of the beam path to fluctuations in angle at z=0

        self.getTransferMatrix()

        self.angleSens = self.transferMatrix[0,1]

        return self.angleSens

    def getLateralSensitivity(self):

        # returns the sensitivity of an element (QPD) at the end of the beam path to fluctuations in lateral displacement at z=0
        
        self.getTransferMatrix()

        self.lateralSens = self.transferMatrix[0,0]

        return self.lateralSens



